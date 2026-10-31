from __future__ import annotations

import csv
import math
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean
from xml.etree import ElementTree as ET


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
SUPPLEMENTARY_DATA_9 = INPUTS / "Supplementary Data 9.xlsx"
RPKM_MATRIX = INPUTS / "RPKM_veter_matrix_clean.tsv"
LIBRARY_CITY = INPUTS / "library_coordinate_coverage_audit.csv"
OUT = HERE / "results"
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

EXPECTED_COUNTS = {
    "Alphatorquevirus": 40,
    "Betatorquevirus": 17,
    "Gammatorquevirus": 10,
    "Torbevirus": 22,
    "Deltapolyomavirus": 3,
    "Gemykibivirus": 3,
    "Orthopicobirnavirus": 3,
    "Cyclovirus": 2,
    "Galvornvirus": 2,
    "Iteradensovirus": 2,
    "Betapolyomavirus": 1,
    "Gemygorvirus": 1,
    "Gemyvongvirus": 1,
    "Orthorubulavirus": 1,
    "Protoparvovirus": 1,
    "Sadetorquevirus": 1,
}


def excel_column_index(reference: str) -> int:
    match = re.match(r"([A-Z]+)", reference)
    if not match:
        raise ValueError(reference)
    index = 0
    for character in match.group(1):
        index = index * 26 + ord(character) - ord("A") + 1
    return index - 1


def shared_strings(archive: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    return ["".join(item.itertext()) for item in root.findall(f"{NS}si")]


def iter_sheet_rows(path: Path, sheet: str, width: int):
    with zipfile.ZipFile(path) as archive:
        strings = shared_strings(archive)
        with archive.open(f"xl/worksheets/{sheet}") as handle:
            for _, element in ET.iterparse(handle, events=("end",)):
                if element.tag != f"{NS}row":
                    continue
                values = [""] * width
                for cell in element.findall(f"{NS}c"):
                    index = excel_column_index(cell.attrib["r"])
                    if index >= width:
                        continue
                    value_node = cell.find(f"{NS}v")
                    cell_type = cell.attrib.get("t")
                    if value_node is None:
                        value = ""
                    elif cell_type == "s":
                        value = strings[int(value_node.text)]
                    else:
                        value = value_node.text
                    values[index] = value
                yield int(element.attrib.get("r", "0")), values
                element.clear()


def human_evidence(row: dict[str, str]) -> str:
    host = row["Host"].strip()
    title = row["stitle"].strip()
    source = row["isolation source"].strip()
    if host.lower().startswith("homo sapiens"):
        return "Reference host annotated as Homo sapiens"
    lowered = title.lower()
    if "human" in lowered:
        return "Reference title explicitly annotated as human-associated"
    if "jc polyomavirus" in lowered or "jc virus" in lowered:
        return "JC polyomavirus reference"
    if row["Virus genus"].strip() == "Torbevirus" and source.lower() in {
        "endotracheal aspirate",
        "saliva",
    }:
        return "Human respiratory/oral specimen annotation"
    return ""


def parse_human_reference_assignments():
    rows = iter_sheet_rows(SUPPLEMENTARY_DATA_9, "sheet2.xml", 28)
    _, header = next(rows)
    header = [value.strip() for value in header]
    records = []
    for source_row, values in rows:
        row = {header[index]: values[index].strip() for index in range(len(header))}
        genus = row["Virus genus"]
        if genus not in EXPECTED_COUNTS:
            continue
        evidence = human_evidence(row)
        records.append(
            {
                "source_row": source_row,
                "query_id": row["qseqid"],
                "viral_genus": genus,
                "reference_accession": row["sseqid Genbank no."],
                "reference_host": row["Host"],
                "isolation_source": row["isolation source"],
                "reference_title": row["stitle"],
                "human_association_evidence": evidence,
            }
        )

    by_genus_query = defaultdict(list)
    for row in records:
        by_genus_query[(row["viral_genus"], row["query_id"])].append(row)

    selected = []
    for (genus, query_id), evidence_rows in sorted(by_genus_query.items()):
        valid = [row for row in evidence_rows if row["human_association_evidence"]]
        if valid:
            chosen = valid[0]
        else:
            continue
        selected.append(chosen)

    observed = Counter(row["viral_genus"] for row in selected)
    if dict(observed) != EXPECTED_COUNTS:
        raise AssertionError(
            "Human-reference curation did not reproduce the manuscript counts.\n"
            f"Observed: {dict(observed)}\nExpected: {EXPECTED_COUNTS}"
        )
    if len(selected) != sum(EXPECTED_COUNTS.values()):
        raise AssertionError(len(selected))
    return selected


def canonical_query_id(matrix_contig: str) -> str:
    match = re.match(
        r"^([A-Za-z]+)(?:_[A-Za-z]+)?_(?:BirdFe|Birdfe|Swabs)(\d+)_k141_(.+)$",
        matrix_contig,
    )
    if not match:
        return ""
    return f"{match.group(1).upper()}_{match.group(2)}_{match.group(3)}"


def read_library_city():
    mapping = {}
    coordinates = {}
    with LIBRARY_CITY.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            mapping[row["Library"]] = row["city"]
            coordinates[row["city"]] = (
                float(row["longitude"]),
                float(row["latitude"]),
            )
    return mapping, coordinates


def read_selected_rpkm(selected):
    selected_queries = {row["query_id"] for row in selected}
    found = {}

    with RPKM_MATRIX.open(encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        libraries = header[1:]
        for line in handle:
            first_tab = line.find("\t")
            if first_tab < 0:
                continue
            contig = line[:first_tab]
            query_id = canonical_query_id(contig)
            if query_id not in selected_queries:
                continue
            values = [float(value) for value in line[first_tab + 1 :].rstrip("\r\n").split("\t")]
            if len(values) != len(libraries):
                raise AssertionError((contig, len(values), len(libraries)))
            if query_id in found:
                raise AssertionError(f"Duplicate canonical query ID: {query_id}")
            found[query_id] = (contig, values)

    missing = sorted(selected_queries - set(found))
    if missing:
        raise AssertionError(f"Missing {len(missing)} query IDs from RPKM matrix: {missing[:10]}")
    return libraries, found


def calculate_profiles(selected, libraries, found, library_city, coordinates):
    missing_libraries = sorted(set(libraries) - set(library_city))
    if missing_libraries:
        raise AssertionError(f"Libraries missing city metadata: {missing_libraries[:10]}")
    if len(libraries) != 985:
        raise AssertionError(len(libraries))

    genera = list(EXPECTED_COUNTS)
    genus_library = {genus: [0.0] * len(libraries) for genus in genera}
    query_audit = []
    assignments_by_query = defaultdict(list)
    for row in selected:
        assignments_by_query[row["query_id"]].append(row)

    for query_id, (matrix_contig, values) in found.items():
        for assignment in assignments_by_query[query_id]:
            genus = assignment["viral_genus"]
            totals = genus_library[genus]
            for index, value in enumerate(values):
                totals[index] += value
            row = dict(assignment)
            row["matrix_contig_id"] = matrix_contig
            row["libraries_detected"] = sum(value > 0 for value in values)
            row["mean_rpkm_all_libraries"] = mean(values)
            row["nonexclusive_genus_assignment"] = (
                "Yes" if len(assignments_by_query[query_id]) > 1 else "No"
            )
            query_audit.append(row)

    city_libraries = defaultdict(list)
    for index, library in enumerate(libraries):
        city_libraries[library_city[library]].append(index)
    if len(city_libraries) != 22:
        raise AssertionError(len(city_libraries))

    city_rows = []
    nationwide_rows = []
    for genus in genera:
        values = genus_library[genus]
        detected = [value > 0 for value in values]
        positive = [value for value in values if value > 0]
        city_detected_count = 0
        for city, indices in sorted(city_libraries.items(), key=lambda item: coordinates[item[0]][0]):
            city_values = [values[index] for index in indices]
            city_positive = [value for value in city_values if value > 0]
            detected_libraries = len(city_positive)
            if detected_libraries:
                city_detected_count += 1
            longitude, latitude = coordinates[city]
            city_rows.append(
                {
                    "viral_genus": genus,
                    "city": city,
                    "longitude": longitude,
                    "latitude": latitude,
                    "libraries_sampled": len(indices),
                    "libraries_detected": detected_libraries,
                    "prevalence": detected_libraries / len(indices),
                    "mean_rpkm_all_libraries": mean(city_values),
                    "mean_rpkm_positive_libraries": mean(city_positive) if city_positive else 0.0,
                }
            )
        nationwide_rows.append(
            {
                "viral_genus": genus,
                "query_sequence_count": EXPECTED_COUNTS[genus],
                "libraries_sampled": len(libraries),
                "libraries_detected": sum(detected),
                "nationwide_prevalence": sum(detected) / len(libraries),
                "cities_sampled": len(city_libraries),
                "cities_detected": city_detected_count,
                "mean_rpkm_all_libraries": mean(values),
                "mean_rpkm_positive_libraries": mean(positive) if positive else 0.0,
                "maximum_library_rpkm": max(values),
            }
        )

    nationwide_rows.sort(
        key=lambda row: (
            -row["nationwide_prevalence"],
            -row["mean_rpkm_all_libraries"],
            row["viral_genus"],
        )
    )
    order = {row["viral_genus"]: index for index, row in enumerate(nationwide_rows)}
    city_rows.sort(key=lambda row: (order[row["viral_genus"]], row["longitude"]))
    query_audit.sort(key=lambda row: (order[row["viral_genus"]], row["query_id"]))
    return nationwide_rows, city_rows, query_audit, genus_library


def write_csv(path: Path, rows):
    if not rows:
        raise ValueError(path)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    selected = parse_human_reference_assignments()
    library_city, coordinates = read_library_city()
    libraries, found = read_selected_rpkm(selected)
    nationwide, city_rows, query_audit, genus_library = calculate_profiles(
        selected, libraries, found, library_city, coordinates
    )

    write_csv(OUT / "genus_summary.csv", nationwide)
    write_csv(OUT / "city_genus_profiles.csv", city_rows)
    write_csv(OUT / "query_assignments.csv", query_audit)

    library_rows = []
    for genus, values in genus_library.items():
        for library, value in zip(libraries, values):
            library_rows.append(
                {
                    "viral_genus": genus,
                    "library": library,
                    "city": library_city[library],
                    "rpkm": value,
                    "detected": 1 if value > 0 else 0,
                }
            )
    write_csv(OUT / "library_genus_rpkm.csv", library_rows)


if __name__ == "__main__":
    main()
