from __future__ import annotations

import csv
import math
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET


OUTDIR = Path(__file__).resolve().parent
SOURCE = OUTDIR / "inputs" / "Supplementary Data 8.xlsx"
RESULTS = OUTDIR / "results"
DNA_GROUPS = {
    "Adenoviridae",
    "Anelloviridae",
    "Cressdnaviricota",
    "Papillomaviridae",
    "Parvoviridae",
    "Polyomaviridae",
}


def as_number(value):
    if value is None:
        return None
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if math.isfinite(float(value)):
            return float(value)
        return None
    text = str(value).strip()
    if text in {"", "/", "-", "NA", "N/A", "nan", "None"}:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def clean(value):
    return "" if value is None else str(value).strip()


def excel_column_index(reference):
    match = re.match(r"([A-Z]+)", reference)
    if not match:
        raise ValueError(f"Invalid cell reference: {reference}")
    value = 0
    for character in match.group(1):
        value = value * 26 + ord(character) - ord("A") + 1
    return value - 1


def iter_sheet_rows(path):
    namespace = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    with zipfile.ZipFile(path) as archive:
        shared_strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared_strings = ["".join(item.itertext()) for item in root.findall(f"{namespace}si")]
        with archive.open("xl/worksheets/sheet1.xml") as sheet:
            for _, element in ET.iterparse(sheet, events=("end",)):
                if element.tag != f"{namespace}row":
                    continue
                values = [None] * 17
                for cell in element.findall(f"{namespace}c"):
                    index = excel_column_index(cell.attrib["r"])
                    if index >= len(values):
                        continue
                    cell_type = cell.attrib.get("t")
                    value_node = cell.find(f"{namespace}v")
                    if cell_type == "inlineStr":
                        inline = cell.find(f"{namespace}is")
                        value = "" if inline is None else "".join(inline.itertext())
                    elif value_node is None:
                        value = None
                    elif cell_type == "s":
                        value = shared_strings[int(value_node.text)]
                    else:
                        value = value_node.text
                    values[index] = value
                yield int(element.attrib.get("r", "0")), values
                element.clear()


def parse_source():
    family = None
    method = None
    header = None
    records = defaultdict(list)
    invalid_rows = []

    for row_number, row in iter_sheet_rows(SOURCE):
        first = row[0]
        nonempty = sum(value is not None for value in row)

        if nonempty == 1 and isinstance(first, str) and "BLAST" in first:
            match = re.match(r"^(.+?)\s*\((BLAST[^)]*)\)\s*$", first.strip())
            if not match:
                raise ValueError(f"Unrecognized family header on row {row_number}: {first}")
            family, method = match.group(1), match.group(2)
            header = None
            continue

        if family and isinstance(first, str) and first.strip().lower() == "qseqid":
            header = [clean(value).lower() for value in row]
            continue

        if not family or not header or first in (None, ""):
            continue

        item = {
            header[index]: row[index]
            for index in range(min(len(header), len(row)))
            if header[index]
        }
        query_id = clean(item.get("qseqid"))
        if not query_id:
            continue

        pident = as_number(item.get("pident"))
        bitscore = as_number(item.get("bitscore"))
        hit = {
            "source_row": row_number,
            "family": family,
            "alignment_method": method,
            "query_id": query_id,
            "reference_accession": clean(item.get("sseqid")),
            "reference_title": clean(item.get("stitle")),
            "percent_identity": pident,
            "bitscore": bitscore,
            "alignment_length": as_number(item.get("length")),
            "query_coverage_hsp": as_number(item.get("qcovhsp")),
            "subject_coverage_hsp": as_number(item.get("scovhsp")),
        }
        records[(family, query_id)].append(hit)
        if pident is None:
            invalid_rows.append(hit)

    return records, invalid_rows


def choose_hits(records):
    selected = []
    comparison = []
    for (family, query_id), hits in records.items():
        valid_identity = [hit for hit in hits if hit["percent_identity"] is not None]
        valid_bitscore = [hit for hit in valid_identity if hit["bitscore"] is not None]

        first_valid = valid_identity[0] if valid_identity else None
        max_identity = (
            max(
                valid_identity,
                key=lambda hit: (
                    hit["percent_identity"],
                    hit["bitscore"] if hit["bitscore"] is not None else -math.inf,
                ),
            )
            if valid_identity
            else None
        )
        top_bitscore = (
            max(
                valid_bitscore,
                key=lambda hit: (hit["bitscore"], hit["percent_identity"]),
            )
            if valid_bitscore
            else first_valid
        )

                                                                             
                                                                             
                                                                              
                                                               
        chosen = max_identity
        if chosen is None:
            selected.append(
                {
                    "query_id": query_id,
                    "viral_group": family,
                    "genome_type": "DNA" if family in DNA_GROUPS else "RNA",
                    "alignment_method": hits[0]["alignment_method"],
                    "closest_reference_accession": "",
                    "closest_reference_title": hits[0]["reference_accession"],
                    "pairwise_identity_percent": "",
                    "bitscore": "",
                    "alignment_length": "",
                    "query_coverage_percent": "",
                    "subject_coverage_percent": "",
                    "below_80_percent": "",
                    "identity_available": "No",
                    "selection_rule": "Maximum percent identity among retained matches",
                    "n_reported_hits": len(hits),
                    "source_row": hits[0]["source_row"],
                }
            )
        else:
            identity = chosen["percent_identity"]
            selected.append(
                {
                    "query_id": query_id,
                    "viral_group": family,
                    "genome_type": "DNA" if family in DNA_GROUPS else "RNA",
                    "alignment_method": chosen["alignment_method"],
                    "closest_reference_accession": chosen["reference_accession"],
                    "closest_reference_title": chosen["reference_title"],
                    "pairwise_identity_percent": identity,
                    "bitscore": chosen["bitscore"] if chosen["bitscore"] is not None else "",
                    "alignment_length": chosen["alignment_length"] if chosen["alignment_length"] is not None else "",
                    "query_coverage_percent": chosen["query_coverage_hsp"] if chosen["query_coverage_hsp"] is not None else "",
                    "subject_coverage_percent": chosen["subject_coverage_hsp"] if chosen["subject_coverage_hsp"] is not None else "",
                    "below_80_percent": "Yes" if identity < 80 else "No",
                    "identity_available": "Yes",
                    "selection_rule": "Maximum percent identity among retained matches",
                    "n_reported_hits": len(hits),
                    "source_row": chosen["source_row"],
                }
            )

        comparison.append(
            {
                "query_id": query_id,
                "viral_group": family,
                "genome_type": "DNA" if family in DNA_GROUPS else "RNA",
                "n_reported_hits": len(hits),
                "first_valid_identity": "" if first_valid is None else first_valid["percent_identity"],
                "top_bitscore_identity": "" if top_bitscore is None else top_bitscore["percent_identity"],
                "top_bitscore_query_coverage_percent": "" if top_bitscore is None else top_bitscore["query_coverage_hsp"],
                "maximum_identity": "" if max_identity is None else max_identity["percent_identity"],
                "maximum_identity_query_coverage_percent": "" if max_identity is None else max_identity["query_coverage_hsp"],
                "maximum_identity_alignment_length": "" if max_identity is None else max_identity["alignment_length"],
            }
        )
    return selected, comparison


def summarize(selected):
    def type7_quantile(sorted_values, probability):
        if not sorted_values:
            return math.nan
        position = (len(sorted_values) - 1) * probability
        lower = math.floor(position)
        upper = math.ceil(position)
        if lower == upper:
            return sorted_values[lower]
        fraction = position - lower
        return sorted_values[lower] * (1 - fraction) + sorted_values[upper] * fraction

    by_family = defaultdict(list)
    for row in selected:
        by_family[row["viral_group"]].append(row)

    summary = []
    for family, rows in by_family.items():
        values = [float(row["pairwise_identity_percent"]) for row in rows if row["identity_available"] == "Yes"]
        values.sort()
        below = sum(value < 80 for value in values)
        n = len(values)
        q1 = type7_quantile(values, 0.25)
        median = type7_quantile(values, 0.50)
        q3 = type7_quantile(values, 0.75)
        summary.append(
            {
                "viral_group": family,
                "genome_type": rows[0]["genome_type"],
                "alignment_method": rows[0]["alignment_method"],
                "n_sequences_total": len(rows),
                "n_with_identity": n,
                "n_missing_identity": len(rows) - n,
                "n_below_80_percent": below,
                "percent_below_80_of_total": 100 * below / len(rows),
                "percent_below_80_of_numeric": 100 * below / n if n else math.nan,
                "median_identity_percent": median,
                "q1_identity_percent": q1,
                "q3_identity_percent": q3,
                "minimum_identity_percent": min(values) if values else math.nan,
                "maximum_identity_percent": max(values) if values else math.nan,
            }
        )
    return sorted(summary, key=lambda row: (row["median_identity_percent"], row["viral_group"]))


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    records, invalid_rows = parse_source()
    selected, comparison = choose_hits(records)
    summary = summarize(selected)

    write_csv(RESULTS / "sequence_identity_by_family.csv", selected)
    write_csv(RESULTS / "family_sequence_identity_summary.csv", summary)
    write_csv(RESULTS / "hit_selection_comparison.csv", comparison)
    write_csv(RESULTS / "rows_with_non_numeric_identity.csv", invalid_rows)


if __name__ == "__main__":
    main()
