from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
RESULTS = HERE / "results"

SUPP1 = INPUTS / "Supplementary Data 1.xlsx"
SUPP11 = INPUTS / "Supplementary Data 11.xlsx"
RPKM = INPUTS / "RPKM_veter_matrix_clean.tsv"
NETWORK = INPUTS / "virus_network_results_fig5b1.xlsx"
COORDINATES = INPUTS / "city_coordinates.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def find_unique(pattern: str, root: Path) -> Path:
    hits = list(root.rglob(pattern))
    if len(hits) != 1:
        raise RuntimeError(f"Expected one {pattern!r} below {root}, found {len(hits)}: {hits}")
    return hits[0]


def normalize_matrix_id(raw_id: str) -> str | None:
    match = re.match(r"^([A-Z]+).*?(\d+)_k141_(\d+)$", str(raw_id).strip())
    if not match:
        return None
    return f"{match.group(1)}_{match.group(2)}_{match.group(3)}"


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    out = np.empty(n, dtype=float)
    out[order] = adjusted
    return out


def haversine_matrix(longitude: np.ndarray, latitude: np.ndarray) -> np.ndarray:
    lon = np.radians(longitude)
    lat = np.radians(latitude)
    dlon = lon[:, None] - lon[None, :]
    dlat = lat[:, None] - lat[None, :]
    a = np.sin(dlat / 2) ** 2 + np.cos(lat[:, None]) * np.cos(lat[None, :]) * np.sin(dlon / 2) ** 2
    return 6371.0088 * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def symmetric_knn_weights(distance: np.ndarray, k: int) -> np.ndarray:
    n = distance.shape[0]
    weights = np.zeros((n, n), dtype=float)
    for i in range(n):
        neighbours = np.argsort(distance[i])[1 : k + 1]
        weights[i, neighbours] = 1.0
    weights = np.maximum(weights, weights.T)
    np.fill_diagonal(weights, 0.0)
    return weights


def global_moran(x: np.ndarray, weights: np.ndarray) -> float:
    z = x - np.mean(x)
    denominator = float(np.dot(z, z))
    if denominator == 0:
        return float("nan")
    s0 = float(weights.sum())
    return float(len(x) / s0 * np.sum(weights * (z[:, None] * z[None, :])) / denominator)


def global_moran_test(x: np.ndarray, weights: np.ndarray, permutations: int, rng: np.random.Generator):
    observed = global_moran(x, weights)
    expected = -1.0 / (len(x) - 1)
    permuted = np.empty(permutations, dtype=float)
    for b in range(permutations):
        permuted[b] = global_moran(rng.permutation(x), weights)
    p_two_sided = (1 + np.count_nonzero(np.abs(permuted - expected) >= abs(observed - expected))) / (
        permutations + 1
    )
    return observed, expected, float(p_two_sided)


def local_gistar(x: np.ndarray, weights_with_self: np.ndarray) -> np.ndarray:
    n = len(x)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if sd == 0:
        return np.zeros(n, dtype=float)
    w_sum = weights_with_self.sum(axis=1)
    w_sq_sum = np.square(weights_with_self).sum(axis=1)
    denominator = sd * np.sqrt((n * w_sq_sum - np.square(w_sum)) / (n - 1))
    numerator = weights_with_self @ x - mean * w_sum
    return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)


def local_gistar_test(
    x: np.ndarray,
    weights_with_self: np.ndarray,
    permutations: int,
    rng: np.random.Generator,
):
    observed = local_gistar(x, weights_with_self)
    greater = np.zeros(len(x), dtype=int)
    less = np.zeros(len(x), dtype=int)
    for _ in range(permutations):
        permuted = local_gistar(rng.permutation(x), weights_with_self)
        greater += permuted >= observed
        less += permuted <= observed
    p_hot = (greater + 1) / (permutations + 1)
    p_cold = (less + 1) / (permutations + 1)
    p_two_sided = np.minimum(1.0, 2.0 * np.minimum(p_hot, p_cold))
    q_two_sided = bh_adjust(p_two_sided)
    classification = np.where(
        (observed > 0) & (q_two_sided < 0.05),
        "Hot spot",
        np.where((observed < 0) & (q_two_sided < 0.05), "Cold spot", "Not significant"),
    )
    return observed, p_hot, p_cold, p_two_sided, q_two_sided, classification


def load_target_definitions():
    network_path = NETWORK
    host_count = pd.read_excel(network_path, sheet_name="virus_host_count")
    host_count["virus"] = host_count["virus"].astype(str).str.strip()
    multi_host = set(host_count.loc[host_count["host_n"] >= 2, "virus"])

    risk = pd.read_excel(SUPP11, sheet_name=0, header=1)
    risk["Name"] = risk["Name"].astype(str).str.strip()
    risk["priority_category"] = risk["priority_category"].astype(str).str.strip()
    high_risk = set(risk.loc[risk["priority_category"].isin(["Very high", "High"]), "Name"])
    dual_priority = multi_host & high_risk

    if len(multi_host) != 24:
        raise AssertionError(f"Expected 24 multi-host sequences, found {len(multi_host)}")
    if len(dual_priority) != 18:
        raise AssertionError(f"Expected 18 multi-host Very high/High sequences, found {len(dual_priority)}")
    if risk["Name"].nunique() != 4746:
        raise AssertionError(f"Expected 4,746 risk-scored sequences, found {risk['Name'].nunique()}")

    definitions = risk[["Name", "taxonomy", "calibrated_score_mean", "calibrated_score_lower", "calibrated_score_upper", "bagged_prediction", "priority_category"]].copy()
    definitions["multi_host_sharing"] = definitions["Name"].isin(multi_host)
    definitions["high_zoonotic_risk"] = definitions["Name"].isin(high_risk)
    definitions["dual_priority"] = definitions["Name"].isin(dual_priority)
    definitions = definitions.rename(columns={"Name": "qseqid"})
    return definitions, multi_host, high_risk, dual_priority, network_path


def load_target_rpkm(target_ids: set[str]):
    found: dict[str, np.ndarray] = {}
    raw_ids: dict[str, str] = {}
    duplicate_normalized: set[str] = set()
    all_normalized: set[str] = set()
    total_rows = 0
    with RPKM.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        libraries = header[1:]
        for row in reader:
            total_rows += 1
            normalized = normalize_matrix_id(row[0])
            if normalized is None:
                continue
            if normalized in all_normalized:
                duplicate_normalized.add(normalized)
            all_normalized.add(normalized)
            if normalized in target_ids:
                found[normalized] = np.asarray(row[1:], dtype=float)
                raw_ids[normalized] = row[0]
    if total_rows != 17439 or len(libraries) != 985 or len(set(libraries)) != 985:
        raise AssertionError("Unexpected RPKM dimensions")
    if duplicate_normalized:
        raise AssertionError(f"Duplicate normalized RPKM IDs: {sorted(duplicate_normalized)[:5]}")
    missing = target_ids - found.keys()
    return libraries, found, raw_ids, total_rows, missing


def load_metadata(libraries: list[str]):
    metadata = pd.read_excel(SUPP1, sheet_name=0, header=1)
    metadata["library_key"] = metadata["Library"].astype(str).str.casefold()
    metadata["city"] = (
        metadata["Location of Sampling (Sampling site, Province)"]
        .astype(str)
        .str.split(",", n=1)
        .str[0]
        .str.strip()
    )
    metadata["host_order"] = metadata["Inferred host order"].fillna("Unassigned").astype(str).str.strip()
    if len(metadata) != 985 or metadata["library_key"].nunique() != 985:
        raise AssertionError("Unexpected metadata dimensions or duplicate library IDs")
    lookup = metadata.set_index("library_key")
    missing = [library for library in libraries if library.casefold() not in lookup.index]
    if missing:
        raise AssertionError(f"Unmapped libraries: {missing[:10]}")
    ordered = lookup.loc[[library.casefold() for library in libraries]].reset_index(drop=True)
    if ordered["city"].nunique() != 22:
        raise AssertionError("Expected 22 cities")
    return ordered


def load_coordinates(cities: list[str]):
    coordinate_path = COORDINATES
    coordinates = pd.read_csv(coordinate_path)
    city_col = next(column for column in coordinates.columns if column.lower() == "city")
    lon_col = next(column for column in coordinates.columns if "longitude" in column.lower())
    lat_col = next(column for column in coordinates.columns if "latitude" in column.lower())
    coordinates = coordinates[[city_col, lon_col, lat_col]].rename(
        columns={city_col: "city", lon_col: "longitude", lat_col: "latitude"}
    )
    coordinates["city"] = coordinates["city"].astype(str).str.strip()
    coordinates = coordinates.drop_duplicates("city")
    missing = set(cities) - set(coordinates["city"])
    if missing:
        raise AssertionError(f"Missing coordinates for: {sorted(missing)}")
    return coordinates.set_index("city").loc[cities].reset_index(), coordinate_path


def group_matrix(ids: set[str], rpkm_rows: dict[str, np.ndarray]) -> np.ndarray:
    return np.vstack([rpkm_rows[qseqid] for qseqid in sorted(ids)])


def expected_incidence_richness(detection_counts: np.ndarray, n_libraries: int, target_n: int = 33) -> float:
                                                                                        
    if target_n > n_libraries:
        raise ValueError("target_n cannot exceed the city library count")
    log_total = math.lgamma(n_libraries + 1) - math.lgamma(target_n + 1) - math.lgamma(n_libraries - target_n + 1)
    expected = 0.0
    for detected_n in np.asarray(detection_counts, dtype=int):
        if detected_n <= 0:
            continue
        absent_n = n_libraries - detected_n
        if absent_n < target_n:
            expected += 1.0
        else:
            log_absent = math.lgamma(absent_n + 1) - math.lgamma(target_n + 1) - math.lgamma(absent_n - target_n + 1)
            expected += 1.0 - math.exp(log_absent - log_total)
    return float(expected)


def summarize_cities(
    libraries: list[str],
    metadata: pd.DataFrame,
    definitions: pd.DataFrame,
    multi_host: set[str],
    high_risk: set[str],
    dual_priority: set[str],
    rpkm_rows: dict[str, np.ndarray],
    raw_ids: dict[str, str],
    coordinates: pd.DataFrame,
):
    groups = {
        "multi_host": multi_host,
        "high_risk": high_risk,
        "dual_priority": dual_priority,
    }
    matrices = {name: group_matrix(ids, rpkm_rows) for name, ids in groups.items()}

    library_summary = pd.DataFrame(
        {
            "library": libraries,
            "city": metadata["city"].values,
            "host_order": metadata["host_order"].values,
        }
    )
    for name, matrix in matrices.items():
        library_summary[f"{name}_sequence_richness"] = np.count_nonzero(matrix > 0, axis=0)
        library_summary[f"{name}_total_rpkm"] = np.sum(matrix, axis=0)

    city_rows = []
    for city, group in library_summary.groupby("city", sort=True):
        row: dict[str, float | int | str] = {"city": city, "n_libraries": len(group)}
        city_mask = metadata["city"].to_numpy() == city
        for name, matrix in matrices.items():
            city_matrix = matrix[:, city_mask]
            detected_sequences = np.any(city_matrix > 0, axis=1)
            sequence_detection_counts = np.count_nonzero(city_matrix > 0, axis=1)
            per_library_richness = np.count_nonzero(city_matrix > 0, axis=0)
            per_library_abundance = city_matrix.sum(axis=0)
            positive = per_library_richness > 0
            row[f"{name}_sequence_richness"] = int(detected_sequences.sum())
            row[f"{name}_standardized_richness_33"] = expected_incidence_richness(
                sequence_detection_counts, n_libraries=len(group), target_n=33
            )
            row[f"{name}_positive_libraries"] = int(positive.sum())
            row[f"{name}_prevalence_percent"] = float(100 * positive.mean())
            row[f"{name}_mean_sequences_per_library"] = float(per_library_richness.mean())
            row[f"{name}_mean_rpkm_all_libraries"] = float(per_library_abundance.mean())
            row[f"{name}_median_rpkm_all_libraries"] = float(np.median(per_library_abundance))
        city_rows.append(row)
    city_summary = pd.DataFrame(city_rows).merge(coordinates, on="city", how="left", validate="one_to_one")

                                                                     
                                                                             
                                                                              
                                                                            
    rank_inputs = [
        "multi_host_standardized_richness_33",
        "high_risk_standardized_richness_33",
        "dual_priority_standardized_richness_33",
    ]
    percentile_ranks = city_summary[rank_inputs].rank(pct=True, method="average")
    city_summary["surveillance_priority_score"] = 100 * percentile_ranks.mean(axis=1)
    city_summary["surveillance_priority_rank"] = city_summary["surveillance_priority_score"].rank(
        ascending=False, method="min"
    ).astype(int)

    sequence_rows = []
    definition_lookup = definitions.set_index("qseqid")
    for qseqid in sorted(multi_host | high_risk):
        values = rpkm_rows[qseqid]
        detected = values > 0
        cities_detected = sorted(set(metadata.loc[detected, "city"]))
        host_orders_detected = sorted(set(metadata.loc[detected, "host_order"]))
        info = definition_lookup.loc[qseqid]
        sequence_rows.append(
            {
                "qseqid": qseqid,
                "rpkm_matrix_row_id": raw_ids[qseqid],
                "taxonomy": info["taxonomy"],
                "priority_category": info["priority_category"],
                "calibrated_score_mean": info["calibrated_score_mean"],
                "multi_host_sharing": qseqid in multi_host,
                "high_zoonotic_risk": qseqid in high_risk,
                "dual_priority": qseqid in dual_priority,
                "n_positive_libraries": int(detected.sum()),
                "n_cities_detected": len(cities_detected),
                "cities_detected": "; ".join(cities_detected),
                "n_study_host_orders_detected": len(host_orders_detected),
                "study_host_orders_detected": "; ".join(host_orders_detected),
                "mean_rpkm_all_libraries": float(values.mean()),
                "maximum_rpkm": float(values.max()),
            }
        )
    sequence_summary = pd.DataFrame(sequence_rows)
    return city_summary, library_summary, sequence_summary


def spatial_statistics(city_summary: pd.DataFrame, k: int, permutations: int, seed: int):
    cities = city_summary["city"].tolist()
    distance = haversine_matrix(city_summary["longitude"].to_numpy(), city_summary["latitude"].to_numpy())
    weights = symmetric_knn_weights(distance, k=k)
    weights_self = weights.copy()
    np.fill_diagonal(weights_self, 1.0)

    edge_rows = []
    for i in range(len(cities)):
        for j in range(i + 1, len(cities)):
            if weights[i, j] > 0:
                edge_rows.append(
                    {"city_1": cities[i], "city_2": cities[j], "distance_km": float(distance[i, j]), "weight": 1}
                )
    edges = pd.DataFrame(edge_rows)

    outcomes = [
        "multi_host_standardized_richness_33",
        "high_risk_standardized_richness_33",
        "dual_priority_standardized_richness_33",
    ]
    global_rows = []
    local_rows = []
    rng = np.random.default_rng(seed)
    for outcome in outcomes:
        x = city_summary[outcome].to_numpy(dtype=float)
        moran, expected, moran_p = global_moran_test(x, weights, permutations, rng)
        global_rows.append(
            {
                "outcome": outcome,
                "n_cities": len(cities),
                "knn_k": k,
                "permutations": permutations,
                "moran_i": moran,
                "expected_moran_i": expected,
                "permutation_p_two_sided": moran_p,
            }
        )
        gi, p_hot, p_cold, p_two, q_two, classification = local_gistar_test(
            x, weights_self, permutations, rng
        )
        for index, city in enumerate(cities):
            local_rows.append(
                {
                    "outcome": outcome,
                    "city": city,
                    "value": x[index],
                    "knn_k": k,
                    "neighbour_count_symmetric": int(weights[index].sum()),
                    "gistar_z": gi[index],
                    "p_hot_one_sided": p_hot[index],
                    "p_cold_one_sided": p_cold[index],
                    "p_two_sided": p_two[index],
                    "q_two_sided_bh_within_outcome": q_two[index],
                    "classification": classification[index],
                }
            )
    global_stats = pd.DataFrame(global_rows)
    global_stats["q_two_sided_bh_across_outcomes"] = bh_adjust(global_stats["permutation_p_two_sided"].to_numpy())
    local_stats = pd.DataFrame(local_rows).merge(
        city_summary[["city", "longitude", "latitude"]], on="city", how="left", validate="many_to_one"
    )
    return global_stats, local_stats, edges, distance, weights


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    definitions, multi_host, high_risk, dual_priority, network_path = load_target_definitions()
    target_ids = multi_host | high_risk
    libraries, rpkm_rows, raw_ids, total_rows, missing_target_ids = load_target_rpkm(target_ids)
    missing_multi_host = multi_host & missing_target_ids
    if missing_multi_host:
        raise AssertionError(f"Multi-host sequences missing from RPKM: {sorted(missing_multi_host)}")
    mapped_high_risk = high_risk & rpkm_rows.keys()
    mapped_dual_priority = dual_priority & rpkm_rows.keys()
    metadata = load_metadata(libraries)
    cities = sorted(metadata["city"].unique())
    coordinates, coordinate_path = load_coordinates(cities)
    city_summary, library_summary, sequence_summary = summarize_cities(
        libraries,
        metadata,
        definitions,
        multi_host,
        mapped_high_risk,
        mapped_dual_priority,
        rpkm_rows,
        raw_ids,
        coordinates,
    )

    global_stats, local_stats, edges, distance, weights = spatial_statistics(
        city_summary, k=4, permutations=99999, seed=123
    )
    sensitivity_rows = []
    for k in [3, 5]:
        sensitivity_global, sensitivity_local, _, _, _ = spatial_statistics(
            city_summary, k=k, permutations=9999, seed=123 + k
        )
        sensitivity_global["analysis"] = f"k={k} sensitivity"
        sensitivity_rows.append(sensitivity_global)
        sensitivity_local = sensitivity_local.copy()
        sensitivity_local["analysis"] = f"k={k} sensitivity"
        sensitivity_rows.append(sensitivity_local)

    city_summary = city_summary.sort_values("surveillance_priority_rank").reset_index(drop=True)
    sequence_summary = sequence_summary.sort_values(
        ["dual_priority", "multi_host_sharing", "calibrated_score_mean", "n_positive_libraries"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)

    missing_table = definitions.loc[definitions["qseqid"].isin(missing_target_ids)].copy()
    missing_table["reason"] = "No one-to-one normalized ID in the final vertebrate-virus RPKM matrix; excluded from city detection metrics"
    missing_table.to_csv(RESULTS / "risk_sequences_not_in_final_rpkm.csv", index=False)

    city_summary.to_csv(RESULTS / "city_surveillance_hotspot_metrics.csv", index=False)
    library_summary.to_csv(RESULTS / "library_priority_virus_metrics.csv", index=False)
    sequence_summary.to_csv(RESULTS / "priority_sequence_distribution.csv", index=False)
    global_stats.to_csv(RESULTS / "global_spatial_autocorrelation.csv", index=False)
    local_stats.to_csv(RESULTS / "local_getis_ord_gistar.csv", index=False)
    edges.to_csv(RESULTS / "spatial_knn_edges.csv", index=False)
    pd.DataFrame(distance, index=cities, columns=cities).to_csv(RESULTS / "city_distance_km.csv")
    pd.DataFrame(weights, index=cities, columns=cities).to_csv(RESULTS / "spatial_knn_weights.csv")

    sensitivity_frames = []
    for frame in sensitivity_rows:
        frame = frame.copy()
        frame["table_type"] = "global" if "moran_i" in frame.columns else "local"
        sensitivity_frames.append(frame)
    pd.concat(sensitivity_frames, ignore_index=True, sort=False).to_csv(
        RESULTS / "knn_sensitivity_statistics.csv", index=False
    )

    top = city_summary.head(8)
    significant_hotspots = local_stats.loc[local_stats["classification"].eq("Hot spot")]
    summary = {
        "n_libraries": len(libraries),
        "n_cities": len(cities),
        "n_risk_scored_sequences": int(definitions["qseqid"].nunique()),
        "n_multi_host_sequences": len(multi_host),
        "n_high_risk_sequences": len(high_risk),
        "n_dual_priority_sequences": len(dual_priority),
        "n_high_risk_sequences_mapped_to_rpkm": len(mapped_high_risk),
        "n_dual_priority_sequences_mapped_to_rpkm": len(mapped_dual_priority),
        "n_risk_sequences_excluded_from_city_metrics": len(missing_target_ids),
        "n_target_sequences_mapped_to_rpkm": len(rpkm_rows),
        "n_significant_hotspot_city_outcome_pairs": int(len(significant_hotspots)),
        "top_surveillance_cities": top[
            [
                "city",
                "surveillance_priority_rank",
                "surveillance_priority_score",
                "multi_host_sequence_richness",
                "multi_host_standardized_richness_33",
                "multi_host_prevalence_percent",
                "multi_host_mean_sequences_per_library",
                "high_risk_sequence_richness",
                "high_risk_standardized_richness_33",
                "high_risk_prevalence_percent",
                "high_risk_mean_sequences_per_library",
                "dual_priority_sequence_richness",
                "dual_priority_standardized_richness_33",
                "dual_priority_prevalence_percent",
                "dual_priority_mean_sequences_per_library",
            ]
        ].to_dict("records"),
        "significant_hotspots": significant_hotspots[
            ["outcome", "city", "value", "gistar_z", "q_two_sided_bh_within_outcome"]
        ].to_dict("records"),
    }
    (RESULTS / "key_results.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
