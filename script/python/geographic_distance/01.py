#!/usr/bin/env python3
                                                                      

                                                                              
                                                                           
                                                                        
                                                            
   

from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
SOURCE_DIR = HERE / "inputs"
METADATA = SOURCE_DIR / "Supplementary Data 1.xlsx"
COORDINATES = SOURCE_DIR / "city_coordinates.csv"
ALL_MATRIX = SOURCE_DIR / "RPKM_allvirus_matrix_clean.tsv"

INPUT_DIR = HERE / "derived_inputs"
DATA_DIR = HERE / "results"

GROUPS = ("all_viruses",)
GROUP_LABELS = {
    "all_viruses": "All viruses",
}
PERMUTATIONS = 5000
SEED = 123


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def expected_detection_lookup(n: int, m: int) -> np.ndarray:
                                                                                   
    lookup = np.zeros(n + 1, dtype=float)
    denominator = math.comb(n, m)
    for k in range(1, n + 1):
        lookup[k] = 1.0 if n - k < m else 1.0 - math.comb(n - k, m) / denominator
    return lookup


def average_ranks(values: np.ndarray) -> np.ndarray:
    return pd.Series(np.asarray(values, dtype=float)).rank(method="average").to_numpy()


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x0 = np.asarray(x, dtype=float) - np.mean(x)
    y0 = np.asarray(y, dtype=float) - np.mean(y)
    denominator = math.sqrt(float(np.dot(x0, x0) * np.dot(y0, y0)))
    return float(np.dot(x0, y0) / denominator)


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return pearson(average_ranks(x), average_ranks(y))


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 1.0
    for reverse_rank in range(len(p) - 1, -1, -1):
        original = order[reverse_rank]
        rank = reverse_rank + 1
        running = min(running, p[original] * len(p) / rank)
        adjusted[original] = running
    return np.clip(adjusted, 0.0, 1.0)


def haversine_matrix(latitude: np.ndarray, longitude: np.ndarray) -> np.ndarray:
    radius_km = 6371.0088
    lat = np.radians(np.asarray(latitude, dtype=float))
    lon = np.radians(np.asarray(longitude, dtype=float))
    dlat = lat[:, None] - lat[None, :]
    dlon = lon[:, None] - lon[None, :]
    a = (
        np.sin(dlat / 2.0) ** 2
        + np.cos(lat[:, None]) * np.cos(lat[None, :]) * np.sin(dlon / 2.0) ** 2
    )
    return 2.0 * radius_km * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def bray_curtis(profile: np.ndarray) -> np.ndarray:
    n_city = profile.shape[1]
    result = np.zeros((n_city, n_city), dtype=float)
    for i in range(n_city):
        for j in range(i + 1, n_city):
            denominator = np.sum(profile[:, i] + profile[:, j])
            value = np.sum(np.abs(profile[:, i] - profile[:, j])) / denominator
            result[i, j] = result[j, i] = float(value)
    return result


def mantel_spearman(
    first: np.ndarray,
    second: np.ndarray,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
                                                                            
    upper = np.triu_indices_from(first, k=1)
    first_rank = average_ranks(first[upper])
    second_rank_values = average_ranks(second[upper])
    second_rank_matrix = np.zeros_like(second, dtype=float)
    second_rank_matrix[upper] = second_rank_values
    second_rank_matrix[(upper[1], upper[0])] = second_rank_values
    observed = pearson(first_rank, second_rank_values)
    exceed = 0
    for _ in range(permutations):
        order = rng.permutation(first.shape[0])
        permuted = second_rank_matrix[np.ix_(order, order)][upper]
        exceed += abs(pearson(first_rank, permuted)) >= abs(observed) - 1e-15
    return observed, (exceed + 1.0) / (permutations + 1.0)


def permutation_spearman(
    first: np.ndarray,
    second: np.ndarray,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    first_rank = average_ranks(first)
    second_rank = average_ranks(second)
    observed = pearson(first_rank, second_rank)
    exceed = 0
    for _ in range(permutations):
        exceed += abs(pearson(first_rank, rng.permutation(second_rank))) >= abs(observed) - 1e-15
    return observed, (exceed + 1.0) / (permutations + 1.0)


def load_and_align_metadata() -> tuple[pd.DataFrame, list[str], np.ndarray, list[str]]:
    metadata = pd.read_excel(METADATA, sheet_name=0, header=1)
    required = {"Library", "Location of Sampling (Sampling site, Province)"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Missing metadata columns: {sorted(missing)}")
    metadata = metadata.copy()
    metadata["library_key"] = metadata["Library"].astype(str).str.strip().str.casefold()
    metadata["city"] = (
        metadata["Location of Sampling (Sampling site, Province)"]
        .astype(str)
        .str.split(",", n=1)
        .str[0]
        .str.strip()
    )
    if metadata["library_key"].duplicated().any():
        raise ValueError("Duplicate library identifiers in Supplementary Data 1")
    header = pd.read_csv(ALL_MATRIX, sep="\t", nrows=0).columns.astype(str).tolist()
    if header[0] != "Contig":
        raise ValueError("All-virus matrix first column is not Contig")
    samples = header[1:]
    city_names = sorted(metadata["city"].unique().tolist())
    if len(metadata) != 985 or len(samples) != 985 or len(city_names) != 22:
        raise ValueError(
            f"Expected 985 metadata rows, 985 matrix libraries, and 22 cities; "
            f"found {len(metadata)}, {len(samples)}, and {len(city_names)}"
        )
    library_city = dict(zip(metadata["library_key"], metadata["city"]))
    unmatched = [sample for sample in samples if sample.strip().casefold() not in library_city]
    if unmatched:
        raise ValueError(f"Matrix libraries absent from metadata: {unmatched[:5]}")
    city_index = {city: index for index, city in enumerate(city_names)}
    sample_city = np.asarray(
        [city_index[library_city[sample.strip().casefold()]] for sample in samples],
        dtype=np.int16,
    )
    return metadata, samples, sample_city, city_names


def load_city_coordinates(metadata: pd.DataFrame, city_names: list[str]) -> pd.DataFrame:
    coordinates = pd.read_csv(COORDINATES)
    required = {"city", "longitude", "latitude"}
    missing = required.difference(coordinates.columns)
    if missing:
        raise ValueError(f"Missing coordinate columns: {sorted(missing)}")
    coordinates = coordinates.copy()
    coordinates["city"] = coordinates["city"].astype(str).str.strip()
    coordinates[["longitude", "latitude"]] = coordinates[
        ["longitude", "latitude"]
    ].apply(pd.to_numeric, errors="raise")
    if coordinates["city"].duplicated().any():
        raise ValueError("Duplicate city identifiers in coordinate data")
    city = coordinates.set_index("city").reindex(city_names).reset_index()
    if city[["longitude", "latitude"]].isna().any().any():
        raise ValueError("Coordinates are not available for every city")
    return city


def matrix_passes(
    samples: list[str],
    sample_city: np.ndarray,
    city_names: list[str],
    chunk_size: int = 1000,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, np.ndarray], dict[str, int]]:
    n_samples = len(samples)
    n_city = len(city_names)
    sample_counts = np.bincount(sample_city, minlength=n_city)
    standard_n = int(sample_counts.min())
    lookups = [expected_detection_lookup(int(n), standard_n) for n in sample_counts]

    observed = {group: np.zeros(n_city, dtype=float) for group in GROUPS}
    expected = {group: np.zeros(n_city, dtype=float) for group in GROUPS}
    sample_totals = {group: np.zeros(n_samples, dtype=float) for group in GROUPS}
    row_counts = {group: 0 for group in GROUPS}
    seen_ids: set[str] = set()
    columns = ["Contig", *samples]

    for chunk in pd.read_csv(
        ALL_MATRIX,
        sep="\t",
        usecols=columns,
        chunksize=chunk_size,
        dtype={"Contig": str},
    ):
        ids = chunk["Contig"].astype(str).to_numpy()
        duplicate = seen_ids.intersection(ids)
        if duplicate:
            raise ValueError(f"Duplicate all-virus vOTU: {next(iter(duplicate))}")
        seen_ids.update(ids)
        values = chunk[samples].to_numpy(dtype=float, copy=False)
        if not np.isfinite(values).all() or np.any(values < 0):
            raise ValueError("All-virus matrix contains invalid abundance values")
        for group in GROUPS:
            group_values = values
            row_counts[group] += len(ids)
            sample_totals[group] += group_values.sum(axis=0)
            present = group_values > 0
            for city_idx in range(n_city):
                occupancy = present[:, sample_city == city_idx].sum(axis=1)
                observed[group][city_idx] += np.count_nonzero(occupancy)
                expected[group][city_idx] += lookups[city_idx][occupancy].sum()

    profile_chunks: dict[str, list[np.ndarray]] = {group: [] for group in GROUPS}
    for chunk in pd.read_csv(
        ALL_MATRIX,
        sep="\t",
        usecols=columns,
        chunksize=chunk_size,
        dtype={"Contig": str},
    ):
        values = chunk[samples].to_numpy(dtype=float, copy=False)
        for group in GROUPS:
            group_values = values
            totals = sample_totals[group]
            normalized = np.divide(
                group_values,
                totals[None, :],
                out=np.zeros_like(group_values, dtype=float),
                where=totals[None, :] > 0,
            )
            city_profile = np.zeros((group_values.shape[0], n_city), dtype=float)
            for city_idx in range(n_city):
                selected = (sample_city == city_idx) & (totals > 0)
                if not np.any(selected):
                    raise ValueError(f"No non-zero {group} libraries in {city_names[city_idx]}")
                city_profile[:, city_idx] = normalized[:, selected].mean(axis=1)
            profile_chunks[group].append(city_profile)

    profiles = {
        group: np.concatenate(group_chunks, axis=0)
        for group, group_chunks in profile_chunks.items()
    }
    return observed, expected, profiles, row_counts


def write_distance(path: Path, values: np.ndarray, city_names: list[str]) -> None:
    pd.DataFrame(values, index=city_names, columns=city_names).to_csv(
        path, index_label="city", float_format="%.10g"
    )


def main() -> None:
    for directory in (INPUT_DIR, DATA_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    metadata, samples, sample_city, city_names = load_and_align_metadata()

    city_coordinates = load_city_coordinates(metadata, city_names)
    city_coordinates["n_libraries"] = np.bincount(sample_city, minlength=len(city_names))
    city_coordinates["n_libraries_with_city_coordinate"] = city_coordinates["n_libraries"]
    city_coordinates["coordinate_source"] = (
        "Sampling-location information; one coordinate pair per city"
    )
    city_coordinates.to_csv(INPUT_DIR / "city_coordinates.csv", index=False, float_format="%.10g")

    geographic = haversine_matrix(
        city_coordinates["latitude"].to_numpy(),
        city_coordinates["longitude"].to_numpy(),
    )
    write_distance(INPUT_DIR / "geographic_distance_km.csv", geographic, city_names)

    observed, expected, profiles, row_counts = matrix_passes(samples, sample_city, city_names)
    sample_counts = np.bincount(sample_city, minlength=len(city_names))
    standard_n = int(sample_counts.min())

    metrics = city_coordinates[
        ["city", "longitude", "latitude", "n_libraries"]
    ].copy()
    for group in GROUPS:
        metrics[f"expected_richness_at_{standard_n}_{group}"] = expected[group]
    metrics.to_csv(INPUT_DIR / "city_standardized_richness.csv", index=False, float_format="%.10g")

    rng = np.random.default_rng(SEED)
    statistical_rows: list[dict[str, object]] = []
    pairwise_rows: list[dict[str, object]] = []
    bray_by_group: dict[str, np.ndarray] = {}
    upper = np.triu_indices(len(city_names), k=1)

    for group in GROUPS:
        richness = expected[group]
        richness_difference = np.abs(richness[:, None] - richness[None, :])
        log_difference = np.abs(np.log1p(richness)[:, None] - np.log1p(richness)[None, :])
        bray = bray_curtis(profiles[group])
        bray_by_group[group] = bray
        write_distance(INPUT_DIR / f"bray_curtis_{group}.csv", bray, city_names)

                                                                            
                                                                             
        for analysis, outcome, family in (
            (
                "Geographic distance vs absolute standardized-richness difference",
                richness_difference,
                "Primary",
            ),
            (
                "Geographic distance vs absolute log1p standardized-richness difference",
                log_difference,
                "Sensitivity",
            ),
            (
                "Geographic distance vs Bray-Curtis dissimilarity",
                bray,
                "Primary",
            ),
        ):
            rho, p_value = mantel_spearman(geographic, outcome, PERMUTATIONS, rng)
            statistical_rows.append(
                {
                    "test_family": family,
                    "viral_group": GROUP_LABELS[group],
                    "analysis": analysis,
                    "n_cities": len(city_names),
                    "n_city_pairs": len(upper[0]),
                    "statistic": "Mantel Spearman rho",
                    "rho": rho,
                    "permutation_p": p_value,
                    "permutations": PERMUTATIONS,
                    "BH_q_across_primary_tests": np.nan,
                }
            )

        for axis in ("latitude", "longitude"):
            rho, p_value = permutation_spearman(
                city_coordinates[axis].to_numpy(), richness, PERMUTATIONS, rng
            )
            statistical_rows.append(
                {
                    "test_family": "Exploratory",
                    "viral_group": GROUP_LABELS[group],
                    "analysis": f"{axis.capitalize()} vs standardized richness",
                    "n_cities": len(city_names),
                    "n_city_pairs": np.nan,
                    "statistic": "Spearman rho",
                    "rho": rho,
                    "permutation_p": p_value,
                    "permutations": PERMUTATIONS,
                    "BH_q_across_primary_tests": np.nan,
                }
            )

        for i, j in zip(*upper):
            pairwise_rows.append(
                {
                    "viral_group": GROUP_LABELS[group],
                    "city_1": city_names[i],
                    "city_2": city_names[j],
                    "geographic_distance_km": geographic[i, j],
                    f"absolute_expected_richness_difference_at_{standard_n}": richness_difference[i, j],
                    f"absolute_log1p_expected_richness_difference_at_{standard_n}": log_difference[i, j],
                    "Bray_Curtis_dissimilarity": bray[i, j],
                }
            )

    statistics = pd.DataFrame(statistical_rows)
    primary = statistics["test_family"].eq("Primary")
    statistics.loc[primary, "BH_q_across_primary_tests"] = bh_adjust(
        statistics.loc[primary, "permutation_p"].to_numpy()
    )
    statistics.to_csv(
        DATA_DIR / "geographic_association_statistics.csv",
        index=False,
        float_format="%.10g",
    )
    pd.DataFrame(pairwise_rows).to_csv(
        DATA_DIR / "pairwise_geographic_virome_data.csv",
        index=False,
        float_format="%.10g",
    )


if __name__ == "__main__":
    main()
