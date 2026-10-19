#!/usr/bin/env python3
                                                                    

                                                                        
                                                                   
   

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


GROUPS = ("all_viruses", "bacteriophages", "vertebrate_eukaryotic_viruses")


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    input_dir = here / "inputs"
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, default=input_dir / "Supplementary Data 1.xlsx")
    parser.add_argument("--all-matrix", type=Path, default=input_dir / "RPKM_allvirus_matrix_clean.tsv")
    parser.add_argument("--vertebrate-matrix", type=Path, default=input_dir / "RPKM_veter_matrix_clean.tsv")
    parser.add_argument("--environment", type=Path, default=here / "city_environment_indicators.csv")
    parser.add_argument("--output", type=Path, default=here / "results")
    parser.add_argument("--chunk-size", type=int, default=500)
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=123)
    return parser.parse_args()


def clean_city(location: str) -> str:
    return str(location).split(",", 1)[0].strip()


def load_metadata(path: Path) -> pd.DataFrame:
    metadata = pd.read_excel(path, sheet_name=0, header=1)
    required = {"Library", "Location of Sampling (Sampling site, Province)"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing columns: {sorted(missing)}")
    metadata = metadata.copy()
    metadata["library_key"] = metadata["Library"].astype(str).str.strip().str.casefold()
    metadata["city"] = metadata["Location of Sampling (Sampling site, Province)"].map(clean_city)
    if metadata["library_key"].duplicated().any():
        duplicates = metadata.loc[metadata["library_key"].duplicated(False), "Library"].tolist()
        raise ValueError(f"Duplicate library identifiers after case normalization: {duplicates[:10]}")
    if metadata["city"].nunique() != 22:
        raise ValueError(f"Expected 22 cities, found {metadata['city'].nunique()}")
    return metadata


def read_matrix_header(path: Path) -> list[str]:
    return pd.read_csv(path, sep="\t", nrows=0).columns.astype(str).tolist()


def align_samples(metadata: pd.DataFrame, matrix_path: Path) -> tuple[list[str], np.ndarray, list[str], dict[str, int]]:
    header = read_matrix_header(matrix_path)
    if not header or header[0] != "Contig":
        raise ValueError(f"First matrix column must be Contig: {matrix_path}")
    samples = header[1:]
    city_names = sorted(metadata["city"].unique().tolist())
    city_to_idx = {city: i for i, city in enumerate(city_names)}
    metadata_map = dict(zip(metadata["library_key"], metadata["city"]))
    sample_cities = []
    unmatched = []
    for sample in samples:
        city = metadata_map.get(sample.strip().casefold())
        if city is None:
            unmatched.append(sample)
        else:
            sample_cities.append(city_to_idx[city])
    if unmatched:
        raise ValueError(f"Matrix samples absent from metadata: {unmatched[:10]}")
    matrix_keys = {sample.strip().casefold() for sample in samples}
    metadata_only = sorted(set(metadata["library_key"]) - matrix_keys)
    if metadata_only:
        raise ValueError(f"Metadata libraries absent from matrix: {metadata_only[:10]}")
    return samples, np.asarray(sample_cities, dtype=np.int16), city_names, city_to_idx


def expected_detection_lookup(n: int, m: int) -> np.ndarray:
    lookup = np.zeros(n + 1, dtype=float)
    denominator = math.comb(n, m)
    for k in range(1, n + 1):
        if n - k < m:
            lookup[k] = 1.0
        else:
            lookup[k] = 1.0 - (math.comb(n - k, m) / denominator)
    return lookup


def rank_average(values: np.ndarray) -> np.ndarray:
    return pd.Series(np.asarray(values, dtype=float)).rank(method="average").to_numpy()


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x - x.mean()
    y = y - y.mean()
    denom = math.sqrt(float(np.dot(x, x) * np.dot(y, y)))
    return float(np.dot(x, y) / denom) if denom > 0 else float("nan")


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return pearson(rank_average(x), rank_average(y))


def permutation_spearman(x: np.ndarray, y: np.ndarray, permutations: int, rng: np.random.Generator) -> tuple[float, float]:
    rx = rank_average(x)
    ry = rank_average(y)
    observed = pearson(rx, ry)
    exceed = 0
    for _ in range(permutations):
        permuted = pearson(rx, rng.permutation(ry))
        exceed += abs(permuted) >= abs(observed) - 1e-15
    return observed, (exceed + 1.0) / (permutations + 1.0)


def bh_adjust(pvalues: list[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 1.0
    for rank_index in range(len(p) - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p[original_index] * len(p) / rank)
        adjusted[original_index] = running
    return np.clip(adjusted, 0.0, 1.0)


def bray_curtis(profile: np.ndarray) -> np.ndarray:
                                                                              
    n_city = profile.shape[1]
    distance = np.zeros((n_city, n_city), dtype=float)
    for i in range(n_city):
        for j in range(i + 1, n_city):
            numerator = np.abs(profile[:, i] - profile[:, j]).sum()
            denominator = (profile[:, i] + profile[:, j]).sum()
            value = float(numerator / denominator) if denominator > 0 else 0.0
            distance[i, j] = distance[j, i] = value
    return distance


def pcoa(distance: np.ndarray, dimensions: int = 3) -> tuple[np.ndarray, np.ndarray]:
    n = distance.shape[0]
    centering = np.eye(n) - np.ones((n, n)) / n
    gower = -0.5 * centering @ (distance ** 2) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gower)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    positive = eigenvalues > 0
    positive_values = eigenvalues[positive]
    positive_vectors = eigenvectors[:, positive]
    use = min(dimensions, len(positive_values))
    coordinates = positive_vectors[:, :use] * np.sqrt(positive_values[:use])
    explained = positive_values[:use] / positive_values.sum() if positive_values.sum() else np.zeros(use)
    return coordinates, explained


def hat_matrix(design: np.ndarray) -> np.ndarray:
    return design @ np.linalg.pinv(design)


def permanova_continuous(distance: np.ndarray, predictor: np.ndarray, permutations: int, rng: np.random.Generator) -> tuple[float, float, float]:
    n = len(predictor)
    centering = np.eye(n) - np.ones((n, n)) / n
    gower = -0.5 * centering @ (distance ** 2) @ centering
    intercept = np.ones((n, 1))
    h0 = hat_matrix(intercept)

    def statistic(x: np.ndarray) -> tuple[float, float]:
        design = np.column_stack((np.ones(n), x))
        h1 = hat_matrix(design)
        ss_model = float(np.trace((h1 - h0) @ gower))
        ss_residual = float(np.trace((np.eye(n) - h1) @ gower))
        pseudo_f = (ss_model / 1.0) / (ss_residual / (n - 2))
        r2 = ss_model / (ss_model + ss_residual)
        return pseudo_f, r2

    observed, r2 = statistic(np.asarray(predictor, dtype=float))
    exceed = 0
    for _ in range(permutations):
        permuted, _ = statistic(rng.permutation(predictor))
        exceed += permuted >= observed - 1e-15
    pvalue = (exceed + 1.0) / (permutations + 1.0)
    return observed, r2, pvalue


def upper_triangle(matrix: np.ndarray) -> np.ndarray:
    return matrix[np.triu_indices_from(matrix, k=1)]


def euclidean_distance(features: np.ndarray) -> np.ndarray:
    diff = features[:, None, :] - features[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=2))


def mantel_spearman(distance_a: np.ndarray, distance_b: np.ndarray, permutations: int, rng: np.random.Generator) -> tuple[float, float]:
    observed = spearman(upper_triangle(distance_a), upper_triangle(distance_b))
    exceed = 0
    n = distance_a.shape[0]
    for _ in range(permutations):
        order = rng.permutation(n)
        permuted = distance_b[np.ix_(order, order)]
        value = spearman(upper_triangle(distance_a), upper_triangle(permuted))
        exceed += abs(value) >= abs(observed) - 1e-15
    return observed, (exceed + 1.0) / (permutations + 1.0)


def load_environment(path: Path, city_names: list[str]) -> pd.DataFrame:
    environment = pd.read_csv(path)
    required = {
        "city",
        "built_up_km2",
        "population_density_per_km2",
        "green_space_percent",
        "road_density_m_per_km2",
    }
    missing = required.difference(environment.columns)
    if missing:
        raise ValueError(f"Environment table is missing columns: {sorted(missing)}")
    if environment["city"].duplicated().any():
        raise ValueError("Environment table has duplicate city names")
    unknown = sorted(set(environment["city"]) - set(city_names))
    if unknown:
        raise ValueError(f"Environment table contains cities absent from metadata: {unknown}")
    included_city_names = [city for city in city_names if city in set(environment["city"])]
    if len(included_city_names) < 20:
        raise ValueError(f"Too few cities have complete environmental data: {len(included_city_names)}")
    environment = environment.set_index("city").reindex(included_city_names)
    if environment[list(required - {"city"})].isna().any().any():
        bad = environment[environment.isna().any(axis=1)].index.tolist()
        raise ValueError(f"Missing environmental data for cities: {bad}")
    numeric = [
        "built_up_km2",
        "population_density_per_km2",
        "green_space_percent",
        "road_density_m_per_km2",
    ]
    environment[numeric] = environment[numeric].apply(pd.to_numeric, errors="raise")
    if (environment[["built_up_km2", "population_density_per_km2", "road_density_m_per_km2"]] <= 0).any().any():
        raise ValueError("Built-up area, population density, and road density must be positive")
    if ((environment["green_space_percent"] < 0) | (environment["green_space_percent"] > 100)).any():
        raise ValueError("Green-space coverage must be between 0 and 100 percent")
    return environment.reset_index()


def urbanization_framework(environment: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, np.ndarray]:
    transformed = np.column_stack(
        (
            np.log1p(environment["built_up_km2"].to_numpy(float)),
            np.log1p(environment["population_density_per_km2"].to_numpy(float)),
            -environment["green_space_percent"].to_numpy(float),
            np.log1p(environment["road_density_m_per_km2"].to_numpy(float)),
        )
    )
    z = (transformed - transformed.mean(axis=0)) / transformed.std(axis=0, ddof=1)
    urbanization_index = z.mean(axis=1)
    u, singular, vt = np.linalg.svd(z, full_matrices=False)
    scores = u[:, 0] * singular[0]
    loadings = vt[0].copy()
    pressure_mean = z.mean(axis=1)
    if pearson(scores, pressure_mean) < 0:
        scores *= -1
        loadings *= -1
    variance = singular ** 2
    explained = float(variance[0] / variance.sum())
    return urbanization_index, scores, loadings, explained, z


def matrix_pass_one(
    matrix_path: Path,
    euk_ids: set[str],
    samples: list[str],
    sample_city_idx: np.ndarray,
    city_names: list[str],
    chunk_size: int,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, float], dict[str, int], set[str]]:
    n_samples = len(samples)
    n_cities = len(city_names)
    sample_totals = {group: np.zeros(n_samples, dtype=float) for group in GROUPS}
    observed = {group: np.zeros(n_cities, dtype=float) for group in GROUPS}
    rarefied = {group: np.zeros(n_cities, dtype=float) for group in GROUPS}
    sample_counts = np.bincount(sample_city_idx, minlength=n_cities)
    rarefy_n = int(sample_counts.min())
    lookups = [expected_detection_lookup(int(n), rarefy_n) for n in sample_counts]
    row_counts = {group: 0 for group in GROUPS}
    seen_ids: set[str] = set()

    usecols = ["Contig"] + samples
    for chunk in pd.read_csv(matrix_path, sep="\t", usecols=usecols, chunksize=chunk_size, dtype={"Contig": str}):
        ids = chunk["Contig"].astype(str).to_numpy()
        duplicate_ids = seen_ids.intersection(ids)
        if duplicate_ids:
            raise ValueError(f"Duplicate vOTU identifier: {next(iter(duplicate_ids))}")
        seen_ids.update(ids)
        values = chunk[samples].to_numpy(dtype=float, copy=False)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError("Matrix contains missing, infinite, or negative abundance values")
        euk_mask = np.fromiter((identifier in euk_ids for identifier in ids), dtype=bool, count=len(ids))
        group_masks = {
            "all_viruses": np.ones(len(ids), dtype=bool),
            "bacteriophages": ~euk_mask,
            "vertebrate_eukaryotic_viruses": euk_mask,
        }
        for group, mask in group_masks.items():
            group_values = values[mask, :]
            row_counts[group] += int(mask.sum())
            sample_totals[group] += group_values.sum(axis=0)
            present = group_values > 0
            for city_index in range(n_cities):
                columns = sample_city_idx == city_index
                occupancy = present[:, columns].sum(axis=1)
                observed[group][city_index] += np.count_nonzero(occupancy)
                rarefied[group][city_index] += lookups[city_index][occupancy].sum()
    return sample_totals, observed, rarefied, row_counts, seen_ids


def matrix_pass_two(
    matrix_path: Path,
    euk_ids: set[str],
    samples: list[str],
    sample_city_idx: np.ndarray,
    city_names: list[str],
    sample_totals: dict[str, np.ndarray],
    chunk_size: int,
) -> dict[str, np.ndarray]:
    n_cities = len(city_names)
    profile_chunks: dict[str, list[np.ndarray]] = {group: [] for group in GROUPS}
    usecols = ["Contig"] + samples
    for chunk in pd.read_csv(matrix_path, sep="\t", usecols=usecols, chunksize=chunk_size, dtype={"Contig": str}):
        ids = chunk["Contig"].astype(str).to_numpy()
        values = chunk[samples].to_numpy(dtype=float, copy=False)
        euk_mask = np.fromiter((identifier in euk_ids for identifier in ids), dtype=bool, count=len(ids))
        group_masks = {
            "all_viruses": np.ones(len(ids), dtype=bool),
            "bacteriophages": ~euk_mask,
            "vertebrate_eukaryotic_viruses": euk_mask,
        }
        for group, mask in group_masks.items():
            group_values = values[mask, :]
            denominators = sample_totals[group]
            normalized = np.divide(
                group_values,
                denominators[None, :],
                out=np.zeros_like(group_values, dtype=float),
                where=denominators[None, :] > 0,
            )
            city_profile = np.zeros((group_values.shape[0], n_cities), dtype=float)
            for city_index in range(n_cities):
                columns = (sample_city_idx == city_index) & (denominators > 0)
                if not np.any(columns):
                    raise ValueError(f"No non-zero {group} libraries in city index {city_index}")
                city_profile[:, city_index] = normalized[:, columns].mean(axis=1)
            profile_chunks[group].append(city_profile)
    return {group: np.concatenate(chunks, axis=0) for group, chunks in profile_chunks.items()}


def save_distance(path: Path, distance: np.ndarray, city_names: list[str]) -> None:
    pd.DataFrame(distance, index=city_names, columns=city_names).to_csv(path, index_label="city")


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    metadata = load_metadata(args.metadata)
    all_samples, _, all_city_names, _ = align_samples(metadata, args.all_matrix)
    vertebrate_header = read_matrix_header(args.vertebrate_matrix)
    if [x.casefold() for x in vertebrate_header[1:]] != [x.casefold() for x in all_samples]:
        raise ValueError("All-virus and vertebrate-virus matrix sample headers do not align")
    euk_ids = set(pd.read_csv(args.vertebrate_matrix, sep="\t", usecols=["Contig"], dtype=str)["Contig"])

    environment = load_environment(args.environment, all_city_names)
    city_names = environment["city"].tolist()
    city_to_idx = {city: index for index, city in enumerate(city_names)}
    library_to_city = dict(zip(metadata["library_key"], metadata["city"]))
    samples = [sample for sample in all_samples if library_to_city[sample.strip().casefold()] in city_to_idx]
    sample_city_idx = np.asarray(
        [city_to_idx[library_to_city[sample.strip().casefold()]] for sample in samples], dtype=np.int16
    )
    urbanization_index, pc1, loadings, explained, standardized_environment = urbanization_framework(environment)
    environment["urbanization_index"] = urbanization_index
    environment["urbanization_pc1"] = pc1

    sample_totals, observed, rarefied, row_counts, all_ids = matrix_pass_one(
        args.all_matrix,
        euk_ids,
        samples,
        sample_city_idx,
        city_names,
        args.chunk_size,
    )
    missing_euk = euk_ids - all_ids
    if missing_euk:
        raise ValueError(f"Vertebrate matrix has {len(missing_euk)} vOTUs absent from all-virus matrix")
    if row_counts["vertebrate_eukaryotic_viruses"] != len(euk_ids):
        raise ValueError("Vertebrate vOTU count does not match subgroup classification")

    profiles = matrix_pass_two(
        args.all_matrix,
        euk_ids,
        samples,
        sample_city_idx,
        city_names,
        sample_totals,
        args.chunk_size,
    )

    sample_counts = np.bincount(sample_city_idx, minlength=len(city_names))
    city_metrics = environment.copy()
    city_metrics["n_libraries"] = sample_counts
    for group in GROUPS:
        city_metrics[f"n_nonzero_libraries_{group}"] = np.bincount(
            sample_city_idx[sample_totals[group] > 0], minlength=len(city_names)
        )
        city_metrics[f"observed_richness_{group}"] = observed[group]
        city_metrics[f"expected_richness_at_33_{group}"] = rarefied[group]
    city_metrics.to_csv(args.output / "city_level_metrics.csv", index=False, float_format="%.10g")

    loading_names = [
        "log1p_built_up_km2",
        "log1p_population_density_per_km2",
        "negative_green_space_percent",
        "log1p_road_density_m_per_km2",
    ]
    pd.DataFrame({"transformed_indicator": loading_names, "PC1_loading": loadings}).assign(
        PC1_explained_variance=explained
    ).to_csv(args.output / "urbanization_pca_loadings.csv", index=False, float_format="%.10g")

    primary_rows = []
    sensitivity_rows = []
    pcoa_rows = []
    environmental_distance = euclidean_distance(standardized_environment)
    for group in GROUPS:
        richness = rarefied[group]
        rho, richness_p = permutation_spearman(urbanization_index, richness, args.permutations, rng)
        primary_rows.append(
            {
                "viral_group": group,
                "analysis": "urbanization_index_vs_expected_richness_at_33",
                "statistic_name": "Spearman_rho",
                "statistic": rho,
                "effect_size_R2": np.nan,
                "permutation_p": richness_p,
                "permutations": args.permutations,
            }
        )

        distance = bray_curtis(profiles[group])
        save_distance(args.output / f"bray_curtis_{group}.csv", distance, city_names)
        pseudo_f, r2, permanova_p = permanova_continuous(distance, urbanization_index, args.permutations, rng)
        primary_rows.append(
            {
                "viral_group": group,
                "analysis": "urbanization_index_vs_Bray_Curtis_PERMANOVA",
                "statistic_name": "pseudo_F",
                "statistic": pseudo_f,
                "effect_size_R2": r2,
                "permutation_p": permanova_p,
                "permutations": args.permutations,
            }
        )
        mantel_rho, mantel_p = mantel_spearman(distance, environmental_distance, args.permutations, rng)
        primary_rows.append(
            {
                "viral_group": group,
                "analysis": "environmental_distance_vs_Bray_Curtis_Mantel",
                "statistic_name": "Mantel_Spearman_rho",
                "statistic": mantel_rho,
                "effect_size_R2": np.nan,
                "permutation_p": mantel_p,
                "permutations": args.permutations,
            }
        )

        pc1_rho, pc1_richness_p = permutation_spearman(pc1, richness, args.permutations, rng)
        sensitivity_rows.append(
            {
                "viral_group": group,
                "analysis": "PCA_PC1_vs_expected_richness_at_33",
                "statistic_name": "Spearman_rho",
                "statistic": pc1_rho,
                "effect_size_R2": np.nan,
                "permutation_p": pc1_richness_p,
                "permutations": args.permutations,
            }
        )
        pc1_f, pc1_r2, pc1_permanova_p = permanova_continuous(distance, pc1, args.permutations, rng)
        sensitivity_rows.append(
            {
                "viral_group": group,
                "analysis": "PCA_PC1_vs_Bray_Curtis_PERMANOVA",
                "statistic_name": "pseudo_F",
                "statistic": pc1_f,
                "effect_size_R2": pc1_r2,
                "permutation_p": pc1_permanova_p,
                "permutations": args.permutations,
            }
        )

        coordinates, pcoa_explained = pcoa(distance, dimensions=3)
        for city_index, city in enumerate(city_names):
            row = {"city": city, "viral_group": group}
            for axis in range(coordinates.shape[1]):
                row[f"PCoA{axis + 1}"] = coordinates[city_index, axis]
                row[f"PCoA{axis + 1}_explained"] = pcoa_explained[axis]
            pcoa_rows.append(row)

    primary = pd.DataFrame(primary_rows)
    primary.to_csv(args.output / "primary_statistics.csv", index=False, float_format="%.10g")
    pd.DataFrame(sensitivity_rows).to_csv(
        args.output / "pca_sensitivity_statistics.csv", index=False, float_format="%.10g"
    )
    pd.DataFrame(pcoa_rows).to_csv(args.output / "pcoa_coordinates.csv", index=False, float_format="%.10g")

    indicator_columns = [
        "built_up_km2",
        "population_density_per_km2",
        "green_space_percent",
        "road_density_m_per_km2",
    ]
    secondary_rows = []
    for group in GROUPS:
        for indicator in indicator_columns:
            rho, pvalue = permutation_spearman(
                environment[indicator].to_numpy(float), rarefied[group], args.permutations, rng
            )
            secondary_rows.append(
                {
                    "viral_group": group,
                    "indicator": indicator,
                    "Spearman_rho": rho,
                    "permutation_p": pvalue,
                }
            )
    adjusted = bh_adjust([row["permutation_p"] for row in secondary_rows])
    for row, qvalue in zip(secondary_rows, adjusted):
        row["BH_FDR_q"] = qvalue
    pd.DataFrame(secondary_rows).to_csv(
        args.output / "indicator_richness_correlations.csv", index=False, float_format="%.10g"
    )


if __name__ == "__main__":
    main()
