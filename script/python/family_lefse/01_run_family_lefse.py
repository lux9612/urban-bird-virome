                                                                             

                                                                        
                                                                           
                                                                       
                                                                              
                                                                    
                                                                           
   

from __future__ import annotations

import math
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
RESULTS = HERE / "results"
RPKM = INPUTS / "RPKM_veter_matrix_clean.tsv"
TAXONOMY = INPUTS / "vOTU_taxonomy.xlsx"
METADATA = INPUTS / "Supplementary Data 1.xlsx"

KW_ALPHA = 0.05
PAIRWISE_ALPHA = 0.05
DEFAULT_LDA_THRESHOLD = 2.0
HIGH_EFFECT_LDA_THRESHOLD = 4.0
LDA_BOOTSTRAPS = 5000
LDA_FRACTION = 0.67
LEFSE_NORMALIZATION = 1_000_000.0
SEED = 123


def normalize_id(value: object) -> str:
    return str(value).strip().casefold()


def average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + 1 + end) / 2.0
        start = end
    return ranks


def gammaincc(a: float, x: float) -> float:
                                                                             
    if x < 0 or a <= 0:
        raise ValueError("Invalid gamma arguments")
    if x == 0:
        return 1.0
    eps, fpmin, maximum = 3e-14, 1e-300, 10000
    gln = math.lgamma(a)
    if x < a + 1.0:
        ap = a
        term = total = 1.0 / a
        for _ in range(maximum):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * eps:
                break
        lower = total * math.exp(-x + a * math.log(x) - gln)
        return max(0.0, min(1.0, 1.0 - lower))
    b = x + 1.0 - a
    c = 1.0 / fpmin
    d = 1.0 / b
    h = d
    for i in range(1, maximum + 1):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < fpmin:
            d = fpmin
        c = b + an / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return max(0.0, min(1.0, math.exp(-x + a * math.log(x) - gln) * h))


def chi_square_sf(statistic: float, degrees_of_freedom: int) -> float:
    return gammaincc(degrees_of_freedom / 2.0, statistic / 2.0)


def kruskal_wallis(values: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    ranks = average_ranks(values)
    n = len(values)
    groups = np.unique(labels)
    statistic = 12.0 / (n * (n + 1.0)) * sum(
        ranks[labels == group].sum() ** 2 / np.sum(labels == group) for group in groups
    ) - 3.0 * (n + 1.0)
    _, counts = np.unique(values, return_counts=True)
    correction = 1.0 - np.sum(counts**3 - counts) / (n**3 - n)
    if correction <= 0:
        return 0.0, 1.0
    statistic /= correction
    return float(statistic), chi_square_sf(float(statistic), len(groups) - 1)


def mann_whitney_two_sided(first: np.ndarray, second: np.ndarray) -> float:
    values = np.concatenate([first, second])
    ranks = average_ranks(values)
    n1, n2 = len(first), len(second)
    u1 = ranks[:n1].sum() - n1 * (n1 + 1.0) / 2.0
    mean = n1 * n2 / 2.0
    _, counts = np.unique(values, return_counts=True)
    n = n1 + n2
    tie_term = np.sum(counts**3 - counts)
    variance = n1 * n2 / 12.0 * ((n + 1.0) - tie_term / (n * (n - 1.0)))
    if variance <= 0:
        return 1.0
    z = max(0.0, abs(u1 - mean) - 0.5) / math.sqrt(variance)
    return math.erfc(z / math.sqrt(2.0))


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty(len(ranked), dtype=float)
    output[order] = np.clip(adjusted, 0, 1)
    return output


def lda_first_direction(matrix: np.ndarray, labels: np.ndarray) -> np.ndarray:
    classes = np.unique(labels)
    overall = matrix.mean(axis=0)
    within = np.zeros((matrix.shape[1], matrix.shape[1]), dtype=float)
    between = np.zeros_like(within)
    for city in classes:
        block = matrix[labels == city]
        mean = block.mean(axis=0)
        centered = block - mean
        within += centered.T @ centered
        delta = (mean - overall)[:, None]
        between += len(block) * (delta @ delta.T)
    within /= max(len(matrix) - len(classes), 1)
    between /= max(len(matrix), 1)
    eigenvalues, eigenvectors = np.linalg.eigh(within)
    tolerance = max(float(eigenvalues.max()), 1.0) * 1e-12
    keep = eigenvalues > tolerance
    if not np.any(keep):
        return np.ones(matrix.shape[1], dtype=float) / math.sqrt(matrix.shape[1])
    whitening = eigenvectors[:, keep] / np.sqrt(eigenvalues[keep])[None, :]
    reduced = whitening.T @ between @ whitening
    values, vectors = np.linalg.eigh(reduced)
    direction = whitening @ vectors[:, np.argmax(values)]
    norm = np.linalg.norm(direction)
    return direction / norm if norm else np.ones(matrix.shape[1], dtype=float) / math.sqrt(matrix.shape[1])


def lefse_lda_scores(
    matrix: np.ndarray,
    labels: np.ndarray,
    feature_names: list[str],
    bootstraps: int = LDA_BOOTSTRAPS,
    fraction: float = LDA_FRACTION,
    seed: int = SEED,
) -> dict[str, float]:
                                                                              
    rng = np.random.default_rng(seed)
    work = matrix.copy()
    classes = sorted(np.unique(labels).tolist())
    for column in range(work.shape[1]):
        for city in classes:
            mask = labels == city
            block = work[mask, column]
            if len(np.unique(block)) <= max(len(block) * 0.5, 4):
                noise_scale = np.maximum(np.abs(block) * 0.05, 0.01)
                work[mask, column] = np.abs(block + rng.normal(0.0, noise_scale))

    pairs = list(combinations(classes, 2))
    pair_contributions = np.zeros((bootstraps, len(pairs), work.shape[1]), dtype=float)
    sample_count = int(len(work) * fraction)
    completed = 0
    for bootstrap in range(bootstraps):
        for _ in range(1000):
            indices = rng.integers(0, len(work), size=sample_count)
            sampled_labels = labels[indices]
            if all(np.sum(sampled_labels == city) >= 2 for city in classes):
                break
        sampled = work[indices]
        direction = lda_first_direction(sampled, sampled_labels)
        projected = sampled @ direction
        for pair_index, (first, second) in enumerate(pairs):
            first_mask, second_mask = sampled_labels == first, sampled_labels == second
            if not first_mask.any() or not second_mask.any():
                continue
            effect = abs(projected[first_mask].mean() - projected[second_mask].mean())
            coefficients = np.abs(direction * effect)
            mean_difference = np.abs(sampled[first_mask].mean(axis=0) - sampled[second_mask].mean(axis=0))
            pair_contributions[bootstrap, pair_index] = (mean_difference + coefficients) * 0.5
        completed += 1
    if completed != bootstraps:
        raise RuntimeError("Could not complete all LDA bootstrap iterations")
    scores = np.log10(1.0 + pair_contributions.mean(axis=0).max(axis=0))
    return dict(zip(feature_names, scores.astype(float)))


def load_abundance() -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    taxonomy = pd.read_excel(TAXONOMY, sheet_name=0)
    expected_columns = {"seq_name", "taxonomy", "classify"}
    if set(taxonomy.columns) != expected_columns:
        raise ValueError(f"Unexpected taxonomy columns: {taxonomy.columns.tolist()}")
    if taxonomy["seq_name"].duplicated().any():
        raise ValueError("Duplicate vOTU identifiers in taxonomy workbook")
    vertebrate = taxonomy.loc[taxonomy["classify"].eq("veter")].copy()
    if len(vertebrate) != 17_439 or vertebrate["taxonomy"].nunique() != 32:
        raise ValueError("Expected 17,439 vertebrate-associated vOTUs in 32 families")
    mapping = vertebrate.set_index("seq_name")["taxonomy"]

    rpkm = pd.read_csv(RPKM, sep="\t", index_col=0)
    if rpkm.shape != (17_439, 985) or set(rpkm.index) != set(mapping.index):
        raise ValueError("Unexpected RPKM shape or incomplete taxonomy match")
    family_rpkm = rpkm.groupby(mapping.loc[rpkm.index].to_numpy()).sum().sort_index()
    totals = family_rpkm.sum(axis=0)
    relative = family_rpkm.div(totals.replace(0, np.nan), axis=1).fillna(0.0) * LEFSE_NORMALIZATION

    metadata = pd.read_excel(METADATA, sheet_name=0, header=1)
    if metadata["Library"].duplicated().any():
        raise ValueError("Duplicate Library identifiers in Supplementary Data 1")
    metadata["library_key"] = metadata["Library"].map(normalize_id)
    metadata["city"] = metadata["Location of Sampling (Sampling site, Province)"].astype("string").str.strip()
    lookup = {normalize_id(column): column for column in relative.columns}
    if set(metadata["library_key"]) != set(lookup):
        missing_meta = sorted(set(lookup) - set(metadata["library_key"]))
        missing_matrix = sorted(set(metadata["library_key"]) - set(lookup))
        raise ValueError(f"Library mismatch: matrix-only={missing_meta[:5]}, metadata-only={missing_matrix[:5]}")
    metadata["matrix_library"] = metadata["library_key"].map(lookup)
    if metadata["city"].isna().any() or metadata["city"].nunique() != 22:
        raise ValueError("Expected complete assignments to 22 cities")
    ordered_libraries = metadata["matrix_library"].tolist()
    relative = relative[ordered_libraries]
    city = pd.Series(metadata["city"].astype(str).to_numpy(), index=ordered_libraries, name="city")
    return relative, metadata, city


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    relative, metadata, city = load_abundance()
    family_names = relative.index.tolist()
    labels = city.loc[relative.columns].to_numpy(str)
    matrix = relative.T.to_numpy(float)

    kw_rows = []
    for column, family in enumerate(family_names):
        statistic, p_value = kruskal_wallis(matrix[:, column], labels)
        means = {name: float(matrix[labels == name, column].mean()) for name in sorted(np.unique(labels))}
        medians = {name: float(np.median(matrix[labels == name, column])) for name in sorted(np.unique(labels))}
        enriched_city = max(means, key=means.get)
        enriched_mask = labels == enriched_city
        other_mask = ~enriched_mask
        pairwise_p = []
        pairwise_direction = []
        for other_city in sorted(set(means) - {enriched_city}):
            first = matrix[labels == enriched_city, column]
            second = matrix[labels == other_city, column]
            pairwise_p.append(mann_whitney_two_sided(first, second))
            pairwise_direction.append(float(np.median(first)) > float(np.median(second)))
        kw_rows.append(
            {
                "family": family,
                "kw_statistic": statistic,
                "kw_p": p_value,
                "enriched_city": enriched_city,
                "enriched_city_mean_per_million": means[enriched_city],
                "enriched_city_median_per_million": medians[enriched_city],
                "enriched_city_mean_percent": means[enriched_city] / 10_000.0,
                "enriched_city_median_percent": medians[enriched_city] / 10_000.0,
                "enriched_city_prevalence_rpkm_gt_0": float(np.mean(matrix[enriched_mask, column] > 0)),
                "other_cities_mean_per_million": float(matrix[other_mask, column].mean()),
                "mean_fold_enrichment_vs_other_cities": (
                    float(matrix[enriched_mask, column].mean() / matrix[other_mask, column].mean())
                    if matrix[other_mask, column].mean() > 0
                    else float("inf")
                ),
                "overall_mean_per_million": float(matrix[:, column].mean()),
                "overall_prevalence_rpkm_gt_0": float(np.mean(matrix[:, column] > 0)),
                "minimum_pairwise_p_vs_other_cities": float(min(pairwise_p)),
                "maximum_pairwise_p_vs_other_cities": float(max(pairwise_p)),
                "pairwise_all_others_p_lt_0_05_and_median_higher": bool(
                    all(p < PAIRWISE_ALPHA for p in pairwise_p) and all(pairwise_direction)
                ),
            }
        )
    results = pd.DataFrame(kw_rows)
    results["kw_q_bh"] = benjamini_hochberg(results["kw_p"].to_numpy(float))

    selected_names = results.loc[results["kw_p"] < KW_ALPHA, "family"].tolist()
    if selected_names:
        selected_indices = [family_names.index(name) for name in selected_names]
        scores = lefse_lda_scores(matrix[:, selected_indices], labels, selected_names)
    else:
        scores = {}
    results["lda_log10_score"] = results["family"].map(scores)
    results["lefse_default_signature"] = (
        (results["kw_p"] < KW_ALPHA) & (results["lda_log10_score"] > DEFAULT_LDA_THRESHOLD)
    )
    results["lefse_high_effect_signature"] = (
        (results["kw_q_bh"] < 0.05) & (results["lda_log10_score"] > HIGH_EFFECT_LDA_THRESHOLD)
    )
    results["strict_pairwise_consistent_signature"] = (
        results["lefse_default_signature"] & results["pairwise_all_others_p_lt_0_05_and_median_higher"]
    )
    results = results.sort_values(
        ["lefse_high_effect_signature", "lda_log10_score", "kw_p"],
        ascending=[False, False, True],
        na_position="last",
    )
    results.to_csv(RESULTS / "family_lefse_all_features.csv", index=False, encoding="utf-8-sig")
    default_signatures = results.loc[results["lefse_default_signature"]].copy()
    high_effect_signatures = results.loc[results["lefse_high_effect_signature"]].copy()
    default_signatures.to_csv(
        RESULTS / "family_lefse_default_signatures.csv", index=False, encoding="utf-8-sig"
    )
    high_effect_signatures.to_csv(
        RESULTS / "family_lefse_signatures.csv", index=False, encoding="utf-8-sig"
    )

    city_means = pd.DataFrame(
        {name: relative.loc[:, city.index[city.eq(name)]].mean(axis=1) for name in sorted(city.unique())}
    )
    city_means.rename_axis("family").reset_index().to_csv(
        RESULTS / "family_city_mean_relative_abundance_per_million.csv", index=False, encoding="utf-8-sig"
    )
    city_counts = city.value_counts().sort_index().rename("library_count").rename_axis("city").reset_index()
    city_counts.to_csv(RESULTS / "city_library_counts.csv", index=False, encoding="utf-8-sig")

    lefse_input = relative.copy()
    with (RESULTS / "family_lefse_input.tsv").open("w", encoding="utf-8", newline="") as handle:
        handle.write("class\t" + "\t".join(labels) + "\n")
        handle.write("sample\t" + "\t".join(relative.columns.astype(str)) + "\n")
        for family, row in lefse_input.iterrows():
            handle.write(family + "\t" + "\t".join(f"{value:.12g}" for value in row.to_numpy(float)) + "\n")

    genus_status = pd.DataFrame(
        [
            {
                "requested_level": "family",
                "mapping_available": True,
                "annotated_votu_count": 17_439,
                "feature_count": 32,
                "analysis_status": "completed",
            },
            {
                "requested_level": "genus",
                "mapping_available": False,
                "annotated_votu_count": 0,
                "feature_count": 0,
                "analysis_status": "mapping unavailable",
            },
        ]
    )
    genus_status.to_csv(RESULTS / "taxonomic_level_availability.csv", index=False, encoding="utf-8-sig")


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
