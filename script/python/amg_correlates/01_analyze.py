from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
RESULTS = HERE / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

SEED = 123
PERMUTATIONS = 100_000
BOOTSTRAPS = 20_000


def rankdata(values: np.ndarray) -> np.ndarray:
    return pd.Series(np.asarray(values, dtype=float)).rank(method="average").to_numpy(dtype=float)


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x - x.mean()
    y = y - y.mean()
    denominator = float(np.sqrt(np.sum(x * x) * np.sum(y * y)))
    return float(np.sum(x * y) / denominator) if denominator else float("nan")


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return pearson(rankdata(x), rankdata(y))


def spearman_permutation(
    x: np.ndarray,
    y: np.ndarray,
    rng: np.random.Generator,
    n_perm: int = PERMUTATIONS,
    strata: np.ndarray | None = None,
) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    if strata is not None:
        strata = np.asarray(strata, dtype=object)[valid]

    rx = rankdata(x).copy()
    ry = rankdata(y).copy()
    if strata is None:
        rx -= rx.mean()
        ry -= ry.mean()
        denominator = float(np.sqrt(np.sum(rx * rx) * np.sum(ry * ry)))
        observed = float(np.sum(rx * ry) / denominator) if denominator else float("nan")
        extreme = 0
        remaining = n_perm
        while remaining:
            batch = min(5_000, remaining)
            order = np.argsort(rng.random((batch, len(ry))), axis=1)
            correlations = (ry[order] @ rx) / denominator
            extreme += int(np.sum(np.abs(correlations) >= abs(observed) - 1e-15))
            remaining -= batch
    else:
        groups = [np.flatnonzero(strata == group) for group in pd.unique(strata)]
        rx_centered = rx.copy()
        ry_centered = ry.copy()
        for idx in groups:
            rx_centered[idx] -= rx_centered[idx].mean()
            ry_centered[idx] -= ry_centered[idx].mean()
        denominator = float(
            np.sqrt(np.sum(rx_centered * rx_centered) * np.sum(ry_centered * ry_centered))
        )
        observed = (
            float(np.sum(rx_centered * ry_centered) / denominator)
            if denominator
            else float("nan")
        )
        extreme = 0
        remaining = n_perm
        while remaining:
            batch = min(250, remaining)
            permuted = np.tile(ry_centered, (batch, 1))
            for idx in groups:
                order = np.argsort(rng.random((batch, len(idx))), axis=1)
                permuted[:, idx] = ry_centered[idx][order]
            correlations = (permuted @ rx_centered) / denominator
            extreme += int(np.sum(np.abs(correlations) >= abs(observed) - 1e-15))
            remaining -= batch
    return observed, float((extreme + 1) / (n_perm + 1))


def bootstrap_spearman_ci(
    x: np.ndarray,
    y: np.ndarray,
    rng: np.random.Generator,
    n_boot: int = BOOTSTRAPS,
) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    estimates: list[float] = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(x), len(x))
        value = spearman(x[idx], y[idx])
        if math.isfinite(value):
            estimates.append(value)
    return tuple(np.quantile(estimates, [0.025, 0.975]).astype(float))


def bh_adjust(values: list[float]) -> list[float]:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result.tolist()


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    n1 = a + b
    successes = a + c
    total = a + b + c + d
    denominator = math.comb(total, n1)

    def probability(x: int) -> float:
        return math.comb(successes, x) * math.comb(total - successes, n1 - x) / denominator

    low = max(0, n1 - (total - successes))
    high = min(n1, successes)
    observed_probability = probability(a)
    p_value = sum(
        probability(x)
        for x in range(low, high + 1)
        if probability(x) <= observed_probability + 1e-15
    )
    odds_ratio = (a * d) / (b * c) if b * c else float("inf")
    return float(odds_ratio), float(min(1.0, p_value))


def rank_sum_permutation(
    x: np.ndarray,
    y: np.ndarray,
    rng: np.random.Generator,
    n_perm: int = PERMUTATIONS,
) -> tuple[float, float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    pooled = np.concatenate([x, y])
    ranks = rankdata(pooled)
    n_x = len(x)
    rank_sum_x = float(ranks[:n_x].sum())
    u_x = rank_sum_x - n_x * (n_x + 1) / 2
    rank_biserial = 2 * u_x / (len(x) * len(y)) - 1
    observed = abs(ranks[:n_x].mean() - ranks[n_x:].mean())
    total_rank_sum = ranks.sum()
    extreme = 0
    remaining = n_perm
    while remaining:
        batch = min(2_000, remaining)
        chosen = np.argsort(rng.random((batch, len(pooled))), axis=1)[:, :n_x]
        chosen_sum = ranks[chosen].sum(axis=1)
        other_sum = total_rank_sum - chosen_sum
        differences = np.abs(chosen_sum / n_x - other_sum / len(y))
        extreme += int(np.sum(differences >= observed - 1e-15))
        remaining -= batch
    return float(u_x), float(rank_biserial), float((extreme + 1) / (n_perm + 1))


def normalize_city(value: object) -> str:
    text = str(value).strip()
    aliases = {"Urumqi": "Ürümqi", "Urumchi": "Ürümqi"}
    return aliases.get(text, text)


def load_sample_amg() -> pd.DataFrame:
    df = pd.read_excel(INPUTS / "Supplementary Data 14.xlsx", sheet_name="Sample_AMG_abundance")
    df.columns = [str(x).strip() for x in df.columns]
    df = df.rename(columns={"Total_AMG": "total_amg", "Sample": "sample", "Region": "city"})
    df["sample"] = df["sample"].astype(str).str.strip()
    df["city"] = df["city"].map(normalize_city)
    df["total_amg"] = pd.to_numeric(df["total_amg"], errors="raise")
    if len(df) != 985 or df["sample"].nunique() != 985:
        raise ValueError("Expected 985 unique libraries in Supplementary Data 14.")
    return df[["sample", "city", "total_amg"]]


def derive_non_amg_phage_metrics(samples: list[str]) -> tuple[pd.DataFrame, dict[str, int]]:
    matrix_path = INPUTS / "RPKM_allvirus_matrix_clean.tsv"
    vertebrate_path = INPUTS / "RPKM_veter_matrix_clean.tsv"
    lifestyle_path = INPUTS / "phage_type.xlsx"

    euk_ids = set(pd.read_csv(vertebrate_path, sep="\t", usecols=["Contig"])["Contig"].astype(str))
    lifestyle = pd.read_excel(lifestyle_path)
    lifestyle.columns = [str(x).strip() for x in lifestyle.columns]
    amg_ids = set(lifestyle["Contig"].astype(str).str.strip())

    matrix_samples = pd.read_csv(matrix_path, sep="\t", nrows=0).columns[1:].tolist()
    matrix_by_normalized = {name.casefold(): name for name in matrix_samples}
    requested = [name.casefold() for name in samples]
    if set(matrix_by_normalized) != set(requested):
        raise ValueError("RPKM matrix and Supplementary Data 14 sample IDs do not match.")
    ordered = [matrix_by_normalized[name] for name in requested]

    total_phage_abundance = np.zeros(len(samples), dtype=float)
    total_phage_richness = np.zeros(len(samples), dtype=np.int64)
    non_amg_abundance = np.zeros(len(samples), dtype=float)
    non_amg_richness = np.zeros(len(samples), dtype=np.int64)
    amg_carrier_abundance = np.zeros(len(samples), dtype=float)
    found_amg: set[str] = set()
    found_phage: set[str] = set()

    for chunk in pd.read_csv(matrix_path, sep="\t", chunksize=1_000):
        chunk = chunk[[chunk.columns[0], *ordered]]
        contigs = chunk.iloc[:, 0].astype(str)
        numeric = chunk.iloc[:, 1:].to_numpy(dtype=float, copy=False)
        phage_mask = ~contigs.isin(euk_ids).to_numpy()
        amg_mask = contigs.isin(amg_ids).to_numpy()
        non_amg_mask = phage_mask & ~amg_mask

        if phage_mask.any():
            values = numeric[phage_mask]
            total_phage_abundance += values.sum(axis=0)
            total_phage_richness += (values > 0).sum(axis=0)
            found_phage.update(contigs[phage_mask].tolist())
        if non_amg_mask.any():
            values = numeric[non_amg_mask]
            non_amg_abundance += values.sum(axis=0)
            non_amg_richness += (values > 0).sum(axis=0)
        if amg_mask.any():
            amg_carrier_abundance += numeric[amg_mask].sum(axis=0)
            found_amg.update(contigs[amg_mask].tolist())

    if len(found_phage) != 34_875:
        raise ValueError(f"Expected 34,875 phage vOTUs; found {len(found_phage)}.")
    if found_amg != amg_ids:
        raise ValueError(f"AMG carrier coverage mismatch: {len(found_amg)}/{len(amg_ids)}.")

    output = pd.DataFrame(
        {
            "sample": samples,
            "total_phage_abundance": total_phage_abundance,
            "total_phage_richness": total_phage_richness,
            "non_amg_phage_abundance": non_amg_abundance,
            "non_amg_phage_richness": non_amg_richness,
            "amg_carrier_abundance": amg_carrier_abundance,
        }
    )
    qc = {
        "libraries": len(samples),
        "phage_votus": len(found_phage),
        "amg_carrier_votus": len(found_amg),
        "non_amg_phage_votus": len(found_phage - amg_ids),
    }
    return output, qc


def city_test_rows(
    city: pd.DataFrame,
    outcome: str,
    predictors: list[str],
    role: str,
    family: str,
    rng: np.random.Generator,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for predictor in predictors:
        subset = city[[outcome, predictor]].dropna()
        rho, p_value = spearman_permutation(
            subset[outcome].to_numpy(), subset[predictor].to_numpy(), rng
        )
        ci_low, ci_high = bootstrap_spearman_ci(
            subset[outcome].to_numpy(), subset[predictor].to_numpy(), rng
        )
        rows.append(
            {
                "analysis_level": "City",
                "analysis_role": role,
                "outcome": outcome,
                "predictor": predictor,
                "n": len(subset),
                "method": f"Two-sided Spearman; {PERMUTATIONS:,} city-label permutations",
                "effect": rho,
                "effect_name": "Spearman rho",
                "ci95_low": ci_low,
                "ci95_high": ci_high,
                "permutation_p": p_value,
                "correction_family": family,
            }
        )
    return rows


def main() -> None:
    rng = np.random.default_rng(SEED)
    sample = load_sample_amg()
    metrics, qc = derive_non_amg_phage_metrics(sample["sample"].tolist())
    sample = sample.merge(metrics, on="sample", how="left", validate="one_to_one")
    sample["amg_detected"] = (sample["total_amg"] > 0).astype(int)
    sample["log1p_total_amg"] = np.log1p(sample["total_amg"])

    city = (
        sample.groupby("city", sort=True)
        .agg(
            libraries=("sample", "size"),
            amg_positive_libraries=("amg_detected", "sum"),
            amg_prevalence=("amg_detected", "mean"),
            mean_total_amg=("total_amg", "mean"),
            median_total_amg=("total_amg", "median"),
            mean_non_amg_phage_abundance=("non_amg_phage_abundance", "mean"),
            median_non_amg_phage_abundance=("non_amg_phage_abundance", "median"),
            mean_non_amg_phage_richness=("non_amg_phage_richness", "mean"),
            median_non_amg_phage_richness=("non_amg_phage_richness", "median"),
        )
        .reset_index()
    )

    positive = sample.loc[sample["amg_detected"] == 1].copy()
    positive_city = (
        positive.groupby("city", sort=True)
        .agg(
            positive_libraries=("sample", "size"),
            positive_mean_log1p_amg=("log1p_total_amg", "mean"),
            positive_median_total_amg=("total_amg", "median"),
            positive_mean_non_amg_phage_abundance=("non_amg_phage_abundance", "mean"),
            positive_mean_non_amg_phage_richness=("non_amg_phage_richness", "mean"),
        )
        .reset_index()
    )
    city = city.merge(positive_city, on="city", how="left")

    env = pd.read_csv(INPUTS / "city_environment_indicators.csv")
    env["city"] = env["city"].map(normalize_city)
    urban = pd.read_csv(INPUTS / "city_level_metrics.csv")
    urban["city"] = urban["city"].map(normalize_city)
    env = env.merge(urban[["city", "urbanization_index"]], on="city", how="left")
    city = city.merge(env, on="city", how="left", validate="one_to_one")

    city_predictors = [
        "mean_non_amg_phage_abundance",
        "mean_non_amg_phage_richness",
        "urbanization_index",
        "built_up_km2",
        "population_density_per_km2",
        "green_space_percent",
        "road_density_m_per_km2",
    ]
    tests = city_test_rows(
        city,
        "amg_prevalence",
        city_predictors,
        "Primary occurrence component",
        "Seven city-level occurrence tests",
        rng,
    )

    intensity_city = city.loc[city["positive_libraries"].fillna(0) >= 3].copy()
    tests.extend(
        city_test_rows(
            intensity_city,
            "positive_mean_log1p_amg",
            [
                "positive_mean_non_amg_phage_abundance",
                "positive_mean_non_amg_phage_richness",
            ],
            "Positive-library intensity sensitivity",
            "Two city-level positive-intensity tests",
            rng,
        )
    )

    library_rows: list[dict[str, object]] = []
    library_specs = [
        (
            sample,
            "amg_detected",
            "Library occurrence component",
            "Two library-level occurrence tests",
        ),
        (
            positive.groupby("city").filter(lambda x: len(x) >= 2),
            "total_amg",
            "Positive-library intensity sensitivity",
            "Two library-level positive-intensity tests",
        ),
    ]
    for frame, outcome, role, family in library_specs:
        for predictor in ["non_amg_phage_abundance", "non_amg_phage_richness"]:
            rho, p_value = spearman_permutation(
                frame[outcome].to_numpy(),
                frame[predictor].to_numpy(),
                rng,
                strata=frame["city"].to_numpy(),
            )
            library_rows.append(
                {
                    "analysis_level": "Library",
                    "analysis_role": role,
                    "outcome": outcome,
                    "predictor": predictor,
                    "n": len(frame),
                    "method": f"Ranked variables centered within city; {PERMUTATIONS:,} within-city permutations",
                    "effect": rho,
                    "effect_name": "Within-city rank correlation",
                    "ci95_low": np.nan,
                    "ci95_high": np.nan,
                    "permutation_p": p_value,
                    "correction_family": family,
                }
            )
    tests.extend(library_rows)
    tests_df = pd.DataFrame(tests)
    tests_df["bh_q"] = np.nan
    for family, idx in tests_df.groupby("correction_family").groups.items():
        tests_df.loc[idx, "bh_q"] = bh_adjust(tests_df.loc[idx, "permutation_p"].tolist())

    gz = sample.loc[sample["city"] == "Guangzhou"].copy()
    sh = sample.loc[sample["city"] == "Shanghai"].copy()
    if len(gz) != 40 or len(sh) != 40:
        raise ValueError("Expected 40 libraries in Guangzhou and Shanghai.")

    odds_ratio, fisher_p = fisher_exact_two_sided(
        int(gz["amg_detected"].sum()),
        int((1 - gz["amg_detected"]).sum()),
        int(sh["amg_detected"].sum()),
        int((1 - sh["amg_detected"]).sum()),
    )
    comparisons: list[dict[str, object]] = [
        {
            "comparison": "AMG detection prevalence",
            "test": "Two-sided Fisher exact test",
            "guangzhou_n": len(gz),
            "shanghai_n": len(sh),
            "guangzhou_value": float(gz["amg_detected"].mean()),
            "shanghai_value": float(sh["amg_detected"].mean()),
            "descriptive_statistic": "Proportion",
            "effect": odds_ratio,
            "effect_name": "Uncorrected sample odds ratio",
            "p_value": fisher_p,
        }
    ]
    for label, column, positive_only in [
        ("Positive-library total AMG abundance", "total_amg", True),
        ("Non-AMG phage abundance", "non_amg_phage_abundance", False),
        ("Non-AMG phage richness", "non_amg_phage_richness", False),
    ]:
        x = gz.loc[gz["amg_detected"] == 1, column].to_numpy() if positive_only else gz[column].to_numpy()
        y = sh.loc[sh["amg_detected"] == 1, column].to_numpy() if positive_only else sh[column].to_numpy()
        _, effect, p_value = rank_sum_permutation(x, y, rng)
        comparisons.append(
            {
                "comparison": label,
                "test": f"Two-sided permutation rank-sum ({PERMUTATIONS:,} permutations)",
                "guangzhou_n": len(x),
                "shanghai_n": len(y),
                "guangzhou_value": float(np.median(x)),
                "shanghai_value": float(np.median(y)),
                "descriptive_statistic": "Median",
                "guangzhou_mean": float(np.mean(x)),
                "shanghai_mean": float(np.mean(y)),
                "effect": effect,
                "effect_name": "Rank-biserial correlation; positive favors Guangzhou",
                "p_value": p_value,
            }
        )
    comparisons_df = pd.DataFrame(comparisons)
    comparisons_df["bh_q_across_four"] = bh_adjust(comparisons_df["p_value"].tolist())

    sample.to_csv(RESULTS / "sample_metrics_optimized.csv", index=False, encoding="utf-8-sig")
    city.to_csv(RESULTS / "city_summary_optimized.csv", index=False, encoding="utf-8-sig")
    tests_df.to_csv(RESULTS / "correlation_tests_optimized.csv", index=False, encoding="utf-8-sig")
    comparisons_df.to_csv(
        RESULTS / "guangzhou_shanghai_optimized.csv", index=False, encoding="utf-8-sig"
    )

    key = {
        "settings": {"seed": SEED, "permutations": PERMUTATIONS, "bootstraps": BOOTSTRAPS},
        "qc": qc,
        "guangzhou": {
            "positive": int(gz["amg_detected"].sum()),
            "n": len(gz),
            "positive_median_amg": float(gz.loc[gz["amg_detected"] == 1, "total_amg"].median()),
            "median_non_amg_richness": float(gz["non_amg_phage_richness"].median()),
        },
        "shanghai": {
            "positive": int(sh["amg_detected"].sum()),
            "n": len(sh),
            "positive_median_amg": float(sh.loc[sh["amg_detected"] == 1, "total_amg"].median()),
            "median_non_amg_richness": float(sh["non_amg_phage_richness"].median()),
        },
    }
    (RESULTS / "key_results_optimized.json").write_text(
        json.dumps(key, ensure_ascii=False, indent=2), encoding="utf-8"
    )



if __name__ == "__main__":
    main()
