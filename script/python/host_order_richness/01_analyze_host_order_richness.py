                                                                       

                                                                             
                                                                              
                                                                      
                                                                     
   

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
METADATA = INPUTS / "Supplementary Data 1.xlsx"
MATRIX = INPUTS / "RPKM_allvirus_matrix_clean.tsv"
RESULTS = HERE / "results"
SEED = 123
PERMUTATIONS = 9999


def normalize_id(value: object) -> str:
    return str(value).strip().casefold()


def expected_incidence_richness(prevalence: np.ndarray, n: int, m: int) -> float:
                                                                            
    prevalence = prevalence.astype(np.float64, copy=False)
    absence = np.ones(prevalence.shape, dtype=np.float64)
    for i in range(m):
        numerator = n - prevalence - i
        numerator[numerator < 0] = 0
        absence *= numerator / (n - i)
    return float(np.sum(1.0 - absence))


def spearman_permutation(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    xr = pd.Series(x).rank(method="average").to_numpy(float)
    yr = pd.Series(y).rank(method="average").to_numpy(float)
    observed = float(np.corrcoef(xr, yr)[0, 1])
    rng = np.random.default_rng(SEED)
    extreme = 0
    for _ in range(PERMUTATIONS):
        permuted = rng.permutation(yr)
        rho = float(np.corrcoef(xr, permuted)[0, 1])
        extreme += abs(rho) >= abs(observed) - 1e-15
    return observed, (extreme + 1) / (PERMUTATIONS + 1)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_excel(METADATA, header=1)
    metadata = metadata.rename(
        columns={
            "Library": "library",
            "Location of Sampling (Sampling site, Province)": "location",
            "Inferred host order": "host_order",
        }
    )
    metadata["city"] = metadata["location"].str.split(",", n=1).str[0].str.strip()
    metadata["library_key"] = metadata["library"].map(normalize_id)

    assert len(metadata) == 985, f"Expected 985 metadata rows, found {len(metadata)}"
    assert metadata["library_key"].is_unique, "Metadata library IDs are not unique"

    matrix_columns = pd.read_csv(MATRIX, sep="\t", nrows=0).columns.tolist()
    assert matrix_columns[0] == "Contig"
    matrix_lookup = {normalize_id(column): column for column in matrix_columns[1:]}
    assert len(matrix_lookup) == 985, f"Expected 985 matrix libraries, found {len(matrix_lookup)}"
    missing = sorted(set(metadata["library_key"]) - set(matrix_lookup))
    extra = sorted(set(matrix_lookup) - set(metadata["library_key"]))
    assert not missing and not extra, f"Metadata/matrix mismatch: missing={missing}, extra={extra}"
    metadata["matrix_column"] = metadata["library_key"].map(matrix_lookup)

    assigned = metadata.dropna(subset=["host_order"]).copy()
    group_samples: dict[tuple[str, str], list[str]] = {
        (city, order): group["matrix_column"].tolist()
        for (city, order), group in assigned.groupby(["city", "host_order"], sort=True)
    }
    richness = {key: 0 for key in group_samples}
    city_samples = {
        city: group["matrix_column"].tolist()
        for city, group in metadata.groupby("city", sort=True)
    }
    city_richness = {city: 0 for city in city_samples}

    passeriformes = assigned.loc[assigned["host_order"].eq("Passeriformes")].copy()
    passeriformes_groups = {
        city: group["matrix_column"].tolist()
        for city, group in passeriformes.groupby("city", sort=True)
    }
    rarefaction_depth = min(map(len, passeriformes_groups.values()))
    expected = {city: 0.0 for city in passeriformes_groups}

    use_columns = ["Contig", *matrix_columns[1:]]
    v_otu_count = 0
    for chunk in pd.read_csv(
        MATRIX,
        sep="\t",
        usecols=use_columns,
        dtype={column: np.float32 for column in matrix_columns[1:]},
        chunksize=2000,
    ):
        v_otu_count += len(chunk)
        values = chunk.drop(columns="Contig")
        detected = values.gt(0)
        for key, samples in group_samples.items():
            richness[key] += int(detected[samples].any(axis=1).sum())
        for city, samples in city_samples.items():
            city_richness[city] += int(detected[samples].any(axis=1).sum())
        for city, samples in passeriformes_groups.items():
            prevalence = detected[samples].sum(axis=1).to_numpy(np.int16)
            expected[city] += expected_incidence_richness(
                prevalence, len(samples), rarefaction_depth
            )

    assert v_otu_count == 52314, f"Expected 52,314 vOTUs, found {v_otu_count}"

    rows = []
    for (city, order), samples in group_samples.items():
        rows.append(
            {
                "city": city,
                "host_order": order,
                "library_count": len(samples),
                "observed_votu_richness": richness[(city, order)],
            }
        )
    city_order = pd.DataFrame(rows).sort_values(
        ["city", "observed_votu_richness", "host_order"], ascending=[True, False, True]
    )
    city_order.to_csv(RESULTS / "city_host_order_richness.csv", index=False)

    city_total = pd.DataFrame(
        {
            "city": list(city_samples),
            "library_count": [len(city_samples[city]) for city in city_samples],
            "observed_votu_richness": [city_richness[city] for city in city_samples],
        }
    ).sort_values("observed_votu_richness")
    city_total["richness_rank_low_to_high"] = np.arange(1, len(city_total) + 1)
    city_total.to_csv(RESULTS / "city_total_richness.csv", index=False)

    all_cities = sorted(metadata["city"].unique())
    passer_rows = []
    for city in all_cities:
        n = len(passeriformes_groups.get(city, []))
        passer_rows.append(
            {
                "city": city,
                "passeriformes_library_count": n,
                "observed_votu_richness": richness.get((city, "Passeriformes"), np.nan),
                f"expected_votu_richness_at_{rarefaction_depth}_libraries": expected.get(city, np.nan),
                "estimable": n > 0,
                "missing_reason": "" if n > 0 else "No library was assigned to Passeriformes",
            }
        )
    passer_summary = pd.DataFrame(passer_rows)
    passer_summary.to_csv(RESULTS / "passeriformes_city_summary.csv", index=False)

    yining = metadata.loc[metadata["city"].eq("Yining")]
    yining_composition = (
        yining.groupby("host_order", dropna=False)
        .size()
        .rename("library_count")
        .reset_index()
        .sort_values("library_count", ascending=False)
    )
    yining_composition["proportion_of_yining_libraries"] = (
        yining_composition["library_count"] / len(yining)
    )
    yining_composition.to_csv(RESULTS / "yining_host_order_composition.csv", index=False)

    represented = passer_summary.loc[passer_summary["estimable"]].copy()
    rho, p_value = spearman_permutation(
        represented["passeriformes_library_count"].to_numpy(float),
        represented["observed_votu_richness"].to_numpy(float),
    )
    stats = pd.DataFrame(
        [
            {
                "analysis": "Passeriformes sample count vs observed vOTU richness",
                "cities": len(represented),
                "statistic": "Spearman rho",
                "estimate": rho,
                "permutation_p": p_value,
                "permutations": PERMUTATIONS,
            }
        ]
    )
    stats.to_csv(RESULTS / "sampling_intensity_test.csv", index=False)


if __name__ == "__main__":
    main()
