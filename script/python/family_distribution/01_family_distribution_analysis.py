                                                                        

                                                                            
                                                                                
   

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "inputs"
OUTPUT = ROOT / "results"

RPKM_PATH = INPUT / "RPKM_veter_matrix_clean.tsv"
TAXONOMY_PATH = INPUT / "vOTU_taxonomy.xlsx"
METADATA_PATH = INPUT / "Supplementary Data 1.xlsx"


def key(value: object) -> str:
    return str(value).strip().casefold()


def concentration_indices(values: pd.Series) -> dict[str, float]:
                                                                       
    x = values.to_numpy(dtype=float)
    total = float(x.sum())
    n = len(x)
    if total <= 0 or n < 2:
        return {
            "simpson": np.nan,
            "berger_parker": np.nan,
            "levins_ba": np.nan,
            "top_share": np.nan,
        }
    p = x / total
    d = float(np.square(p).sum())
    return {
        "simpson": d,
        "berger_parker": float(p.max()),
        "levins_ba": float(((1.0 / d) - 1.0) / (n - 1.0)),
        "top_share": float(p.max()),
    }


def distribution_class(row: pd.Series) -> str:
                                                                          
    if (
        row["city_occupancy_fraction"] >= 0.75
        and row["city_standardized_levins_breadth"] >= 0.50
        and row["order_occupancy_fraction"] >= 0.50
    ):
        return "Pan-urban ubiquitous within surveyed network"
    if (
        row["city_occupancy_fraction"] <= 0.25
        or row["city_berger_parker_dominance"] >= 0.50
    ):
        return "Geographically restricted"
    if row["city_occupancy_fraction"] >= 0.75:
        return "Widespread but uneven"
    return "Intermediate"


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)

    taxonomy = pd.read_excel(TAXONOMY_PATH)
    required_taxonomy = {"seq_name", "taxonomy", "classify"}
    missing = required_taxonomy.difference(taxonomy.columns)
    if missing:
        raise ValueError(f"Missing taxonomy columns: {sorted(missing)}")
    taxonomy = taxonomy.loc[taxonomy["classify"].eq("veter")].copy()
    taxonomy["seq_name"] = taxonomy["seq_name"].astype(str)
    taxonomy["family"] = taxonomy["taxonomy"].astype(str).str.strip()
    if len(taxonomy) != 17_439 or taxonomy["seq_name"].nunique() != 17_439:
        raise ValueError("Expected 17,439 unique vertebrate-associated vOTUs")
    if taxonomy["family"].nunique() != 32:
        raise ValueError("Expected 32 vertebrate-associated viral families")

    rpkm = pd.read_csv(RPKM_PATH, sep="\t", index_col=0)
    rpkm.index = rpkm.index.astype(str)
    if rpkm.shape != (17_439, 985):
        raise ValueError(f"Unexpected RPKM matrix shape: {rpkm.shape}")
    if rpkm.index.has_duplicates or rpkm.columns.has_duplicates:
        raise ValueError("RPKM identifiers must be unique")
    if np.any(rpkm.to_numpy(dtype=float) < 0):
        raise ValueError("RPKM values must be non-negative")

    family_map = taxonomy.set_index("seq_name")["family"]
    if set(rpkm.index) != set(family_map.index):
        raise ValueError("RPKM and taxonomy vOTU identifiers do not match")
    family_rpkm = rpkm.groupby(family_map.loc[rpkm.index].to_numpy()).sum().sort_index()
    family_rpkm.index.name = "family"

    library_totals = family_rpkm.sum(axis=0)
    library_ra = (
        family_rpkm.div(library_totals.replace(0, np.nan), axis=1)
        .fillna(0.0)
        .mul(100.0)
    )
    positive_sums = library_ra.sum(axis=0)[library_totals > 0]
    if not np.allclose(positive_sums.to_numpy(), 100.0, atol=1e-9):
        raise ValueError("Non-zero library profiles do not sum to 100%")

    metadata = pd.read_excel(METADATA_PATH, sheet_name=0, header=1)
    required_metadata = {
        "Library",
        "Location of Sampling (Sampling site, Province)",
        "Inferred host order",
    }
    missing = required_metadata.difference(metadata.columns)
    if missing:
        raise ValueError(f"Missing metadata columns: {sorted(missing)}")
    metadata = metadata.copy()
    metadata["library_key"] = metadata["Library"].map(key)
    metadata["city"] = (
        metadata["Location of Sampling (Sampling site, Province)"]
        .astype(str)
        .str.split(",", n=1)
        .str[0]
        .str.strip()
    )
    metadata["host_order"] = metadata["Inferred host order"].astype("string").str.strip()

    matrix_lookup = {key(column): column for column in family_rpkm.columns}
    if set(metadata["library_key"]) != set(matrix_lookup):
        raise ValueError("Metadata and RPKM library identifiers do not match")
    metadata["matrix_library"] = metadata["library_key"].map(matrix_lookup)
    metadata = metadata.set_index("matrix_library").loc[family_rpkm.columns].copy()
    metadata["matrix_library"] = metadata.index.astype(str)
    metadata = metadata.reset_index(drop=True)
    if len(metadata) != 985 or metadata["matrix_library"].nunique() != 985:
        raise ValueError("Metadata did not align one-to-one with 985 libraries")

    assigned = metadata.dropna(subset=["host_order"]).copy()
    city_groups = {
        name: group["matrix_library"].tolist()
        for name, group in metadata.groupby("city", sort=True)
    }
    order_groups = {
        name: group["matrix_library"].tolist()
        for name, group in assigned.groupby("host_order", sort=True)
    }
    if len(city_groups) != 22 or len(order_groups) != 22:
        raise ValueError("Expected 22 cities and 22 assigned avian orders")
    well_sampled_order_groups = {
        name: libraries for name, libraries in order_groups.items() if len(libraries) >= 10
    }
    if len(well_sampled_order_groups) != 7:
        raise ValueError("Expected seven avian orders represented by at least 10 libraries")

    def mean_profile(groups: dict[str, list[str]]) -> pd.DataFrame:
        return pd.DataFrame(
            {
                name: library_ra[libraries].mean(axis=1)
                for name, libraries in groups.items()
            }
        ).sort_index()

    def detected(groups: dict[str, list[str]]) -> pd.DataFrame:
        return pd.DataFrame(
            {
                name: (family_rpkm[libraries] > 0).any(axis=1)
                for name, libraries in groups.items()
            }
        ).sort_index()

    city_mean = mean_profile(city_groups)
    order_mean = mean_profile(order_groups)
    city_detect = detected(city_groups)
    order_detect = detected(order_groups)
    well_order_detect = detected(well_sampled_order_groups)

    records: list[dict[str, object]] = []
    for family in family_rpkm.index:
        city_index = concentration_indices(city_mean.loc[family])
        order_index = concentration_indices(order_mean.loc[family])
        records.append(
            {
                "family": family,
                "votu_count": int((family_map == family).sum()),
                "libraries_detected_rpkm_gt_0": int((family_rpkm.loc[family] > 0).sum()),
                "library_prevalence_fraction": float((family_rpkm.loc[family] > 0).mean()),
                "mean_library_relative_abundance_percent": float(library_ra.loc[family].mean()),
                "cities_detected_n_of_22": int(city_detect.loc[family].sum()),
                "city_occupancy_fraction": float(city_detect.loc[family].mean()),
                "city_simpson_concentration": city_index["simpson"],
                "city_berger_parker_dominance": city_index["berger_parker"],
                "city_standardized_levins_breadth": city_index["levins_ba"],
                "top_city": str(city_mean.loc[family].idxmax()),
                "top_city_share_of_family_profile_percent": float(city_index["top_share"] * 100.0),
                "avian_orders_detected_n_of_22": int(order_detect.loc[family].sum()),
                "order_occupancy_fraction": float(order_detect.loc[family].mean()),
                "order_simpson_concentration": order_index["simpson"],
                "order_berger_parker_dominance": order_index["berger_parker"],
                "order_standardized_levins_breadth": order_index["levins_ba"],
                "top_avian_order": str(order_mean.loc[family].idxmax()),
                "top_order_share_of_family_profile_percent": float(order_index["top_share"] * 100.0),
                "well_sampled_orders_detected_n_of_7": int(well_order_detect.loc[family].sum()),
            }
        )

    metrics = pd.DataFrame(records)
    metrics["distribution_class"] = metrics.apply(distribution_class, axis=1)
    class_order = {
        "Pan-urban ubiquitous within surveyed network": 0,
        "Widespread but uneven": 1,
        "Intermediate": 2,
        "Geographically restricted": 3,
    }
    metrics["_class_order"] = metrics["distribution_class"].map(class_order)
    metrics = (
        metrics.sort_values(
            ["_class_order", "mean_library_relative_abundance_percent"],
            ascending=[True, False],
        )
        .drop(columns="_class_order")
        .reset_index(drop=True)
    )

    metrics.to_csv(
        OUTPUT / "family_indices.csv", index=False, encoding="utf-8-sig"
    )
    city_mean.rename_axis("family").reset_index().to_csv(
        OUTPUT / "city_mean_relative_abundance_percent.csv",
        index=False,
        encoding="utf-8-sig",
    )
    order_mean.rename_axis("family").reset_index().to_csv(
        OUTPUT / "avian_order_mean_relative_abundance_percent.csv",
        index=False,
        encoding="utf-8-sig",
    )
    pd.DataFrame(
        {
            "city": list(city_groups),
            "library_count": [len(x) for x in city_groups.values()],
        }
    ).to_csv(OUTPUT / "city_library_counts.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(
        {
            "avian_order": list(order_groups),
            "library_count": [len(x) for x in order_groups.values()],
        }
    ).to_csv(
        OUTPUT / "avian_order_library_counts.csv", index=False, encoding="utf-8-sig"
    )


if __name__ == "__main__":
    main()
