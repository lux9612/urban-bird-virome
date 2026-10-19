#!/usr/bin/env python3
                                                                               

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


def read_ucdb(path: Path, columns: list[str]) -> pd.DataFrame:
    last_error = None
    for encoding in ("utf-8-sig", "cp1252", "latin1"):
        try:
            frame = pd.read_csv(
                path,
                encoding=encoding,
                usecols=lambda name: name.replace("ï»¿", "").lstrip("\ufeff") in columns,
            )
            break
        except UnicodeDecodeError as error:
            last_error = error
    else:
        raise last_error
    frame.columns = [name.replace("ï»¿", "") for name in frame.columns]
    missing = set(columns).difference(frame.columns)
    if missing:
        raise ValueError(f"{path.name} is missing columns: {sorted(missing)}")
    return frame[columns]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    here = Path(__file__).resolve().parent
    source = here / "official_source_csv"
    mapping_path = here / "city_ucdb_mapping.csv"
    files = {
        "general": next(source.glob("*GENERAL_CHARACTERISTICS*.csv")),
        "ghsl": next(source.glob("*GHSL_GLOBE*.csv")),
        "greenness": next(source.glob("*GREENNESS*.csv")),
        "infrastructure": next(source.glob("*INFRASTRUCTURES*.csv")),
    }

    mapping = pd.read_csv(mapping_path)
    if len(mapping) != 22 or mapping["city"].nunique() != 22:
        raise ValueError("City mapping must contain the 22 unique sampled cities")
    included = mapping[mapping["include_in_environment_analysis"].eq("Yes")].copy()
    excluded = mapping[~mapping["include_in_environment_analysis"].eq("Yes")].copy()
    if len(included) != 21 or excluded["city"].tolist() != ["Nyingchi"]:
        raise ValueError("Expected 21 included cities and the documented Nyingchi exclusion")
    included["ID_UC_G0"] = included["ID_UC_G0"].astype(int)
    taian_mapping = included.loc[included["city"].eq("Tai'an"), "ID_UC_G0"]
    if len(taian_mapping) != 1 or int(taian_mapping.iloc[0]) != 9731:
        raise ValueError("Tai'an must map to official GHS entity ID 9731")

    general = read_ucdb(
        files["general"],
        [
            "ID_UC_G0",
            "GC_UCN_MAI_2025",
            "GC_UCN_LIS_2025",
            "GC_CNT_GAD_2025",
            "GC_UCA_KM2_2025",
            "GC_POP_TOT_2025",
        ],
    )
    ghsl = read_ucdb(
        files["ghsl"],
        ["ID_UC_G0", "GC_UCN_MAI_2025", "GC_CNT_GAD_2025", "GH_BUS_TOT_2025", "GH_XST_D30_2025"],
    ).rename(columns={"GC_UCN_MAI_2025": "GHSL_name", "GC_CNT_GAD_2025": "GHSL_country"})
    greenness = read_ucdb(
        files["greenness"],
        ["ID_UC_G0", "GC_UCN_MAI_2025", "GC_CNT_GAD_2025", "GR_SHB_HGR_2025"],
    ).rename(columns={"GC_UCN_MAI_2025": "greenness_name", "GC_CNT_GAD_2025": "greenness_country"})
    infrastructure = read_ucdb(
        files["infrastructure"],
        ["ID_UC_G0", "GC_UCN_MAI_2025", "GC_CNT_GAD_2025", "IN_ROA_DEN_2024"],
    ).rename(columns={"GC_UCN_MAI_2025": "infrastructure_name", "GC_CNT_GAD_2025": "infrastructure_country"})

    merged = included.merge(general, on="ID_UC_G0", how="left", validate="one_to_one")
    merged = merged.merge(ghsl, on="ID_UC_G0", how="left", validate="one_to_one")
    merged = merged.merge(greenness, on="ID_UC_G0", how="left", validate="one_to_one")
    merged = merged.merge(infrastructure, on="ID_UC_G0", how="left", validate="one_to_one")
    if merged.isna().any().any():
        columns = merged.columns[merged.isna().any()].tolist()
        raise ValueError(f"Missing values after exact ID joins: {columns}")
    country_columns = ["GC_CNT_GAD_2025", "GHSL_country", "greenness_country", "infrastructure_country"]
    if not merged[country_columns].eq("China").all().all():
        raise ValueError("At least one mapped urban centre is not classified as China across source themes")
    name_columns = ["GC_UCN_MAI_2025", "GHSL_name", "greenness_name", "infrastructure_name"]
    if not merged[name_columns].nunique(axis=1).eq(1).all():
        raise ValueError("GHS main names are inconsistent across source themes")

    result = pd.DataFrame(
        {
            "city": merged["city"],
            "ID_UC_G0": merged["ID_UC_G0"],
            "GHS_main_name": merged["GC_UCN_MAI_2025"],
            "GHS_name_list": merged["GC_UCN_LIS_2025"],
            "urban_centre_area_km2": pd.to_numeric(merged["GC_UCA_KM2_2025"], errors="raise"),
            "urban_centre_population_2025": pd.to_numeric(merged["GC_POP_TOT_2025"], errors="raise"),
            "built_up_km2": pd.to_numeric(merged["GH_BUS_TOT_2025"], errors="raise") / 1_000_000.0,
            "population_density_per_km2": pd.to_numeric(merged["GH_XST_D30_2025"], errors="raise"),
            "green_space_percent": pd.to_numeric(merged["GR_SHB_HGR_2025"], errors="raise"),
            "road_density_m_per_km2": pd.to_numeric(merged["IN_ROA_DEN_2024"], errors="raise") * 1_000_000.0,
            "matching_evidence": merged["matching_evidence"],
        }
    ).sort_values("city")

    if (result[["built_up_km2", "population_density_per_km2", "road_density_m_per_km2"]] <= 0).any().any():
        raise ValueError("A positive urbanization indicator is non-positive")
    if not result["green_space_percent"].between(0, 100).all():
        raise ValueError("Green-space percentage is outside [0,100]")
    taian = result.loc[result["city"].eq("Tai'an")]
    if len(taian) != 1 or int(taian.iloc[0]["ID_UC_G0"]) != 9731:
        raise ValueError("Extracted Tai'an record is not official GHS entity ID 9731")
    result.to_csv(here / "city_environment_indicators.csv", index=False, float_format="%.10g")
    excluded[["city", "include_in_environment_analysis", "matching_evidence"]].to_csv(
        here / "environment_exclusions.csv", index=False
    )

    checksum_lines = [f"{sha256(path)}  official_source_csv/{path.name}" for path in files.values()]
    (here / "source_checksums.sha256").write_text("\n".join(checksum_lines) + "\n", encoding="ascii")


if __name__ == "__main__":
    main()
