                                                                               

from __future__ import annotations

import hashlib
import heapq
import math
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
RESULTS = HERE / "results"
MATRIX = INPUTS / "RPKM_veter_matrix_clean.tsv"
METADATA = INPUTS / "Supplementary Data 1.xlsx"
SEED = 123
PERMUTATIONS = 5000
PRIMARY_JACCARD = 0.06
MIN_SHARED = 30
ALPHA = 0.01


def clean_city(value: object) -> str:
    return str(value).split(",", 1)[0].strip()


def normalize_id(value: object) -> str:
    return str(value).strip().casefold()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def benjamini_hochberg(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty(len(values), dtype=float)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def load_metadata_and_header() -> tuple[list[str], np.ndarray, list[str], pd.DataFrame]:
    metadata = pd.read_excel(METADATA, sheet_name=0, header=1)
    required = {"Library", "Location of Sampling (Sampling site, Province)"}
    if not required.issubset(metadata.columns):
        raise ValueError(f"Missing metadata columns: {sorted(required - set(metadata.columns))}")
    metadata = metadata.copy()
    metadata["library_key"] = metadata["Library"].map(normalize_id)
    metadata["city"] = metadata["Location of Sampling (Sampling site, Province)"].map(clean_city)
    if metadata["library_key"].duplicated().any():
        raise ValueError("Duplicate normalized Library identifiers")
    if metadata["city"].nunique() != 22:
        raise ValueError(f"Expected 22 cities, found {metadata['city'].nunique()}")

    header = pd.read_csv(MATRIX, sep="\t", nrows=0).columns.astype(str).tolist()
    if header[0] != "Contig" or len(header) != 986:
        raise ValueError(f"Unexpected matrix header: first={header[0]!r}, columns={len(header)}")
    samples = header[1:]
    metadata_lookup = dict(zip(metadata["library_key"], metadata["city"]))
    missing = [sample for sample in samples if normalize_id(sample) not in metadata_lookup]
    if missing:
        raise ValueError(f"Matrix libraries absent from metadata: {missing[:10]}")
    if set(map(normalize_id, samples)) != set(metadata["library_key"]):
        raise ValueError("The 985 matrix and metadata library sets are not identical")
    city_names = sorted(metadata["city"].unique().tolist())
    city_to_index = {city: index for index, city in enumerate(city_names)}
    sample_city_indices = np.array(
        [city_to_index[metadata_lookup[normalize_id(sample)]] for sample in samples], dtype=np.int16
    )
    return samples, sample_city_indices, city_names, metadata


def build_city_presence(samples: list[str], sample_city_indices: np.ndarray, city_names: list[str]) -> tuple[list[str], np.ndarray]:
    identifiers: list[str] = []
    blocks: list[np.ndarray] = []
    for chunk in pd.read_csv(MATRIX, sep="\t", chunksize=500):
        identifiers.extend(chunk.iloc[:, 0].astype(str).tolist())
        values = chunk.iloc[:, 1:].to_numpy(dtype=float, copy=False)
        if np.any(values < 0):
            raise ValueError("Negative RPKM values are invalid for presence-absence analysis")
        presence = np.zeros((len(chunk), len(city_names)), dtype=bool)
        for city_index in range(len(city_names)):
            presence[:, city_index] = np.any(values[:, sample_city_indices == city_index] > 0, axis=1)
        blocks.append(presence)
    matrix = np.vstack(blocks)
    if matrix.shape != (17_439, 22):
        raise ValueError(f"Expected a 17,439 × 22 city-presence matrix, found {matrix.shape}")
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("Duplicate vOTU identifiers in all-virus RPKM matrix")
    return identifiers, matrix


def pairwise_network(city_names: list[str], presence: np.ndarray) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    universe = presence.shape[0]
    richness = presence.sum(axis=0).astype(int)
    rows = []
    for first in range(len(city_names)):
        for second in range(first + 1, len(city_names)):
            shared = int(np.sum(presence[:, first] & presence[:, second]))
            union = int(richness[first] + richness[second] - shared)
            jaccard = shared / union if union else 0.0
            null = rng.hypergeometric(
                ngood=int(richness[first]),
                nbad=int(universe - richness[first]),
                nsample=int(richness[second]),
                size=PERMUTATIONS,
            )
            empirical_p = (int(np.sum(null >= shared)) + 1.0) / (PERMUTATIONS + 1.0)
            rows.append(
                {
                    "city_1": city_names[first],
                    "city_2": city_names[second],
                    "richness_city_1": int(richness[first]),
                    "richness_city_2": int(richness[second]),
                    "shared_votu_count": shared,
                    "union_votu_count": union,
                    "jaccard_similarity": jaccard,
                    "permutation_p": empirical_p,
                    "null_mean_shared": float(null.mean()),
                    "null_q99_shared": float(np.quantile(null, 0.99, method="higher")),
                }
            )
    frame = pd.DataFrame(rows)
    frame["permutation_q_bh"] = benjamini_hochberg(frame["permutation_p"].to_numpy(float))
    frame["primary_edge"] = (
        (frame["permutation_p"] < ALPHA)
        & (frame["jaccard_similarity"] >= PRIMARY_JACCARD)
        & (frame["shared_votu_count"] >= MIN_SHARED)
    )
    return frame


def adjacency_from_edges(city_names: list[str], edges: pd.DataFrame, distance_mode: str = "inverse") -> dict[str, dict[str, float]]:
    adjacency = {city: {} for city in city_names}
    for row in edges.itertuples(index=False):
        similarity = float(row.jaccard_similarity)
        if distance_mode == "inverse":
            distance = 1.0 / similarity
        elif distance_mode == "one_minus":
            distance = 1.0 - similarity
        elif distance_mode == "unweighted":
            distance = 1.0
        else:
            raise ValueError(distance_mode)
        adjacency[row.city_1][row.city_2] = distance
        adjacency[row.city_2][row.city_1] = distance
    return adjacency


def connected_components(adjacency: dict[str, dict[str, float]]) -> list[list[str]]:
    remaining = set(adjacency)
    components = []
    while remaining:
        seed = min(remaining)
        stack = [seed]
        component = []
        remaining.remove(seed)
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbor in adjacency[node]:
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    stack.append(neighbor)
        components.append(sorted(component))
    return sorted(components, key=lambda group: (-len(group), group))


def dijkstra(adjacency: dict[str, dict[str, float]], source: str) -> dict[str, float]:
    distances = {node: math.inf for node in adjacency}
    distances[source] = 0.0
    queue = [(0.0, source)]
    while queue:
        distance, node = heapq.heappop(queue)
        if distance > distances[node] + 1e-12:
            continue
        for neighbor, length in adjacency[node].items():
            candidate = distance + length
            if candidate + 1e-12 < distances[neighbor]:
                distances[neighbor] = candidate
                heapq.heappush(queue, (candidate, neighbor))
    return distances


def betweenness(adjacency: dict[str, dict[str, float]]) -> dict[str, float]:
    nodes = sorted(adjacency)
    centrality = {node: 0.0 for node in nodes}
    for source in nodes:
        stack: list[str] = []
        predecessors = {node: [] for node in nodes}
        paths = {node: 0.0 for node in nodes}
        paths[source] = 1.0
        distances = {node: math.inf for node in nodes}
        distances[source] = 0.0
        queue = [(0.0, source)]
        while queue:
            distance, node = heapq.heappop(queue)
            if distance > distances[node] + 1e-12:
                continue
            stack.append(node)
            for neighbor, length in adjacency[node].items():
                candidate = distance + length
                if candidate + 1e-12 < distances[neighbor]:
                    distances[neighbor] = candidate
                    heapq.heappush(queue, (candidate, neighbor))
                    paths[neighbor] = paths[node]
                    predecessors[neighbor] = [node]
                elif abs(candidate - distances[neighbor]) <= 1e-12:
                    paths[neighbor] += paths[node]
                    predecessors[neighbor].append(node)
        dependency = {node: 0.0 for node in nodes}
        while stack:
            target = stack.pop()
            if paths[target] > 0:
                coefficient = (1.0 + dependency[target]) / paths[target]
                for predecessor in predecessors[target]:
                    dependency[predecessor] += paths[predecessor] * coefficient
            if target != source:
                centrality[target] += dependency[target]
    normalization = 1.0 / ((len(nodes) - 1) * (len(nodes) - 2)) if len(nodes) > 2 else 0.0
    return {node: value * normalization for node, value in centrality.items()}


def centralities(city_names: list[str], edges: pd.DataFrame, distance_mode: str = "inverse") -> pd.DataFrame:
    adjacency = adjacency_from_edges(city_names, edges, distance_mode)
    n = len(city_names)
    between = betweenness(adjacency)
    rows = []
    edge_lookup = {}
    for row in edges.itertuples(index=False):
        edge_lookup[frozenset((row.city_1, row.city_2))] = float(row.jaccard_similarity)
    for city in city_names:
        distances = dijkstra(adjacency, city)
        finite = [value for node, value in distances.items() if node != city and math.isfinite(value)]
        reachable = len(finite)
        closeness = (reachable / (n - 1)) * (reachable / sum(finite)) if finite and sum(finite) > 0 else 0.0
        degree = len(adjacency[city])
        strength = sum(edge_lookup[frozenset((city, neighbor))] for neighbor in adjacency[city])
        rows.append(
            {
                "city": city,
                "degree": degree,
                "degree_centrality": degree / (n - 1),
                "jaccard_strength": strength,
                "normalized_jaccard_strength": strength / (n - 1),
                "betweenness_centrality": between[city],
                "closeness_centrality": closeness,
                "reachable_nodes": reachable,
            }
        )
    return pd.DataFrame(rows)


def classify_nodes(frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    metrics = ["degree_centrality", "betweenness_centrality", "closeness_centrality"]
    high_counts = np.zeros(len(output), dtype=int)
    low_counts = np.zeros(len(output), dtype=int)
    z_columns = []
    for metric in metrics:
        high = float(output[metric].quantile(0.75))
        low = float(output[metric].quantile(0.25))
        output[f"{metric}_q75"] = high
        output[f"{metric}_q25"] = low
        high_counts += (
            (output[metric].to_numpy(float) >= high) & (output[metric].to_numpy(float) > 0)
        ).astype(int)
        low_counts += (output[metric].to_numpy(float) <= low).astype(int)
        standard_deviation = float(output[metric].std(ddof=0))
        z_name = f"z_{metric}"
        output[z_name] = (
            (output[metric] - float(output[metric].mean())) / standard_deviation
            if standard_deviation > 0
            else 0.0
        )
        z_columns.append(z_name)
    output["high_metric_count"] = high_counts
    output["low_metric_count"] = low_counts
    output["composite_centrality_z"] = output[z_columns].mean(axis=1)
    output["network_role"] = np.select(
        [output["high_metric_count"] >= 2, output["low_metric_count"] >= 2],
        ["hub", "peripheral"],
        default="intermediate",
    )
    return output.sort_values(
        ["network_role", "composite_centrality_z"],
        key=lambda series: series.map({"hub": 0, "intermediate": 1, "peripheral": 2})
        if series.name == "network_role"
        else -series,
    )


def analyse_variant(city_names: list[str], pairs: pd.DataFrame, name: str, jaccard: float, use_q: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    significance = pairs["permutation_q_bh"] < 0.05 if use_q else pairs["permutation_p"] < ALPHA
    mask = significance & (pairs["jaccard_similarity"] >= jaccard) & (pairs["shared_votu_count"] >= MIN_SHARED)
    edges = pairs.loc[mask].copy()
    metrics = classify_nodes(centralities(city_names, edges, "inverse"))
    metrics["variant"] = name
    metrics["edge_count"] = len(edges)
    metrics["component_count"] = len(connected_components(adjacency_from_edges(city_names, edges)))
    return edges, metrics


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    samples, sample_city_indices, city_names, metadata = load_metadata_and_header()
    identifiers, presence = build_city_presence(samples, sample_city_indices, city_names)

    presence_frame = pd.DataFrame(presence.astype(np.uint8), columns=city_names)
    presence_frame.insert(0, "votu_id", identifiers)
    presence_frame.to_csv(RESULTS / "city_votu_presence.tsv", sep="\t", index=False)

    pairs = pairwise_network(city_names, presence)
    pairs.to_csv(RESULTS / "all_city_pair_network_statistics.csv", index=False, encoding="utf-8-sig")
    primary_edges = pairs.loc[pairs["primary_edge"]].copy()
    primary_edges.to_csv(RESULTS / "primary_network_edges.csv", index=False, encoding="utf-8-sig")

    primary_metrics = classify_nodes(centralities(city_names, primary_edges, "inverse"))
    one_minus = classify_nodes(centralities(city_names, primary_edges, "one_minus"))
    unweighted = classify_nodes(centralities(city_names, primary_edges, "unweighted"))
    sensitivity_fields = [
        "city",
        "betweenness_centrality",
        "closeness_centrality",
        "composite_centrality_z",
        "network_role",
    ]
    one_minus = one_minus[sensitivity_fields].rename(
        columns={column: f"one_minus_{column}" for column in sensitivity_fields if column != "city"}
    )
    unweighted = unweighted[sensitivity_fields].rename(
        columns={column: f"unweighted_{column}" for column in sensitivity_fields if column != "city"}
    )
    primary_metrics = primary_metrics.merge(one_minus, on="city", validate="one_to_one").merge(
        unweighted, on="city", validate="one_to_one"
    )

    total_shared = []
    for city in city_names:
        connected = primary_edges.loc[(primary_edges["city_1"] == city) | (primary_edges["city_2"] == city)]
        total_shared.append(
            {
                "city": city,
                "sum_pairwise_shared_votu_counts": int(connected["shared_votu_count"].sum()),
                "union_shared_votu_count": int(
                    np.sum(
                        presence[:, city_names.index(city)]
                        & np.any(
                            presence[:, [city_names.index(other) for other in city_names if other != city]], axis=1
                        )
                    )
                ),
            }
        )
    primary_metrics = primary_metrics.merge(pd.DataFrame(total_shared), on="city", validate="one_to_one")
    primary_metrics.to_csv(RESULTS / "city_network_centrality_and_roles.csv", index=False, encoding="utf-8-sig")

    variants = []
    variant_definitions = [
        ("jaccard_0.05", 0.05, False),
        ("primary_jaccard_0.06", 0.06, False),
        ("jaccard_0.07", 0.07, False),
        ("bh_q_0.05_jaccard_0.06", 0.06, True),
    ]
    for name, threshold, use_q in variant_definitions:
        _, metrics = analyse_variant(city_names, pairs, name, threshold, use_q)
        variants.append(metrics)
    variant_frame = pd.concat(variants, ignore_index=True)
    variant_frame.to_csv(RESULTS / "network_threshold_sensitivity.csv", index=False, encoding="utf-8-sig")
    stability = (
        variant_frame.assign(is_hub=variant_frame["network_role"].eq("hub"))
        .groupby("city", as_index=False)
        .agg(hub_variants=("is_hub", "sum"), variants_tested=("variant", "nunique"))
    )
    stability["robust_hub"] = stability["hub_variants"] >= 3
    stability.to_csv(RESULTS / "hub_classification_stability.csv", index=False, encoding="utf-8-sig")

    city_counts = metadata["city"].value_counts().sort_index().rename("library_count").rename_axis("city").reset_index()
    city_richness = pd.DataFrame({"city": city_names, "city_votu_richness": presence.sum(axis=0).astype(int)})
    city_counts.merge(city_richness, on="city", validate="one_to_one").to_csv(
        RESULTS / "city_library_and_votu_counts.csv", index=False, encoding="utf-8-sig"
    )


if __name__ == "__main__":
    main()
