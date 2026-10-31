from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd

from newick_utils import distance_matrix, read_newick


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
RESULTS = HERE / "results"

TREE_PATH = INPUTS / "aves_order_7orders_all_original_tips.nwk"
BRAY_PATH = INPUTS / "BirdOrder_pairwise_BrayCurtis.tsv"

ORDERS = [
    "Accipitriformes",
    "Charadriiformes",
    "Columbiformes",
    "Falconiformes",
    "Galliformes",
    "Passeriformes",
    "Psittaciformes",
]
CHAR_VARIANTS = [f"Charadriiformes-{index}" for index in range(1, 9)]


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


def correlation(first: np.ndarray, second: np.ndarray, method: str) -> float:
    if method == "spearman":
        first = average_ranks(first)
        second = average_ranks(second)
    elif method != "pearson":
        raise ValueError(method)
    if np.std(first) == 0 or np.std(second) == 0:
        return float("nan")
    return float(np.corrcoef(first, second)[0, 1])


def upper_values(matrix: np.ndarray) -> np.ndarray:
    return matrix[np.triu_indices(matrix.shape[0], k=1)]


def exact_mantel(
    phylogenetic: np.ndarray,
    bray_curtis: np.ndarray,
    method: str,
) -> tuple[float, float, float, int]:
    observed = correlation(
        upper_values(phylogenetic), upper_values(bray_curtis), method
    )
    two_sided = 0
    positive = 0
    count = 0
    for permutation in itertools.permutations(range(phylogenetic.shape[0])):
        permuted = bray_curtis[np.ix_(permutation, permutation)]
        statistic = correlation(
            upper_values(phylogenetic), upper_values(permuted), method
        )
        two_sided += abs(statistic) >= abs(observed) - 1e-15
        positive += statistic >= observed - 1e-15
        count += 1
    return observed, two_sided / count, positive / count, count


def load_bray_curtis() -> tuple[pd.DataFrame, pd.DataFrame]:
    pairs = pd.read_csv(BRAY_PATH, sep="\t")
    required = {
        "Order1",
        "Order2",
        "BrayCurtis_dissimilarity",
        "BrayCurtis_similarity",
    }
    if set(pairs.columns) != required:
        raise ValueError(f"Unexpected Bray-Curtis columns: {pairs.columns.tolist()}")
    if len(pairs) != 21:
        raise ValueError(f"Expected 21 pairwise rows, found {len(pairs)}")
    if set(pairs["Order1"]) | set(pairs["Order2"]) != set(ORDERS):
        raise ValueError("The input does not contain exactly the seven target orders")
    canonical = pairs.apply(
        lambda row: tuple(sorted((row["Order1"], row["Order2"]))), axis=1
    )
    if canonical.nunique() != 21:
        raise ValueError("Duplicate or missing order pairs")
    if not np.allclose(
        pairs["BrayCurtis_dissimilarity"] + pairs["BrayCurtis_similarity"], 1.0
    ):
        raise ValueError("Bray-Curtis dissimilarity and similarity are inconsistent")
    if not (
        pairs["BrayCurtis_dissimilarity"].between(0, 1).all()
        and pairs["BrayCurtis_similarity"].between(0, 1).all()
    ):
        raise ValueError("Bray-Curtis values fall outside [0, 1]")
    matrix = pd.DataFrame(0.0, index=ORDERS, columns=ORDERS)
    for row in pairs.itertuples(index=False):
        matrix.loc[row.Order1, row.Order2] = row.BrayCurtis_dissimilarity
        matrix.loc[row.Order2, row.Order1] = row.BrayCurtis_dissimilarity
    return pairs, matrix


def phylogenetic_matrix(tree: object, char_tip: str) -> pd.DataFrame:
    names = [char_tip if order == "Charadriiformes" else order for order in ORDERS]
    values = np.round(distance_matrix(tree, names), 8)
    return pd.DataFrame(values, index=ORDERS, columns=ORDERS)


def pairwise_rows(
    phylogenetic: pd.DataFrame,
    pairs: pd.DataFrame,
    variant: str,
) -> list[dict[str, object]]:
    output = []
    for row in pairs.itertuples(index=False):
        output.append(
            {
                "charadriiformes_variant": variant,
                "order_1": row.Order1,
                "order_2": row.Order2,
                "phylogenetic_distance_Myr": phylogenetic.loc[
                    row.Order1, row.Order2
                ],
                "bray_curtis_dissimilarity": row.BrayCurtis_dissimilarity,
                "bray_curtis_similarity": row.BrayCurtis_similarity,
            }
        )
    return output


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    pairs, bray = load_bray_curtis()
    tree = read_newick(TREE_PATH)
    all_pairwise: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    leave_rows: list[dict[str, object]] = []
    for variant in CHAR_VARIANTS:
        phylogenetic = phylogenetic_matrix(tree, variant)
        rho, p_two, p_positive, permutations = exact_mantel(
            phylogenetic.to_numpy(), bray.to_numpy(), "spearman"
        )
        pearson_r, pearson_two, pearson_positive, _ = exact_mantel(
            phylogenetic.to_numpy(), bray.to_numpy(), "pearson"
        )
        summary_rows.append(
            {
                "charadriiformes_tip": variant,
                "avian_order_count": 7,
                "order_pair_count": 21,
                "mantel_spearman_rho": rho,
                "two_sided_exact_p": p_two,
                "directional_positive_exact_p": p_positive,
                "mantel_pearson_r": pearson_r,
                "pearson_two_sided_exact_p": pearson_two,
                "pearson_directional_positive_exact_p": pearson_positive,
                "label_permutations": permutations,
            }
        )
        all_pairwise.extend(pairwise_rows(phylogenetic, pairs, variant))
        for omitted in ORDERS:
            retained = [order for order in ORDERS if order != omitted]
            leave_rho, leave_two, leave_positive, leave_count = exact_mantel(
                phylogenetic.loc[retained, retained].to_numpy(),
                bray.loc[retained, retained].to_numpy(),
                "spearman",
            )
            leave_rows.append(
                {
                    "charadriiformes_tip": variant,
                    "omitted_order": omitted,
                    "retained_order_count": 6,
                    "retained_pair_count": 15,
                    "mantel_spearman_rho": leave_rho,
                    "two_sided_exact_p": leave_two,
                    "directional_positive_exact_p": leave_positive,
                    "label_permutations": leave_count,
                }
            )
    summary = pd.DataFrame(summary_rows)
    pairwise = pd.DataFrame(all_pairwise)
    leave_one_out = pd.DataFrame(leave_rows)
    summary.to_csv(
        RESULTS / "Charadriiformes_8variants_Mantel_results.tsv",
        sep="\t",
        index=False,
    )
    pairwise.to_csv(
        RESULTS / "all_variants_pairwise_data.tsv", sep="\t", index=False
    )
    leave_one_out.to_csv(
        RESULTS / "leave_one_order_out_sensitivity.tsv", sep="\t", index=False
    )
    if set(summary["label_permutations"]) != {5040}:
        raise RuntimeError("Seven-order exact enumeration did not produce 5,040 labels")
    if set(leave_one_out["label_permutations"]) != {720}:
        raise RuntimeError("Six-order exact enumeration did not produce 720 labels")


if __name__ == "__main__":
    main()
