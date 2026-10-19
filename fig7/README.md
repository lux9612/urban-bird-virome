# Figure 7

AMG detection rates and abundance among positive libraries.

## Requirements

Python 3.12 or newer.

```bash
pip install -r requirements.txt
python plot_fig7b.py
```

On Windows, `run_windows.bat` creates a local virtual environment, installs the
required packages, and runs the script.

## Files

| Path | Contents |
| ---- | -------- |
| `plot_fig7b.py` | Figure 7B plotting script |
| [Figure7C_Expanded.pdf](Figure7C_Expanded.pdf) | Expanded Figure 7C |
| `data/library_metrics.csv` | AMG abundance for 985 libraries |
| `data/city_summary.csv` | Detection counts, rates, and Wilson 95% confidence intervals for 22 cities |
| `data/statistics.csv` | Overall statistical test results |
| `data/positive_abundance_dunn.csv` | Dunn pairwise comparisons with Benjamini-Hochberg correction |

Outputs: `output/Figure7B.pdf`, `output/Figure7B.svg`, and `output/Figure7B.png`.
