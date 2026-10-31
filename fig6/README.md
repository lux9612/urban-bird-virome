# Figure 6

AMG detection rates and abundance among positive libraries.

## Requirements

Python 3.12 or newer.

```bash
pip install -r requirements.txt
python plot_fig6b.py
```

On Windows, `run_windows.bat` creates a local virtual environment, installs the
required packages, and runs the script.

## Files

| Path | Contents |
| ---- | -------- |
| `plot_fig6b.py` | Figure 6B plotting script |
| `data/library_metrics.csv` | AMG abundance for 985 libraries |
| `data/city_summary.csv` | Detection counts, rates, and Wilson 95% confidence intervals for 22 cities |
| `data/statistics.csv` | Overall statistical test results |
| `data/positive_abundance_dunn.csv` | Dunn pairwise comparisons with Benjamini-Hochberg correction |

The script checks the library count, sample identifiers, positive-library count,
and agreement between the library-level and city-level data. Outputs are written
to `output/Figure6B.pdf`, `output/Figure6B.svg`, and
`output/Figure6B.png`. The output directory is not tracked by Git.

The expanded Figure 6C is provided as
[`Figure6C_Expanded.pdf`](Figure6C_Expanded.pdf).
