"""Plot AMG detection rates and positive-library abundance (Figure 7B)."""
from pathlib import Path
P = Path(__file__).resolve().parent
(P / 'output').mkdir(exist_ok=True)
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.spines.top': False, 'axes.spines.right': False})
d = pd.read_csv(P / 'data/library_metrics.csv')
cs = pd.read_csv(P / 'data/city_summary.csv')
cities = cs.City.tolist()
x = np.arange(len(cities))
red = '#C00000'
blue = '#397C9A'
assert len(d) == 985 and int((d.Total_AMG > 0).sum()) == 266
assert len(cs) == 22 and int(cs.Positive.sum()) == 266
assert not d.Sample.duplicated().any()
for row in cs.itertuples():
    subset = d[d.Region == row.City]
    assert len(subset) == row.N and int((subset.Total_AMG > 0).sum()) == row.Positive
statistics = pd.read_csv(P / 'data/statistics.csv').set_index('Analysis')
posthoc = pd.read_csv(P / 'data/positive_abundance_dunn.csv')
detection_p = statistics.loc['AMG detection by city', 'P']
abundance_p = statistics.loc['AMG abundance among positives', 'P']
minimum_q = posthoc.BH_q.min()
fig, (ax, bx) = plt.subplots(2, 1, figsize=(5.0, 5.4), gridspec_kw={'height_ratios': [1, 1.25]})
fig.subplots_adjust(left=0.14, right=0.985, top=0.86, bottom=0.17, hspace=0.39)
fig.text(0.015, 0.965, '(B)', color=red, fontsize=11)
fig.text(0.14, 0.965, 'AMG detection and positive-library abundance', color=red, fontsize=8)
y = cs.Detection_rate.to_numpy() * 100
lo = cs.Wilson95_low.to_numpy() * 100
hi = cs.Wilson95_high.to_numpy() * 100
ax.bar(x, y, color=blue, width=0.68)
ax.errorbar(x, y, yerr=[y - lo, hi - y], fmt='none', ecolor='#222222', capsize=2, lw=0.7)
for i, r in cs.iterrows():
    ax.text(i, 112, f'{int(r.Positive)}/{int(r.N)}', ha='center', va='bottom', rotation=90, fontsize=6, color=red)
ax.set_ylim(0, 145)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_ylabel('AMG-positive libraries (%)', color=red)
ax.set_xticks([])
ax.set_xlim(-0.65, len(cities) - 0.35)
ax.set_title(f'{len(d)} libraries; chi-square P = {detection_p:.2e}', loc='left', color=red, fontsize=7, pad=24)
vals = [d.loc[(d.Region == c) & (d.Total_AMG > 0), 'Total_AMG'].to_numpy() for c in cities]
rng = np.random.default_rng(20261019)
for i, v in enumerate(vals):
    if len(v):
        if len(v) >= 2:
            bx.boxplot([v], positions=[i], widths=0.6, patch_artist=True, showfliers=False, boxprops={'facecolor': '#C3DAE4', 'edgecolor': '#333333', 'linewidth': 0.65}, medianprops={'color': '#111111', 'linewidth': 0.8}, whiskerprops={'linewidth': 0.6}, capprops={'linewidth': 0.6})
        bx.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=5, color='#222222', alpha=0.68, linewidths=0, zorder=3)
    else:
        bx.text(i, 2, 'No positives', rotation=90, ha='center', va='bottom', fontsize=6, color=red)
bx.set_yscale('log')
bx.set_ylim(0.9, 30000)
bx.set_xlim(-0.65, len(cities) - 0.35)
bx.set_xticks(x, cities, rotation=90)
bx.tick_params(axis='x', labelsize=7)
bx.set_ylabel('Total AMG abundance (RPKM)', color=red)
bx.set_title(f'{int((d.Total_AMG > 0).sum())} positive libraries; Kruskal–Wallis P = {abundance_p:.5f}\nDunn-BH: minimum q = {minimum_q:.4f}', loc='left', color=red, fontsize=7, pad=8)
for a in [ax, bx]:
    a.tick_params(width=0.6)
    a.spines['left'].set_linewidth(0.6)
    a.spines['bottom'].set_linewidth(0.6)
fig.savefig(P / 'output/Figure7B.pdf')
fig.savefig(P / 'output/Figure7B.svg')
fig.savefig(P / 'output/Figure7B.png', dpi=180)
plt.close(fig)
