"""Figures 6-8 for Paper 1 (whole-spine). 2026-11-06."""
import sys, pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '/mnt/project-files/DDD_papers/analyses_v2/code'); import style as S
S.setup()
R = "/mnt/project-files/DDD_papers/analyses_ws/results"; OUT = "/mnt/project-files/DDD_papers/figures_v2/Paper1"
lab = ['C2-3','C3-4','C4-5','C5-6','C6-7','C7-T1'] + [f'T{i}-{i+1}' for i in range(1, 12)] + ['T12-L1'] + [f'L{i}-{i+1}' for i in range(1, 5)] + ['L5-S1']
NAME = {'PS': 'Pfirrmann ≥ IV', 'PDB': 'Posterior bulge', 'ADB': 'Anterior bulge', 'HIZ': 'HIZ', 'OMC': 'Modic (overall)', 'AMC': 'Modic (anterior)', 'SN': "Schmorl's node"}
COL = {'PS': S.FEAT['SS'], 'PDB': S.FEAT['DB'], 'ADB': '#E08A7E', 'HIZ': S.FEAT['HIZ'], 'OMC': S.FEAT['MC'], 'AMC': '#A99BD0', 'SN': S.FEAT['SN']}
FE = ['PS', 'PDB', 'ADB', 'HIZ', 'OMC', 'AMC', 'SN']
def regbars(ax, y0=-0.9):
    for x0, x1, t in [(-0.5, 5.5, 'Cervical'), (5.5, 17.5, 'Thoracic'), (17.5, 22.5, 'Lumbar')]:
        ax.axvline(x1, color='white', lw=1.2)
# ---------- Figure 6: atlas
p = pd.read_csv(f"{R}/ws_prevalence_by_level.csv"); pa = pd.read_csv(f"{R}/ws_prevalence_by_age_level.csv")
L = pd.read_csv(f"{R}/ws_long.csv")
L['ps4'] = np.where(L.feature == 'PS', (L.y >= 4).astype(float), np.nan); L.loc[L.y.isna(), 'ps4'] = np.nan
bins = [0, 35, 45, 55, 65, 100]; ag = ['<35', '35–44', '45–54', '55–64', '≥65']
L['agegrp'] = pd.cut(L.age, bins, labels=ag)
fig = plt.figure(figsize=(S.COL2, 150 * S.MM)); gs = fig.add_gridspec(3, 1, height_ratios=[1.25, 1, 1], hspace=0.55)
ax = fig.add_subplot(gs[0])
M = np.array([[p[(p.feature == f) & (p.level == l)]['prev_adv' if f == 'PS' else 'prev_any'].values[0] for l in lab] for f in FE])
im = ax.imshow(M * 100, aspect='auto', cmap='Blues', vmin=0, vmax=60)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j] * 100; ax.text(j, i, f'{v:.0f}' if v >= 0.5 else '·', ha='center', va='center', fontsize=5, color='white' if v > 35 else S.INK)
ax.set_yticks(range(len(FE)), [NAME[f] for f in FE]); ax.set_xticks(range(23), lab, rotation=90); regbars(ax)
cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01); cb.set_label('Prevalence (%)')
ax.set_title('Prevalence by disc level (n = 2,007)', loc='left'); S.panel(ax, 'A', x=-0.1)
for k, (f, col, title, letter) in enumerate([('PS', 'ps4', 'Pfirrmann ≥ IV by age group (%)', 'B'), ('OMC', None, 'Modic change (overall) by age group (%)', 'C')]):
    ax = fig.add_subplot(gs[k + 1])
    if f == 'PS': g = L[L.feature == 'PS'].groupby(['agegrp', 'level'], observed=True)['ps4'].mean().unstack()[lab]
    else:
        q = L[L.feature == f].copy(); q['b'] = (q.y > 0).astype(float).where(q.y.notna()); g = q.groupby(['agegrp', 'level'], observed=True)['b'].mean().unstack()[lab]
    vmax = 80 if f == 'PS' else 30
    im = ax.imshow(g.values * 100, aspect='auto', cmap='Blues' if f == 'PS' else 'Purples', vmin=0, vmax=vmax)
    ax.set_yticks(range(len(ag)), [f'{a} ({n})' for a, n in zip(ag, L[(L.feature == 'PS') & (L.level == 'C2-3')].agegrp.value_counts().reindex(ag).values)])
    ax.set_xticks(range(23), lab, rotation=90); regbars(ax); fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    ax.set_title(title, loc='left'); S.panel(ax, letter, x=-0.1)
S.save(fig, OUT, 'Fig6_whole_spine_atlas')
# ---------- Figure 7: spatial structure and lumbar course
c = pd.read_csv(f"{R}/ws_corr_by_lag.csv"); rc = pd.read_csv(f"{R}/ws_factor_regioncor.csv"); lc = pd.read_csv(f"{R}/ws_lumbar_course_vs_regions.csv")
fig, axs = plt.subplots(1, 3, figsize=(S.COL2, 66 * S.MM), gridspec_kw=dict(width_ratios=[1, 0.95, 1], wspace=0.95))
ax = axs[0]
for f in FE:
    s = c[c.feature == f]; ax.plot(s.lag, s.r_partial, '-', color=COL[f], lw=1.1, label=NAME[f].replace(' ≥ IV', ''))
ax.axhline(0, color=S.GREY, lw=0.5); ax.set_xlabel('Distance between discs (levels)'); ax.set_ylabel('Mean partial Spearman ρ\n(age, sex)')
ax.set_xlim(1, 22); ax.legend(fontsize=5.5, loc='upper right', handlelength=1.2); ax.set_title('Correlation by distance', loc='left'); S.panel(ax, 'A')
ax = axs[1]; pairs = ['C-T', 'L-T', 'C-L']; mk = {'C-T': 'o', 'L-T': 's', 'C-L': '^'}; pn = {'C-T': 'cervical–thoracic', 'L-T': 'thoracic–lumbar', 'C-L': 'cervical–lumbar'}
for i, f in enumerate(FE):
    for k, pr in enumerate(pairs):
        r = rc[(rc.feature == f) & (rc.pair == pr)]
        if len(r): ax.errorbar(r.r.values[0], i + (k - 1) * 0.22, xerr=1.96 * r.se.values[0], fmt=mk[pr], ms=3, color=COL[f], mfc='white' if pr == 'C-L' else COL[f], lw=0.7, capsize=0)
ax.set_yticks(range(len(FE)), [NAME[f].replace(' ≥ IV', '') for f in FE]); ax.invert_yaxis(); ax.set_xlim(-0.2, 1); ax.axvline(0, color=S.GREY, lw=0.5)
from matplotlib.lines import Line2D
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.36), ncol=3, columnspacing=0.6, handletextpad=0.2, handles=[Line2D([], [], marker=mk[p_], ls='', color=S.INK, mfc='white' if p_ == 'C-L' else S.INK, ms=3) for p_ in pairs], labels=['C–T', 'T–L', 'C–L'], fontsize=5.5)
ax.set_xlabel('Correlation between region factors\n(age- and sex-adjusted)'); ax.set_title('Shared between regions', loc='left'); S.panel(ax, 'B', x=-0.62, y=1.04)
ax = axs[2]; outs = [('PS_C', 'Pfirrmann, C'), ('PS_T', 'Pfirrmann, T'), ('PDB_C', 'Post. bulge, C'), ('PDB_T', 'Post. bulge, T'),
                     ('HIZ_C', 'HIZ, C'), ('OMC_C', 'Modic, C'), ('SN_T', "Schmorl's, T")]
pc = {'core_onset': ('Later lumbar Core onset', S.GREY, 'o'), 'core_rate': ('Faster lumbar Core progression', S.FEAT['Core'], 's')}
for i, (o, on) in enumerate(outs):
    for k, (pr, (nm, col, m)) in enumerate(pc.items()):
        r = lc[(lc.predictor == pr) & (lc.outcome == o)].iloc[0]
        ax.errorbar(r.beta, i + (k - 0.5) * 0.3, xerr=1.96 * r.se, fmt=m, ms=3, color=col, lw=0.7, capsize=0, label=nm if i == 0 else None)
ax.set_yticks(range(len(outs)), [o[1] for o in outs]); ax.invert_yaxis(); ax.axvline(0, color=S.GREY, lw=0.5)
ax.set_xlabel('SD change in outcome per SD of predictor'); ax.legend(fontsize=5.5, loc='lower center', bbox_to_anchor=(0.5, -0.33)); ax.set_title('Lumbar course vs C/T spine', loc='left'); S.panel(ax, 'C', x=-0.5)
S.save(fig, OUT, 'Fig7_spatial_structure')
# ---------- Figure 8: risk factors by region
rf = pd.read_csv(f"{R}/ws_riskfactors_by_region.csv"); tt = pd.read_csv(f"{R}/ws_riskfactors_region_tests.csv")
TERMS = [('age_z', 'Age (per SD)'), ('male', 'Male sex'), ('BMI_age40_z', 'BMI (per SD)'), ('Height_age40_z', 'Height (per SD)'), ('injury', 'Lumbar injury'), ('smoke', 'Ever smoked'), ('workload_z', 'Workload (per SD)')]
RC = {'C': '#8DA0CB', 'T': '#66C2A5', 'L': S.FEAT['Core']}
fig, axs = plt.subplots(1, 3, figsize=(S.COL2, 68 * S.MM), sharey=True, gridspec_kw=dict(wspace=0.38))
for ax, f, letter in zip(axs, ['PS', 'PDB', 'OMC'], 'ABC'):
    for i, (t, tn) in enumerate(TERMS):
        for k, r_ in enumerate(['C', 'T', 'L']):
            q = rf[(rf.feature == f) & (rf.term == t) & (rf.region == r_)]
            if len(q): ax.errorbar(q.beta.values[0], i + (k - 1) * 0.25, xerr=1.96 * q.se.values[0], fmt='o', ms=2.8, color=RC[r_], lw=0.7, capsize=0,
                                    label={'C': 'Cervical', 'T': 'Thoracic', 'L': 'Lumbar'}[r_] if i == 0 else None)
        pt = tt[(tt.feature == f) & (tt.term == t)]
        if len(pt) and pt.p.values[0] < 0.01: ax.text(1.02, i, 'p$_{het}$=' + (f'{pt.p.values[0]:.0e}' if pt.p.values[0] < 0.001 else f'{pt.p.values[0]:.3f}'), transform=ax.get_yaxis_transform(), ha='left', va='center', fontsize=5, color=S.INK, clip_on=False)
    ax.axvline(0, color=S.GREY, lw=0.5); ax.set_title({'PS': 'Pfirrmann grade', 'PDB': 'Posterior bulge', 'OMC': 'Modic change'}[f], loc='left'); S.panel(ax, letter, x=-0.1 if letter != 'A' else -0.45, y=1.08)
    ax.set_xlabel('SD difference in regional score')
axs[0].set_yticks(range(len(TERMS)), [t[1] for t in TERMS]); axs[0].invert_yaxis(); axs[0].legend(fontsize=5.5, loc='lower left')
S.save(fig, OUT, 'Fig8_regional_risk_factors')
print('ok')
