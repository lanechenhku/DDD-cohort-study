"""v2 figures (2026-10-23) drawn from the re-analysis results; shared style with figures_v2."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import style as S
from common_v2 import *
S.setup()
OUT1 = "/mnt/project-files/DDD_papers/figures_v2/Paper1"; OUT2 = "/mnt/project-files/DDD_papers/figures_v2/Paper2"
feats = ["Core", "HIZ", "MC", "DH"]
TM = {"zage": "Age", "zage2": "Age²", "male": "Sex (male)", "z_core_onset": "Core onset age", "z_BMI_age40": "BMI", "z_Height_age40": "Height",
      "z_packyr": "Smoking", "z_injury": "Lumbar injury", "z_sport": "Sports", "zage:male": "Age × sex (male)", "zage:z_BMI_age40": "Age × BMI",
      "zage:z_Height_age40": "Age × height", "zage:z_packyr": "Age × smoking", "zage:z_injury": "Age × injury", "zage:z_sport": "Age × sports"}
PG = {"EUR.BMI.PRS": "BMI (EUR)", "EAS.BMI.PRS": "BMI (EAS)", "Multiethnic.BMI.PRS": "BMI (multi-ancestry)", "EUR.height.PRS": "Height (EUR)",
      "EUR.LDL_C.PRS": "LDL-C (EUR)", "EAS.LDL_C.PRS": "LDL-C (EAS)", "Mixed.LDL_C.PRS": "LDL-C (mixed)", "EUR.HDL_C.PRS": "HDL-C (EUR)",
      "EAS.HDL_C.PRS": "HDL-C (EAS)", "EUR.T2DM.PRS": "T2D (EUR)", "EAS.T2DM.PRS": "T2D (EAS)", "EUR.lsBMD.PRS": "Lumbar-spine BMD (EUR)",
      "CHI.MD.PRS": "MDD (Chinese)", "EUR.backpain.PRS": "Back pain (EUR)", "EUR.multisitechrpain.PRS": "Multisite chronic pain (EUR)", "EUR.SpineOA.PRS": "Spine OA (EUR)"}

def fig4():
    m = pd.read_csv(p(RES, "p1_wgee_model1.csv")); m = m[m.model == "M1"]; m["t"] = m.term.map(TM)
    terms = list(TM.values())
    m2 = pd.read_csv(p(RES, "p1_wgee_model2.csv")); m3 = pd.read_csv(p(RES, "p1_wgee_model3.csv")); m3 = m3[m3.model == "M3"]
    def lab(t):
        for k, v in PG.items():
            if k in t: return ("Age × " if t.startswith("zage:") else "") + v + " PGS"
    fig = plt.figure(figsize=(S.COL2, 160 * S.MM))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.55, 1], hspace=0.45, wspace=0.10, left=0.16, right=0.99, top=0.93, bottom=0.07)
    y = np.arange(len(terms))[::-1]
    for j, f in enumerate(feats):
        ax = fig.add_subplot(gs[0, j]); d = m[m.outcome == f].set_index("t")
        ax.axvline(0, color=S.GREY, lw=0.6, ls=(0, (3, 2)), zorder=0)
        for t, yy in zip(terms, y):
            if yy % 2 == 0: ax.axhspan(yy - .5, yy + .5, color="#F4F5F7", zorder=-1, lw=0)
            if t not in d.index: ax.text(0, yy, "n/a", ha="center", va="center", fontsize=5.5, color=S.GREY); continue
            r = d.loc[t]; sig = r.p < 0.05; c = S.FEAT[f]
            ax.plot([r.lo, r.hi], [yy, yy], color=c, lw=1.1, solid_capstyle="round")
            ax.plot(r.est, yy, "o", ms=3.6, mfc=c if sig else "white", mec=c, mew=0.9)
        ax.set_title(f, color=S.FEAT[f], fontsize=8); ax.set_ylim(-.6, len(terms) - .4); ax.set_yticks(y)
        ax.set_yticklabels(terms if j == 0 else []); ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
        ax.set_xlabel("β (95% CI), SD units"); ax.set_xlim(-0.6, 0.85); ax.set_xticks([-0.4, 0, 0.4, 0.8])
        if j == 0: S.panel(ax, "A", x=-0.75, y=1.02)
    def forest(ax, d, letter, title):
        d = d.reset_index(drop=True); yy = np.arange(len(d))[::-1]
        ax.axvline(0, color=S.GREY, lw=0.6, ls=(0, (3, 2)), zorder=0)
        for k, r in d.iterrows():
            c = S.FEAT[r.outcome]
            if yy[k] % 2 == 0: ax.axhspan(yy[k] - .5, yy[k] + .5, color="#F4F5F7", zorder=-1, lw=0)
            ax.plot([r.lo, r.hi], [yy[k]] * 2, color=c, lw=1.2); ax.plot(r.est, yy[k], "o", ms=4, mfc=c, mec=c)
            ax.text(1.01, yy[k], f"{r.est:.3f}  p={r.p:.3f}", transform=ax.get_yaxis_transform(), va="center", fontsize=5.6)
        ax.set_yticks(yy); ax.set_yticklabels([f"{lab(t)} → {o}" for t, o in zip(d.term, d.outcome)]); ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False); ax.set_ylim(-.6, len(d) - .4); ax.set_xlabel("β (95% CI), SD units"); ax.set_title(title, loc="left", fontweight="normal")
        S.panel(ax, letter, x=-1.25, y=1.04)
    axB = fig.add_subplot(gs[1, 0:2]); forest(axB, m2[m2.p < 0.05].sort_values("p"), "B", "Model 2 (one PGS at a time), p < 0.05")
    axC = fig.add_subplot(gs[1, 2:4]); forest(axC, m3[m3.term.str.contains("PRS") & (m3.p < 0.05)].sort_values("p"), "C", "Model 3 (PGS jointly), p < 0.05")
    axB.set_position([0.255, 0.07, 0.17, 0.25]); axC.set_position([0.745, 0.07, 0.15, 0.25])
    h = [Line2D([], [], marker="o", ls="-", color=S.FEAT[f], mfc=S.FEAT[f], ms=4, label=f) for f in feats]
    h += [Line2D([], [], marker="o", ls="", color=S.GREY, mfc="white", ms=4, label="p ≥ 0.05")]
    fig.legend(handles=h, loc="upper center", ncol=5, bbox_to_anchor=(0.57, 1.005))
    S.save(fig, OUT1, "Fig4_WGEE_risk_factors_v2")

def fig5():
    c = pd.read_csv(p(RES, "p1_pain_crosssectional.csv")); c = c[c.model == "D: + MRI + MCS"]
    outs = ["LBP in past year", "Sciatica ever", "ODI", "RMDQ", "VAS today", "VAS worst ever", "SF-36 PCS"]
    terms = [("zCore", "Core"), ("zHIZ", "HIZ"), ("zMC", "MC"), ("zDH", "DH"), ("zmcs1", "MCS")]
    col = {**S.FEAT, "MCS": "#3A9A5B"}
    fig, axes = plt.subplots(1, len(outs), figsize=(S.COL2, 52 * S.MM), sharey=True)
    for ax, o in zip(axes, outs):
        d = c[c.outcome == o].set_index("term"); logit = d.scale.iloc[0].startswith("OR")
        ax.axvline(1 if logit else 0, color=S.GREY, lw=.6, ls=(0, (3, 2)))
        for k, (t, nm) in enumerate(terms[::-1]):
            r = d.loc[t]; cc = col[nm]
            ax.plot([r.lo, r.hi], [k, k], color=cc, lw=1.1); ax.plot(r.est, k, "o", ms=3.6, mfc=cc if r.p < 0.05 else "white", mec=cc)
        ax.set_title(o, fontsize=6.8); ax.set_xlabel("OR per SD" if logit else "Δ per SD", fontsize=6)
        if logit: ax.set_xscale("log"); ax.set_xticks([0.6, 1, 1.6]); ax.set_xticklabels(["0.6", "1", "1.6"])
        ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    axes[0].set_yticks(range(len(terms))); axes[0].set_yticklabels([n for _, n in terms[::-1]])
    S.panel(axes[0], "", x=-0.5); fig.subplots_adjust(wspace=0.35, left=0.06, right=0.99, bottom=0.25, top=0.85)
    S.save(fig, OUT1, "Fig5_pain_mental_health")

def figS_onset():
    t = pd.read_csv(p(RES, "p1_turnbull_curves.csv")); s = pd.read_csv(p(RES, "p1_sise_curves.csv"))[["feature", "age", "F_logNe"]]   # v3: SISE re-run (p1_sise.R)
    fig, axes = plt.subplots(2, 4, figsize=(S.COL2, 85 * S.MM), sharex=True, sharey=True); axes = axes.ravel()
    fcol = [c for c in s.columns if c.lower() in ("feature", "feat")][0]; acol = [c for c in s.columns if "age" in c.lower()][0]
    ycol = [c for c in s.columns if c not in (fcol, acol)][0]
    for ax, f in zip(axes, ["SS", "DB", "DN", "DH", "HIZ", "MC", "SN"]):
        d = t[t.feature == f].sort_values("age"); ax.step(d.age, d.F, where="post", color=S.FEAT[f], lw=1.1, label="Turnbull NPMLE")
        e = s[s[fcol].astype(str).str.upper() == f].sort_values(acol)
        if len(e): ax.plot(e[acol], e[ycol], color=S.GREY, lw=1, ls="--", label="SISE")
        ax.axhline(.5, color=S.LIGHT, lw=.6); ax.set_title(f, color=S.FEAT[f]); ax.set_xlim(15, 85); ax.set_ylim(0, 1.02)
    axes[-1].axis("off"); axes[0].legend(loc="lower right", fontsize=5.5)
    for ax in axes[4:7]: ax.set_xlabel("Age (years)")
    axes[0].set_ylabel("Cumulative risk"); axes[4].set_ylabel("Cumulative risk")
    fig.subplots_adjust(hspace=.35, wspace=.12, left=.07, right=.99, top=.93, bottom=.1); S.save(fig, OUT1, "FigS5_onset_Turnbull_vs_SISE")

def figS_dbn():
    fig, axes = plt.subplots(1, 2, figsize=(S.COL2, 70 * S.MM))
    for ax, nm, order in [(axes[0], "feature", ["SS", "DB", "DN", "DH", "HIZ", "MC", "SN"]), (axes[1], "level", ["L1", "L2", "L3", "L4", "L5"])]:
        d = pd.read_csv(p(RES, f"p1_dbn_{nm}.csv")); M = np.zeros((len(order), len(order)))
        for _, r in d.iterrows(): M[order.index(r.frm), order.index(r.to)] = r.bootstrap_strength
        im = ax.imshow(M, cmap="Blues", vmin=0, vmax=1)
        for i in range(len(order)):
            for j in range(len(order)):
                if M[i, j] >= 0.05: ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=5.5, color="white" if M[i, j] > .6 else S.INK)
        ax.set_xticks(range(len(order))); ax.set_xticklabels(order); ax.set_yticks(range(len(order))); ax.set_yticklabels(order)
        ax.set_xlabel("Node at later scan (t+1)"); ax.set_ylabel("Node at earlier scan (t)"); ax.set_title(("Feature" if nm == "feature" else "Disc") + "-level DBN: bootstrap arc strength")
        for s in ax.spines.values(): s.set_visible(False)
    fig.colorbar(im, ax=axes, shrink=.7, label="Proportion of 200 bootstrap networks"); S.save(fig, OUT1, "FigS6_DBN_reestimated")

def figS_p2():
    sens = pd.read_csv(p(RES, "p2_rs1549596_sensitivity.csv")); lev = pd.read_csv(p(RES, "p2_rs1549596_by_level.csv"))
    g = pd.read_csv(p(UP, "Analyses_v2_2026-10-23/data/itga9_locus_genotypes.csv")).rename(columns={"3:37495010_T_C": "g"})[["IID", "g"]]
    ph = pd.read_csv(p(D8, "GWAS_phenotype.20210313.txt"), sep=" ").merge(g, on="IID").dropna(subset=["MC"])
    fig = plt.figure(figsize=(S.COL2, 120 * S.MM)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1], hspace=.45, wspace=.35)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])]
    ax = axes[0]; ph["gl"] = ph.g.map({0: "CC", 1: "CT", 2: "TT"})
    for k, gl in enumerate(["CC", "CT", "TT"]):
        v = ph[ph.gl == gl].MC.values; jit = np.random.default_rng(k).uniform(-.25, .25, len(v))
        ax.scatter(k + jit, v, s=2, color=S.FEAT["MC"], alpha=.35, lw=0); ax.plot([k - .3, k + .3], [v.mean()] * 2, color=S.INK, lw=1.2)
        ax.text(k, ph.MC.max() * 1.02, f"n={len(v)}\nmean={v.mean():.2f}", ha="center", fontsize=5.5)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["CC", "CT", "TT"]); ax.set_xlabel("rs1549596 genotype"); ax.set_ylabel("MC.avg"); S.panel(ax, "A", x=-.35)
    ax = axes[1]; y = np.arange(5)[::-1]
    ax.axvline(1, color=S.GREY, lw=.6, ls=(0, (3, 2)))
    for k, r in lev.iterrows():
        ax.plot([r.lo, min(r.hi, 8)], [y[k]] * 2, color=S.FEAT["MC"], lw=1.1); ax.plot(r.OR, y[k], "o", ms=3.6, mfc=S.FEAT["MC"] if r.p < .05 else "white", mec=S.FEAT["MC"])
        ax.text(1.02, y[k], f"{int(r.cases)} cases", va="center", fontsize=5.5, transform=ax.get_yaxis_transform())
    ax.set_xscale("log"); ax.set_xticks([0.25, 1, 4]); ax.set_xticklabels(["0.25", "1", "4"]); ax.set_yticks(y); ax.set_yticklabels([f"{l}" for l in lev.level])
    ax.set_xlabel("OR per T allele (MC at level, any scan)"); ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0); S.panel(ax, "B", x=-.25)
    ax = axes[2]; keep = sens[sens.model == "linear"].copy(); keep["z"] = keep.beta / keep.se
    keep = keep[~keep.analysis.str.contains("specificity|Interaction")]
    yy = np.arange(len(keep))[::-1]
    ax.barh(yy, -np.log10(keep.p), color=[S.FEAT["MC"]] * len(keep), height=.6)
    ax.axvline(-np.log10(5e-8), color="#C8443A", lw=.6, ls="--"); ax.axvline(-np.log10(0.05), color=S.GREY, lw=.6, ls=":")
    ax.set_yticks(yy); ax.set_yticklabels(keep.analysis, fontsize=5.6); ax.set_xlabel("−log10 p"); ax.tick_params(axis="y", length=0); S.panel(ax, "C", x=-0.42)
    ax.set_position([0.36, 0.07, 0.6, 0.38]); fig.subplots_adjust(left=.08, right=.92, top=.93); S.save(fig, OUT2, "FigS3_rs1549596_robustness")

if __name__ == "__main__":
    import sys
    for f in (sys.argv[1:] or ["fig4", "fig5", "figS_dbn", "figS_p2"]): globals()[f]()
