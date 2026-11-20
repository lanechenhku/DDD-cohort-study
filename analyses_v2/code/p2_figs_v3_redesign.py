"""Paper 2 figure redesign (v3, 2026-11-06, after Lane: "Figures 2 and 3 are too weak").
Figure 2  rs1549596 and Modic changes in the cohort: age, level, sex, specificity.
Figure 3  External and biological context: large-scale look-ups, ITGA9 single-cell expression (HPA), hypothesis schematic.
Supp. Figure S7  ITGA9 eQTL look-up of rs1549596 and colocalisation (coloc.abf) with the MC.avg signal.
Inputs: results/p2_fig2_mc_prevalence_by_carrier.csv (p2_fig2_data.py), p2_rs1549596_sensitivity.csv, p2_rs1549596_by_level.csv,
p2_phenotype_descriptives.csv, p2_rs1549596_spine_lookups.csv, p2_go2_rs1549596.csv, p2_ITGA9_HPA_singlecell.tsv,
p2_rs1549596_ITGA9_eQTLCatalogue.csv, p2_rs1549596_ITGA9_GTEx_dyneqtl.csv, p2_itga9_coloc.csv."""
import os, sys, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
sys.path.insert(0, os.path.dirname(__file__)); import style as S
from common_v2 import RES
S.setup(); OUT = os.path.join(RES, "..", "..", "figures_v2", "Paper2")
R = lambda f, **k: pd.read_csv(os.path.join(RES, f), **k)
NC, CA = "#A9B3BF", S.FEAT["MC"]          # non-carrier, carrier
def fp(v): return f"{v:.1e}".replace("e-0", "×10⁻").replace("e-", "×10⁻") if v < 1e-3 else f"{v:.3f}" if v < 0.01 else f"{v:.2f}"
sup = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
def fp(v):
    if v >= 0.01: return f"{v:.2f}"
    if v >= 1e-3: return f"{v:.3f}"
    m, e = f"{v:.1e}".split("e"); return f"{m}×10$^{{{int(e)}}}$"

# ---------------- Figure 2 ----------------
pv = R("p2_fig2_mc_prevalence_by_carrier.csv"); se = R("p2_rs1549596_sensitivity.csv").set_index("analysis")
bl = R("p2_rs1549596_by_level.csv"); ds = R("p2_phenotype_descriptives.csv").set_index("phenotype")
fig = plt.figure(figsize=(S.COL2, 118 * S.MM))
gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.38, left=0.07, right=0.985, top=0.95, bottom=0.08)
# A: prevalence by age at last scan
ax = fig.add_subplot(gs[0, 0]); a = pv[pv.kind == "age"]; grp = list(dict.fromkeys(a.group)); xs = np.arange(len(grp))
for c, col, lab, off in [(0, NC, "Non-carriers (CC)", -0.08), (1, CA, "Carriers (CT/TT)", 0.08)]:
    s = a[a.carrier == c].set_index("group").loc[grp]
    ax.errorbar(xs + off, s.prev * 100, yerr=[(s.prev - s.lo) * 100, (s.hi - s.prev) * 100], color=col, marker="o", ms=4,
                lw=1.2, capsize=1.5, elinewidth=0.7, label=lab)
    for i, (k, n) in enumerate(zip(s.cases, s.n)):
        ax.text(xs[i] + off, -6.5 if c == 0 else -10.5, f"{k}/{n}", ha="center", fontsize=5, color=col if c else "#6B7480")
ax.set_xticks(xs); ax.set_xticklabels(grp); ax.set_xlabel("Age at last MRI (years)", labelpad=13)
ax.set_ylabel("Participants with any MC (%)"); ax.set_ylim(-12, 60); ax.set_yticks(range(0, 61, 10)); ax.spines["bottom"].set_position(("data", -12))
ax.axhline(0, color=S.LIGHT, lw=0.5, zorder=0); ax.legend(loc="upper left", fontsize=6)
o = se.loc["Any MC at any scan (logistic)"]; lo, hi = np.exp(o.beta - 1.96 * o.se), np.exp(o.beta + 1.96 * o.se)
ax.text(0.03, 0.74, f"Any MC, OR per T allele {o.OR:.2f}\n(95% CI {lo:.2f}–{hi:.2f}; p = {fp(o.p)})", transform=ax.transAxes,
        ha="left", va="top", fontsize=5.8)
ax.set_title("MC by age and carrier status", loc="left", pad=8); S.panel(ax, "A", x=-0.15, y=1.03)
# B: by disc level (spine order top->bottom)
ax = fig.add_subplot(gs[0, 1]); L = pv[pv.kind == "level"]; lev = ["L1/2", "L2/3", "L3/4", "L4/5", "L5/S1"]; ys = np.arange(5)[::-1]
for c, col, off in [(0, NC, 0.18), (1, CA, -0.18)]:
    s = L[L.carrier == c].set_index("group").loc[lev]
    ax.barh(ys + off, s.prev * 100, height=0.34, color=col, lw=0)
    ax.errorbar(s.prev * 100, ys + off, xerr=[(s.prev - s.lo) * 100, (s.hi - s.prev) * 100], fmt="none", ecolor="#444", elinewidth=0.5, capsize=1)
ax.set_yticks(ys); ax.set_yticklabels(lev); ax.set_xlabel("Participants with MC at the level, any scan (%)"); ax.set_xlim(0, 36)
from matplotlib.patches import Patch as _P
ax.legend([_P(color=NC), _P(color=CA)], ["Non-carriers", "Carriers"], loc="upper right", bbox_to_anchor=(1, 1.13), ncol=2, fontsize=5.8, handlelength=1)
for y, (_, r) in zip(ys, bl.iterrows()):
    ax.text(35.8, y, f"OR {r.OR:.2f}  p = {fp(r.p)}", ha="right", va="center", fontsize=5.8,
            color=S.INK if r.p < 0.05 else S.GREY, fontweight="bold" if r.p < 0.05 else "normal")
ax.set_title("MC by lumbar level", loc="left", pad=8); ax.tick_params(axis="y", length=0); S.panel(ax, "B", x=-0.15, y=1.03)
# C: by sex
ax = fig.add_subplot(gs[1, 0]); X = pv[pv.kind == "sex"]; sx = ["Women", "Men"]; xs = np.arange(2)
for c, col, off in [(0, NC, -0.19), (1, CA, 0.19)]:
    s = X[X.carrier == c].set_index("group").loc[sx]
    ax.bar(xs + off, s.prev * 100, width=0.36, color=col, lw=0)
    ax.errorbar(xs + off, s.prev * 100, yerr=[(s.prev - s.lo) * 100, (s.hi - s.prev) * 100], fmt="none", ecolor="#444", elinewidth=0.5, capsize=1.2)
    for i, (k, n) in enumerate(zip(s.cases, s.n)): ax.text(xs[i] + off, 1, f"{k}/{n}", ha="center", fontsize=5, color="white" if c else S.INK, rotation=90, va="bottom")
w, m, ix = se.loc["Women only"], se.loc["Men only"], se.loc["Interaction rs1549596 x sex"]
ax.set_xticks(xs); ax.set_xticklabels([f"Women\nMC.avg β = {w.beta:.2f}, p = {fp(w.p)}", f"Men\nMC.avg β = {m.beta:.2f}, p = {fp(m.p)}"], fontsize=6)
ax.set_ylabel("Participants with any MC (%)"); ax.set_ylim(0, 45)
ax.text(0.5, 0.97, f"rs1549596 × sex interaction p = {fp(ix.p)}", transform=ax.transAxes, ha="center", va="top", fontsize=6)
ax.set_title("MC by sex and carrier status", loc="left", pad=8); S.panel(ax, "C", x=-0.15, y=1.03)
# D: specificity across GWAS phenotypes (standardised beta per T allele)
ax = fig.add_subplot(gs[1, 1])
rows = [("MC.avg", "Primary model, hard calls (replicates GWAS)", "MC", S.FEAT["MC"]), ("Core.slope", "Core.slope (specificity)", "Core", S.FEAT["Core"]),
        ("DH.avg", "DH.avg (specificity)", "DH", S.FEAT["DH"]), ("HIZ.avg", "HIZ.avg (specificity)", "HIZ", S.FEAT["HIZ"])]
for i, (lab, key, ph, col) in enumerate(rows):
    r = se.loc[key]; sd = ds.loc[ph, "sd"]; b, e = r.beta / sd, r.se / sd; y = len(rows) - 1 - i
    ax.errorbar(b, y, xerr=1.96 * e, fmt="o", color=col, ms=4.5, elinewidth=1.1, capsize=0)
    ax.text(0.73, y, f"{b:.2f} SD, p = {fp(r.p)}", va="center", fontsize=5.8, color=S.INK)
ax.axvline(0, color=S.GREY, lw=0.6, ls="--"); ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows][::-1])
ax.set_xlim(-0.2, 1.12); ax.set_xticks([-0.2, 0, 0.2, 0.4, 0.6]); ax.set_xlabel("Effect per T allele (SD units, 95% CI)")
ax.tick_params(axis="y", length=0); ax.set_title("Effect across the four GWAS phenotypes", loc="left", pad=8); S.panel(ax, "D", x=-0.15, y=1.03)
S.save(fig, OUT, "Fig2_rs1549596_MC_in_cohort")

# ---------------- Figure 3 ----------------
fig = plt.figure(figsize=(S.COL2, 140 * S.MM))
gs = fig.add_gridspec(2, 2, height_ratios=[1.3, 0.75], width_ratios=[1.25, 1], hspace=0.22, wspace=0.62, left=0.27, right=0.985, top=0.97, bottom=0.03)
# A: forest of look-ups
ax = fig.add_subplot(gs[0, 0])
F = []
for key, lab in [("Any MC at any scan (logistic)", "Any MC, any scan"), ("Incident MC T1->T2 among MC-free at T1 (logistic)", "Incident MC")]:
    r = se.loc[key]; F.append(("This study (Hong Kong, EAS)", lab, r.beta, r.se, {"Any MC, any scan": "363 cases / 2,101", "Incident MC": "n = 1,317 MC-free at scan 1"}[lab], "DDD", r.p))
sl = R("p2_rs1549596_spine_lookups.csv")
for _, r in sl.iterrows():
    F.append(("Back pain and disc disorders", f"{r.phenotype} ({r.source.split(',')[0].replace(' et al.', '').split(' M13_')[0]}, {r.ancestry})", r.beta_logOR, r.se, r.cases_controls.split(" /")[0] + " cases", "ext", r.p))
go = R("p2_go2_rs1549596.csv")
for _, r in go[go.phenotype.isin(["All OA", "Spine OA", "Knee OA", "Hip OA", "Hand OA", "Total joint replacement"])].iterrows():
    F.append(("Osteoarthritis (GO 2.0, multi-ancestry)", r.phenotype, np.log(r.OR_T), (np.log(r.ci_hi) - np.log(r.ci_lo)) / 3.92, f"{r.ncases:,} cases", "ext", r.p))
F = pd.DataFrame(F, columns=["grp", "lab", "b", "se", "cases", "src", "p"])
y = 0; yt, yl = [], []
for gname, s in F.groupby("grp", sort=False):
    ax.text(-0.02, y, gname, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.2, fontweight="bold"); y -= 1
    for _, r in s.iterrows():
        col = S.FEAT["MC"] if r.src == "DDD" else "#55606E"
        ax.errorbar(np.exp(r.b), y, xerr=[[np.exp(r.b) - np.exp(r.b - 1.96 * r.se)], [np.exp(r.b + 1.96 * r.se) - np.exp(r.b)]], fmt="s" if r.src == "DDD" else "o",
                    color=col, ms=3.6 if r.src == "DDD" else 2.6, elinewidth=0.9, capsize=0)
        yt.append(y); yl.append(f"{r.lab}  [{r.cases}]"); y -= 1
    y -= 0.4
ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=5.6); ax.set_xscale("log"); ax.set_xlim(0.7, 3.6)
ax.set_xticks([0.8, 1, 1.25, 1.5, 2, 3]); ax.set_xticklabels(["0.8", "1", "1.25", "1.5", "2", "3"]); ax.minorticks_off()
ax.axvline(1, color=S.GREY, lw=0.6, ls="--"); ax.set_xlabel("Odds ratio per T allele (95% CI)"); ax.tick_params(axis="y", length=0)
ax.set_ylim(y + 0.6, 1.6); S.panel(ax, "A", x=-1.0, y=1.0)
# B: HPA single-cell
ax = fig.add_subplot(gs[0, 1])
h = R("p2_ITGA9_HPA_singlecell.tsv", sep="\t", header=None, names=["ens", "gene", "cell", "ncpm"]).sort_values("ncpm", ascending=False).head(20)[::-1]
def lin(c):
    if c in ("hematopoietic stem cells", "megakaryocyte-erythroid progenitors", "megakaryocyte progenitors", "neutrophil progenitors", "kupffer cells", "mast cells", "macrophages", "platelets", "microglia"): return "#7461AA"
    if "endothelial" in c: return "#C8443A"
    if c in ("fibro-adipogenic progenitors", "hepatic stellate cells", "peritubular myoid cells", "thymic myoid cells", "adipocytes"): return "#3A9A5B"
    return "#B9BEC6"
ax.barh(range(len(h)), h.ncpm, color=[lin(c) for c in h.cell], height=0.72, lw=0)
ax.set_yticks(range(len(h))); ax.set_yticklabels([c[0].upper() + c[1:] for c in h.cell], fontsize=5.6); ax.tick_params(axis="y", length=0)
ax.set_xlabel("ITGA9 expression (nCPM, HPA)"); ax.xaxis.grid(True); ax.set_axisbelow(True)
from matplotlib.patches import Patch
ax.legend([Patch(color=c) for c in ["#7461AA", "#C8443A", "#3A9A5B", "#B9BEC6"]], ["Haematopoietic / myeloid", "Endothelial", "Mesenchymal / stromal", "Other"],
          loc="lower right", fontsize=5.4, handlelength=1, borderaxespad=0.2)
S.panel(ax, "B", x=-0.75, y=1.0)
# C: hypothesis schematic
ax = fig.add_subplot(gs[1, :]); ax.set_position([0.03, 0.02, 0.955, 0.27]); ax.set_xlim(0, 100); ax.set_ylim(2, 34); ax.axis("off")
def box(x, y, w, hh, title, body, fc, ec, dashed=False):
    ax.add_patch(FancyBboxPatch((x, y), w, hh, boxstyle="round,pad=0.4,rounding_size=1.2", fc=fc, ec=ec, lw=0.9, ls="--" if dashed else "-"))
    ax.text(x + w / 2, y + hh - 1.6, title, ha="center", va="top", fontsize=6.6, fontweight="bold", color=S.INK)
    ax.text(x + w / 2, y + hh - 5.6, body, ha="center", va="top", fontsize=5.5, color=S.INK, linespacing=1.35)
def arr(x1, y1, x2, y2, dashed=False):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=8, lw=0.9, color="#444", ls="--" if dashed else "-"))
Y, Hh = 10, 20
box(0.8, Y, 17.5, Hh, "rs1549596-T", "Intron 1 of ITGA9\nCTCF motif; enhancer\nGH03J037451\nT frequency 6.4% (EAS),\n~30% (EUR)", "#EFEDF6", S.FEAT["MC"])
box(21.8, Y, 17.5, Hh, "Lower ITGA9?", "Weak eQTL in blood and\nmonocytes (T lowers\nexpression), opposite in\nadipose; colocalisation\nnot supported (Fig. S7)", "#F5F5F5", S.GREY, dashed=True)
box(42.8, Y, 17.5, Hh, "Integrin α9β1", "Receptor for tenascin-C,\nosteopontin, VCAM-1,\nfibronectin; supports\nCD34$^+$ HSPC adhesion\nand β-catenin signalling", "#F5F5F5", S.GREY)
box(63.8, Y, 17.5, Hh, "Vertebral marrow", "Altered haematopoietic,\nvascular and stromal niche\nbeneath the endplate;\ninflammatory (type 1) or\nfatty (type 2) remodelling", "#F5F5F5", S.GREY, dashed=True)
box(84.3, Y, 14.8, Hh, "Modic changes", "OR 2.0 per T allele\nmainly L4/5 and L5/S1\nstronger in women\nweak or no effect on\nCore, DH, HIZ, OA, pain", "#EFEDF6", S.FEAT["MC"])
for x in (18.8, 39.8, 60.8): arr(x + 0.2, Y + Hh / 2, x + 2.6, Y + Hh / 2, dashed=(x == 18.8))
arr(81.8, Y + Hh / 2, 83.9, Y + Hh / 2, dashed=True)
ax.text(50, 4.5, "Solid outline: observed in this study or public data.   Dashed outline / arrows: hypothesised links, not yet tested in disc, endplate or marrow tissue.",
        ha="center", fontsize=5.6, color=S.GREY, style="italic")
S.panel(ax, "C", x=0.0, y=1.0)
S.save(fig, OUT, "Fig3_ITGA9_context_mechanism")

# ---------------- Supplementary Figure S7 ----------------
q = R("p2_rs1549596_ITGA9_eQTLCatalogue.csv"); gt = R("p2_rs1549596_ITGA9_GTEx_dyneqtl.csv")
q = q[q.p < 0.05].copy(); q["lab"] = q.study.str.replace("_", " ") + " – " + q.group.str.replace("_", " ") + " (n=" + q.n.astype(str) + ")"
q = q.sort_values("beta")
co = R("p2_itga9_coloc.csv"); co = co[co.trait == "MC"]
fig, axs = plt.subplots(1, 2, figsize=(S.COL2, max(60, 4.2 * len(q) + 20) * S.MM), gridspec_kw=dict(width_ratios=[1.3, 1], wspace=1.05))
ax = axs[0]
for i, (_, r) in enumerate(q.iterrows()):
    ax.errorbar(r.beta, i, xerr=1.96 * r.se, fmt="o", ms=3, color=S.FEAT["MC"] if r.beta < 0 else S.FEAT["DH"], elinewidth=0.8)
ax.set_yticks(range(len(q))); ax.set_yticklabels(q.lab, fontsize=5.4); ax.axvline(0, color=S.GREY, ls="--", lw=0.6)
ax.set_xlabel("Effect of T allele on ITGA9 expression (β, 95% CI)"); ax.tick_params(axis="y", length=0)
ax.set_title(f"eQTL datasets with p < 0.05 ({len(q)} of {len(R('p2_rs1549596_ITGA9_eQTLCatalogue.csv'))})", loc="center", fontsize=6.5)
S.panel(ax, "A", x=-1.0, y=1.04)
ax = axs[1]; co = co.iloc[::-1]
H = ["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]; cols = ["#E6E6E6", "#C9CCD1", "#A9B3BF", "#C8443A", S.FEAT["MC"]]
left = np.zeros(len(co))
for hcol, c in zip(H, cols):
    ax.barh(range(len(co)), co[hcol], left=left, color=c, height=0.6, lw=0, label=hcol.replace("PP.", ""))
    left += co[hcol].values
ax.set_yticks(range(len(co))); ax.set_yticklabels(co.label, fontsize=5.6); ax.set_xlim(0, 1); ax.set_xlabel("Posterior probability (coloc.abf)")
ax.legend(ncol=5, loc="upper center", bbox_to_anchor=(0.4, -0.25), fontsize=5.4, handlelength=1, columnspacing=0.8)
ax.set_title("Colocalisation with MC.avg GWAS", loc="center", fontsize=6.5); ax.tick_params(axis="y", length=0)
S.panel(ax, "B", x=-0.75, y=1.04)
S.save(fig, OUT, "FigS7_ITGA9_eQTL_coloc")
print("done")
