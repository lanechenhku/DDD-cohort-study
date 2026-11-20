"""Paper 1 Figure 1: design of the Hong Kong Spine Cohort (timeline + framework).

Bar spans follow the 2021 figure (Figures/Figure 1.tif); counts from the text.
"""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import style as S
from common import P1

S.setup()
MRI, QUE, GEN = "#5B8DC9", "#D9786F", "#E8A84A"
rows = [  # label, n, start, end, colour, row -- spans from the scan / questionnaire dates
    # (allpheno_v3.xlsx date1-3; DateOf1st/2ndQuestionnaire; 1st-99th percentile for questionnaires), checked 2026-11-06
    ("1st MRI", "3,463", 1990.7, 2017.1, MRI, 0),
    ("2nd MRI", "2,050", 2002.1, 2017.3, MRI, 1),
    ("3rd MRI", "450", 2003.2, 2018.0, MRI, 2),
    ("Whole-spine MRI (at 2nd MRI)", "2,007", 2004.5, 2017.3, "#2E5E9E", 3),
    ("1st questionnaire", "3,095", 2002.1, 2015.6, QUE, 4.3),
    ("2nd quest.", "950", 2009.1, 2017.8, QUE, 5.3),
]
fig = plt.figure(figsize=(S.COL2, 70 * S.MM))
ax = fig.add_axes([0.10, 0.17, 0.56, 0.74])
for lab, n, a, b, c, r in rows:
    ax.add_patch(FancyBboxPatch((a, r - 0.36), b - a, 0.72, boxstyle="round,pad=0,rounding_size=0.25",
                                fc=c, ec="none", alpha=0.92, mutation_aspect=0.08))
    ax.text((a + b) / 2, r, f"{lab}  {n}", ha="center", va="center", color="white", fontsize=6.4, fontweight="bold")
for y, lab in [(1.5, "Imaging"), (4.8, "Clinical, anthropometric\n& lifestyle"), (6.6, "Genetics")]:
    ax.text(1986.6, y, lab, ha="right", va="center", fontsize=6.8, color=S.INK, style="italic")
ax.set_xlim(1988, 2021)
ax.set_ylim(7.3, -0.7)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xticks(range(1990, 2021, 5))
ax.xaxis.grid(True, color=S.LIGHT, lw=0.5)
ax.set_axisbelow(True)
ax.set_xlabel("Calendar year")
ax.text(1990.5, 6.6, "Genome-wide genotyping 2,121  ·  whole-exome sequencing 685  (participants after QC)",
        ha="left", va="center", fontsize=6.2, color=S.INK)
S.panel(ax, "A", x=-0.17, y=1.0)

# B: analysis framework
bx = fig.add_axes([0.68, 0.06, 0.33, 0.90]); bx.axis("off"); bx.set_xlim(0, 10); bx.set_ylim(0, 10)
def box(x, y, w, h, t, c, bold=False):
    bx.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.35",
                                fc="white", ec=c, lw=1.2))
    bx.text(x, y, t, ha="center", va="center", fontsize=7, color=S.INK, fontweight="bold" if bold else "normal")
def arrow(p, q, c):
    bx.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=7, color=c, lw=1.0, shrinkA=0, shrinkB=0))
box(5.0, 8.6, 4.6, 1.4, "Onset & progression\n(Figures 2–3)", "#4C9A6A")
box(5.0, 5.0, 2.8, 1.4, "Disc\ndegeneration", S.INK, bold=True)
box(5.0, 1.4, 5.2, 1.4, "Risk factors (Figures 4, 8)\nnon-genetic · PGS", "#C8443A")
box(9.0, 5.0, 2.0, 1.4, "Pain &\nMCS", "#8C86C4")
box(1.2, 5.0, 2.0, 1.4, "Whole\nspine", "#2E5E9E")
arrow((5.0, 7.9), (5.0, 5.7), "#4C9A6A")
arrow((5.0, 2.1), (5.0, 4.3), "#C8443A")
arrow((6.4, 5.0), (7.97, 5.0), "#8C86C4")
arrow((2.23, 5.0), (3.6, 5.0), "#2E5E9E")
bx.text(9.0, 4.0, "(Figure 5)", ha="center", va="top", fontsize=6, color=S.GREY)
bx.text(1.2, 4.0, "(Figures 6–7)", ha="center", va="top", fontsize=6, color=S.GREY)
S.panel(bx, "B", x=0.0, y=0.93)
S.save(fig, P1, "Fig1_cohort_design")
