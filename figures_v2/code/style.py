"""Shared figure style for the two LDD papers (HK Spine Cohort).

One palette, one font, journal sizes: single column 89 mm, 1.5 column 120 mm,
double column 183 mm. Text 7 pt, panel letters 9 pt bold. Exports vector PDF,
600-dpi TIFF (LZW) and a PNG preview.
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt

MM = 1 / 25.4
COL1, COL15, COL2 = 89 * MM, 120 * MM, 183 * MM

# MRI features: Core family in blues/slate, complications in distinct hues.
FEAT = {
    "SS": "#2C6FB7",
    "DB": "#C8443A",
    "DN": "#6F8FA6",
    "DH": "#E3912B",
    "HIZ": "#D4508A",
    "MC": "#7461AA",
    "SN": "#3A9A5B",
    "Core": "#1F4E8C",
    "Any": "#444444",
}
STATE = {
    "Normal": "#C9CCD1",
    "Core only": FEAT["Core"],
    "Core + HIZ": FEAT["HIZ"],
    "Core + MC": FEAT["MC"],
    "Core + DH": FEAT["DH"],
    "Core + two": "#2A9D8F",
    "All four": "#7A1F3D",
}
GREY = "#7F7F7F"
LIGHT = "#E6E6E6"
INK = "#222222"


def setup():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 7.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 7,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5,
        "legend.title_fontsize": 7,
        "legend.frameon": False,
        "axes.linewidth": 0.6,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "grid.color": LIGHT,
        "grid.linewidth": 0.5,
        "lines.linewidth": 1.0,
        "patch.linewidth": 0.4,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.dpi": 600,
        "figure.dpi": 150,
    })


def panel(ax, letter, x=-0.12, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left", color=INK)


def save(fig, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, name)
    fig.savefig(base + ".pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(base + ".tiff", dpi=600, bbox_inches="tight", pad_inches=0.02,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(base + ".png", dpi=200, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
