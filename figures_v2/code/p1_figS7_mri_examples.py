# Figure S7: MRI examples of the seven LDD features, cropped at native resolution
# from the original Figure 2C (Figures/Figure 2.tif, 10000 x 13355 px).
import sys, os
import numpy as np
from PIL import Image
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from style import setup, FEAT as PAL, save, COL2
src = sys.argv[1]   # folder with C_left.png / C_right.png crops of panel C
out = sys.argv[2]
L = Image.open(os.path.join(src, "C_left.png")).convert("RGB")
R = Image.open(os.path.join(src, "C_right.png")).convert("RGB")
W, H = 330, 750
boxes = {  # feature: (image, x, y) of the inset in the crops
    "SS": (L, 494, 2257), "DB": (L, 495, 1415), "DN": (L, 494, 3098), "DH": (L, 1480, 205),
    "HIZ": (L, 494, 206), "MC": (R, 518, 1860), "SN": (R, 518, 3098)}
names = {"SS": "Schneiderman score\n(nucleus pulposus)", "DB": "Disc bulging\n(annulus fibrosus)",
         "DN": "Disc narrowing\n(nucleus pulposus)", "DH": "Disc herniation\n(annulus fibrosus)",
         "HIZ": "High-intensity zone\n(annulus fibrosus)", "MC": "Modic change\n(vertebral bone)",
         "SN": "Schmorl's node\n(endplate)"}
setup()
fig, axes = plt.subplots(1, 7, figsize=(COL2, 3.3))
for ax, (k, (im, x, y)) in zip(axes, boxes.items()):
    ax.imshow(np.asarray(im.crop((x, y, x + W, y + H))), interpolation="none")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(True); s.set_edgecolor(PAL[k]); s.set_linewidth(1.5)
    ax.set_title(k, color=PAL[k], fontsize=9, fontweight="bold", pad=3)
    ax.set_xlabel(names[k], fontsize=6, labelpad=3)
fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.14, wspace=0.08)
save(fig, out, "FigS7_MRI_feature_examples")
