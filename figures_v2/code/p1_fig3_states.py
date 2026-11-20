"""Paper 1 Figure 3: co-occurrence and progression of Core and complication features.

A  frequency of composite states at each scan (Wilson 95% CI)
B  transitions between consecutive scans (T1->T2 and T2->T3 pooled)
Both from all.emergence.data.20201215.csv; Core = SS or DB present at any level.
Participants with HIZ/MC but no Core (n = 2 at T1, 1 at T3) are omitted.
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import matplotlib.colors as mc
from matplotlib.patches import Patch
import style as S
import derive as Dv
from common import P1

S.setup()


def shade(c, f):
    r = np.array(mc.to_rgb(c)); return tuple(r + (1 - r) * f)


def main():
    e = Dv.emergence()
    st = Dv.states(e)
    fr = Dv.state_freq(st)
    tr = Dv.transitions(st)
    fig = plt.figure(figsize=(S.COL2, 150 * S.MM))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.15], hspace=0.42, left=0.11, right=0.985, top=0.95, bottom=0.07)
    ax = fig.add_subplot(gs[0])
    w = 0.26
    fades = [0.55, 0.25, 0.0]
    ns = fr.groupby("time").n.first()
    for i, s in enumerate(Dv.STATES):
        d = fr[fr.state == s].set_index("time").reindex(["T1", "T2", "T3"])
        xs = i + (np.arange(3) - 1) * w
        ax.bar(xs, d.p * 100, w * 0.9, color=[shade(S.STATE[s], f) for f in fades], edgecolor="none")
        ax.errorbar(xs, d.p * 100, yerr=[(d.p - d.lo) * 100, (d.hi - d.p) * 100], fmt="none",
                    ecolor="#3A3A3A", elinewidth=0.5, capsize=1.2, capthick=0.5)
        for x, t in zip(xs, ["T1", "T2", "T3"]):
            ax.text(x, -1.5, t, ha="center", va="top", fontsize=5.5, color=S.GREY)
    ax.set_xticks(range(len(Dv.STATES)))
    ax.set_xticklabels(Dv.STATES, fontweight="bold")
    ax.tick_params(axis="x", length=0, pad=11)
    for t, s in zip(ax.get_xticklabels(), Dv.STATES):
        t.set_color(S.STATE[s] if s != "Normal" else "#6B7078")
    ax.set_ylabel("Participants (%)")
    ax.set_ylim(0, 56)
    ax.yaxis.grid(True); ax.set_axisbelow(True)
    h = [Patch(color=shade("#555555", f), label=f"{t} (n = {ns[t]:,})") for t, f in zip(["T1", "T2", "T3"], fades)]
    ax.legend(handles=h, title="Scan", loc="upper right")
    S.panel(ax, "A", x=-0.085, y=1.02)

    bx = fig.add_subplot(gs[1])
    P = tr.div(tr.sum(1), axis=0) * 100
    order = Dv.STATES[::-1]
    for k, s in enumerate(order):
        left = 0
        for t in Dv.STATES:
            v = P.loc[s, t]
            if v <= 0:
                continue
            bx.barh(k, v, left=left, color=S.STATE[t], edgecolor="white", lw=0.4, height=0.72)
            if v >= 7:
                bx.text(left + v / 2, k, f"{v:.0f}", ha="center", va="center", fontsize=5.8,
                        color="white" if t not in ("Normal",) else S.INK)
            left += v
        bx.text(101.5, k, f"n = {int(tr.loc[s].sum()):,}", va="center", fontsize=6, color=S.GREY)
    bx.set_yticks(range(len(order)))
    bx.set_yticklabels(order, fontweight="bold")
    for t, s in zip(bx.get_yticklabels(), order):
        t.set_color(S.STATE[s] if s != "Normal" else "#6B7078")
    bx.tick_params(axis="y", length=0)
    bx.spines["left"].set_visible(False)
    bx.set_xlim(0, 100)
    bx.set_xlabel("State at next scan (%)")
    bx.set_ylabel("State at earlier scan")
    h = [Patch(color=S.STATE[s], label=s) for s in Dv.STATES]
    bx.legend(handles=h, ncol=7, loc="lower center", bbox_to_anchor=(0.5, 1.0), handlelength=1.0,
              columnspacing=1.0, title="State at next scan", title_fontsize=6.5)
    S.panel(bx, "B", x=-0.085, y=1.1)
    S.save(fig, P1, "Fig3_states_transitions")
    fr.to_csv(P1 + "/Fig3A_state_frequencies.csv", index=False)
    tr.to_csv(P1 + "/Fig3B_transition_counts.csv")


if __name__ == "__main__":
    main()
