"""
Figure 1 for the manuscript: three-panel schematic.
  A. Information bottleneck: Cm(t) -> {a11, a22, P} vs {CLd, Kp, CLint, kb}
  B. Continuous equivalence family parametrized by kb, collapsing to a point
     once kb is independently known (NSB control)
  C. Structural identifiability -> (necessary, not sufficient) -> practical
     estimability

Journal spec: multipart figure 178 mm wide, line art at 1000 dpi (text-bearing
figure -> use 400 dpi minimum, we target 600 dpi for safety), PDF output.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import matplotlib.lines as mlines

WIDTH_MM = 178
WIDTH_IN = WIDTH_MM / 25.4
HEIGHT_IN = WIDTH_IN * 0.48

plt.rcParams.update({
    "font.size": 8,
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "dejavusans",
})

fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(WIDTH_IN, HEIGHT_IN))
fig.subplots_adjust(wspace=0.5, top=0.84, bottom=0.04, left=0.03, right=0.99)


def box(ax, xy, w, h, text, fc="white", ec="black", fontsize=9.5, lw=1.2, style="round,pad=0.02"):
    b = FancyBboxPatch(xy, w, h, boxstyle=style, linewidth=lw, edgecolor=ec, facecolor=fc, zorder=2)
    ax.add_patch(b)
    ax.text(xy[0]+w/2, xy[1]+h/2, text, ha="center", va="center", fontsize=fontsize, zorder=3)
    return b


def arrow(ax, p0, p1, color="black", lw=1.4, style="-|>", connectionstyle="arc3,rad=0.0"):
    a = FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=11, linewidth=lw,
                         color=color, connectionstyle=connectionstyle, zorder=1)
    ax.add_patch(a)


# ---------------- Panel A: information bottleneck ----------------
axA.set_xlim(0, 10); axA.set_ylim(0, 10); axA.axis("off")
axA.set_title("A", fontsize=11, fontweight="bold", loc="left", pad=6)

box(axA, (1.0, 7.6), 4.0, 1.6, r"$C_m(t)$", fc="#eaf2fb")
arrow(axA, (3.0, 7.6), (3.0, 6.35))
box(axA, (0.3, 4.4), 5.4, 1.9, r"$\{a_{11},\,a_{22},\,P\}$", fc="#dff0d8")
arrow(axA, (3.0, 4.4), (3.0, 3.15), color="#b00000")
box(axA, (0.0, 0.9), 6.0, 2.1, r"$\{CL_d,\,K_p,\,CL_{int},\,k_b\}$", fc="#fbeaea")

axA.text(6.4, 2.7, "3 vs. 4", fontsize=9, color="#b00000", rotation=90, va="center")

# ---------------- Panel B: continuous equivalence family ----------------
axB.set_title("B", fontsize=11, fontweight="bold", loc="left", pad=6)
kb_range = np.linspace(0.2, 3.0, 200)
CLd_curve = 9.0 - 1.6*kb_range
CLint_curve = 4.0 / (0.5 + kb_range) + 1.0
axB.plot(CLd_curve, CLint_curve, color="#b00000", lw=2.2, zorder=2)
sample_idx = [10, 60, 120, 180]
for i in sample_idx:
    axB.plot(CLd_curve[i], CLint_curve[i], "o", color="#b00000", ms=6, zorder=3)

true_i = 90
axB.plot(CLd_curve[true_i], CLint_curve[true_i], "*", color="#1a6b1a", ms=20, zorder=4,
         markeredgecolor="black", markeredgewidth=0.6)
axB.annotate(r"$k_b$ known", xy=(CLd_curve[true_i], CLint_curve[true_i]),
             xytext=(CLd_curve[true_i]+1.0, CLint_curve[true_i]-1.6),
             fontsize=9, color="#1a6b1a",
             arrowprops=dict(arrowstyle="->", color="#1a6b1a", lw=1.0))

axB.set_xlabel(r"$CL_d$", fontsize=10)
axB.set_ylabel(r"$CL_{int}$", fontsize=10, labelpad=2)
axB.set_xticks([]); axB.set_yticks([])
axB.set_ylim(0, 9.5)
for spine in ["top", "right"]:
    axB.spines[spine].set_visible(False)

# ---------------- Panel C: structural vs practical ----------------
axC.set_xlim(0, 10); axC.set_ylim(0, 10); axC.axis("off")
axC.set_title("C", fontsize=11, fontweight="bold", loc="left", pad=6)

box(axC, (0.7, 7.9), 8.6, 1.4, "Independent NSB", fc="#eaf2fb")
arrow(axC, (5.0, 7.9), (5.0, 6.65))
box(axC, (0.7, 4.95), 8.6, 1.55, "Structural identifiability", fc="#dff0d8")
axC.text(5.0, 4.35, "necessary, not sufficient", fontsize=8.5, ha="center",
          style="italic", color="#b00000")
arrow(axC, (5.0, 4.95), (5.0, 3.15), color="#b00000")
box(axC, (0.7, 1.55), 8.6, 1.55, "Practical estimability", fc="#fdf2d0")

fig.savefig("/mnt/user-data/outputs/figure1_schematic.pdf", dpi=600, bbox_inches="tight")
fig.savefig("/mnt/user-data/outputs/figure1_schematic.png", dpi=600, bbox_inches="tight")
print("saved figure1_schematic.pdf and .png")
print(f"figure width target: {WIDTH_MM} mm ({WIDTH_IN:.3f} in)")
