"""
make_all_figures.py
Reproduces every figure in the manuscript in one command:

    python make_all_figures.py

Figure 1 (schematic summary)         -> figures/Fig1.pdf, Fig1.png
Figure 2 (continuous equivalence     -> figures/Fig2.pdf, Fig2.png
          family, numerical check)
Figure 3 (FIM instability)           -> figures/Fig3.pdf, Fig3.png
Figure 4 (single- vs multi-start     -> figures/Fig4.pdf, Fig4.png
          profile likelihood)

Each script also prints, to stdout, the exact numerical values quoted in
the manuscript text (e.g. condition numbers, relative standard errors,
the 3.92e-07 maximum relative difference in the equivalence-family demo),
so the printed output of this script is itself a reproducibility check
against the numbers in Results.

Runtime: a few minutes total (Fig. 4's multi-start profile likelihood is
the slow part).
"""
import subprocess
import sys

SCRIPTS = [
    ("Figure 1 (schematic summary)", "figure1_schematic.py"),
    ("Figures 2-4 (degeneracy demo, FIM instability, profile likelihood)", "figures2to4.py"),
]

if __name__ == "__main__":
    for label, script in SCRIPTS:
        print(f"\n{'='*70}\n{label}\n{'='*70}")
        result = subprocess.run([sys.executable, script])
        if result.returncode != 0:
            print(f"FAILED: {script} exited with code {result.returncode}")
            sys.exit(result.returncode)
    print(f"\n{'='*70}\nAll manuscript figures reproduced in figures/\n{'='*70}")
    print("\nNote: the supplementary mu/beta phase-map figure (Discussion,")
    print("response to practical-identifiability-regime reviewer comment) is")
    print("generated separately by phase_map.py, since it is computationally")
    print("heavier (multi-start profile likelihood over a parameter grid) and")
    print("is not yet referenced by a figure number in the current manuscript draft.")
