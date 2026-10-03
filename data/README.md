`row_0.npz` through `row_4.npz` are the raw per-row outputs of the
exploratory (mu, beta) phase-map sweep in `../phase_map.py` (one row =
one value of mu, all tested values of beta). These are retained for
transparency, as described in the main README's "Status and known gaps"
section: this sweep did not produce a validated result (it uncovered
that `r = Kp*Vc/Vm`, not mu/beta, dominates practical identifiability,
but confirming this under realistic noise requires a sampling schedule
rescaled per r, which was not implemented). Combine with
`python phase_map.py combine` if you want to inspect the full grid.
