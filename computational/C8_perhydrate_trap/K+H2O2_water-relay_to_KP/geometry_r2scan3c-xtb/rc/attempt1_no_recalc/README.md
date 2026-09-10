# Attempt 1: killed at cycle 131, RFO Hessian drifting

Plain `Calc_Hess true` (no periodic refresh), `opt freq`. Converged steadily
through cycle ~111 (RMS gradient 1.85e-4 vs 1e-4 tolerance, MAX gradient
1.87e-3 vs 3e-4), then regressed by cycle 131 (RMS 2.25e-4, MAX 2.38e-3) --
the same RFO-Hessian-drift signature the TS searches showed, milder since a
minimization's BFGS/RFO update is more forgiving than TS eigenvector-
following. `job_cycle131.xyz` is where it was killed; attempt 2 restarts
from here with `Recalc_Hess 50` added rather than from scratch.
