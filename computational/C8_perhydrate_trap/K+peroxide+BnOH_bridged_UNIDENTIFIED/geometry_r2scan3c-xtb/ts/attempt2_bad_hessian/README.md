# Attempt 2: killed, unreliable initial Hessian

Seeded from attempt 1's cycle-127 geometry (`job_seed_cycle127.xyz`) instead
of the scan guess, on the theory that starting closer to the saddle would
help a fresh `Calc_Hess true`. It didn't: **60% of the 439 numerical-Hessian
displacements failed** ("ORCA finished with an error in the QM2 calculation"
-- the GFN-xTB SCF, not the r2SCAN-3c QM region) before being killed at
displacement ~140/439.

Almost all the per-displacement scratch was already gone by the time this was
investigated (killed mid-run, most temp files cleaned up as ORCA proceeded),
so the exact geometric cause per failed displacement isn't recorded. The
reading: attempt 1's 127 cycles left the macrocycle somewhere a numerical
displacement routinely breaks GFN-xTB's SCF -- i.e. attempt 1 wasn't just
*slow*, it plausibly drifted somewhere structurally marginal. That reframes
attempt 1's plateau as a symptom worth taking more seriously than "just
needs more cycles."

**What attempt 3 does differently**: restarts from the scan's own guess
(`../scan/job.009.xyz`), which *did* complete `Calc_Hess true` cleanly in
attempt 1 with no displacement failures -- known-good starting geometry.
Adds `Recalc_Hess 10` (ORCA 6.1 manual, TS Searches: recommended when "the
PES near the TS can be very far from ideal for a Newton-Raphson step") so
the Hessian gets refreshed periodically through the whole search rather than
relying on one initial guess plus 100+ cycles of RFO updates -- both to
prevent the earlier stall and to surface a repeat of this displacement-
failure pattern early (at cycle 10, 20, ...) if the search drifts somewhere
bad again, instead of 127 cycles in.
