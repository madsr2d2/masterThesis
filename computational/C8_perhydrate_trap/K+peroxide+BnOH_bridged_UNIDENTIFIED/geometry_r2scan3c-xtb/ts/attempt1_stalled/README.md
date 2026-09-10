# Attempt 1: stalled, killed at cycle 127

`job.inp` here is exactly what ran (QM/XTB `r2SCAN-3c`/ALPB(water) `OptTS Freq`,
QM region `{146 130 128 131 137 139 136 133 132}`, seeded from `../scan`'s
energy maximum). `job_stalled_cycle127.xyz` is the geometry it had reached
when killed.

**What happened**: the Hessian's negative-eigenvalue count improved from 5
(early cycles) to a stable 2, then settled at 1 from cycle ~119 onward --
outwardly consistent with approaching a genuine first-order saddle. But the
actual convergence criteria never moved: RMS gradient sat at ~4-5e-4 (tol
1e-4) and MAX gradient at ~2.5-3.0e-3 (tol 3e-4) essentially unchanged from
cycle 84 through cycle 127, 40+ cycles with no real progress. `MaxIter` was
441, so it would not have self-stopped soon.

Read live with a `Monitor` watch on `job.out` (`GEOMETRY OPTIMIZATION CYCLE`,
`Hessian has ... negative eigenvalue`, `TS mode is mode number`, the
`Geometry convergence` table) -- the eigenvalue-count trend alone would have
looked like progress; only the actual gradient/energy numbers showed the
plateau.

**Likely cause**: the 130-atom cyclodextrin scaffold outside the QM region
can wobble each cycle in ways that perturb the QM region's local Hessian
update without moving the reactive coordinate -- a known pathology for
QM/XTB on a large, floppy host. The unresolved double-link-atom warning on
atom 128 (bonded to both 120 and 124) is a plausible contributor too.

**What attempt 2 does differently**: reruns the identical `job.inp`
(`Calc_Hess true` means every invocation computes its own fresh Hessian, no
input change needed) but seeded from this attempt's cycle-127 geometry
instead of the scan guess -- closer to the saddle already, so if the
recomputed Hessian is a better local description than the RFO-updated
approximation this run had drifted onto, it may converge instead of
plateauing again. If it stalls the same way, the QM-region/constraint fix
(freezing the outer scaffold) is next, not a third blind restart.
