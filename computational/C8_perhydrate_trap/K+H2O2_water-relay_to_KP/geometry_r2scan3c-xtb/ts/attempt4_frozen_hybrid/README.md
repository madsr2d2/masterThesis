# Attempt 4: frozen host + hybrid Hessian — the workaround that reached the TS

`job.inp` here is attempt4 as run: 20-atom QM region, `AutoFF_QM2_Method GFNFF`,
`Hybrid_Hess {relay atoms}` + `HESS_Modification EV_Reverse` in place of a full
Hessian, the non-QM host frozen by Cartesian `Constraints`, and
`TS_Mode {B 135 128}` / `TS_Active_Atoms`. (Its `job.out` was removed by a
home-directory cleanup; the geometry survives as `../job.xyz` and its
`job.property.json`/`job.hess`.)

It "converged" in one cycle — the frozen host plus hybrid Hessian made
attempt3's c20 frame look like a clean stationary point, which looked hollow
until the exact full `NumFreq` in `../verify_exact_hessian/` confirmed the
geometry is a genuine first-order saddle: **one** imaginary mode at
−576.3 cm⁻¹ (`orca_io.stationary_point` = `"transition state"`). So the
workaround landed on the right TS.

Superseded by the parent `../job.inp`, which runs the same geometry as an
**unconstrained** OptTS with a full exact Hessian (no `Constraints`, no
`Hybrid_Hess`). See `COMPUTATIONAL.md` 2026-09-18.
