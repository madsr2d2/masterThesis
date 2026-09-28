# RC, first pass: 12-atom QM region (superseded)

`job.inp` here is the reactant-complex minimisation as first run: the
**12-atom** QM region (the ketone C=O, both peroxide O's, all four
water/H2O2 hydrogens, both waters' oxygens), `ddCOSMO(Water)`, `opt freq`,
`Calc_Hess true` / `Recalc_Hess 50`, default QM2 topology builder. It
converged (ORCA terminated normally) -- but its `job.out` was later removed in
a home-directory cleanup, so only `job.property.json` survives.

Superseded because the transition state was built with an **enlarged 20-atom**
QM region (adding the carbonyl's two alpha-carbons and the two beta-ether
oxygens, so the QM1 electrophile is the real 1,3-dialkoxy ketone and not
H2C=O). A barrier must have reactant and TS at the same QM region, so the
parent `../job.inp` re-runs the RC at 20 atoms (and `AutoFF_QM2_Method GFNFF`,
the tier convention). See `COMPUTATIONAL.md` 2026-09-19.
