# Attempt 1: aborted in multiscale setup, default QM2 topology builder

`job.inp` here is exactly what ran (ORCA 6.1.1, `! QM/XTB r2SCAN-3c
ddCOSMO(Water) OptTS Freq TightSCF`, `nprocs 8`, `maxcore 4000`,
`AutoFF_QM2_Method` left at its default `XTB`), seeded from the converged
ALPB-level TS (`../alpb_uncorrected/job_converged_alpb.xyz`).

**What happened:** ORCA aborted seconds into multiscale setup, before any QM1
energy or gradient. The default GFN-xTB topology builder perceived a QM–QM2
bond between the carbonyl oxygen O130 and ring carbon C59 at 2.48 Å — a 1,5
non-bonded contact across the folded dialkoxy ring, not a bond — added a third
link atom on top of the two real cuts (C129–C121, C129–C125), and the crude
link-atom pre-optimisation then failed (`Optimization not converged after 3
cycles` → `CANNOT OPEN FILE job_S_Link.ORCAFF.prms.tmp`). The geometry is a
valid stationary point; the fault is in ORCA's distance-based boundary
perception, which is built once at setup.

`job.out` is the raw ORCA log; `job_S_Link_opt0.*` and the `*.prms.tmp` files
are the multiscale-setup scratch. All are homelab-local (gitignored).

**What attempt 2 does differently:** the parent `../job.inp` adds
`AutoFF_QM2_Method GFNFF`. GFN-FF's connectivity sees only the two real C–C
cuts, giving 2 link atoms and clearing setup. Topology only — the QM2 level
stays `XTB2` and its Hirshfeld embedding charges are unchanged. See
`COMPUTATIONAL.md` log 2026-09-17.
