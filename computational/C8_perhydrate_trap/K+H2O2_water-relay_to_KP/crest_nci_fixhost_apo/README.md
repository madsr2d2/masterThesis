# crest_nci_fixhost_apo -- water placement around a FROZEN host

Seed: `solvcluster.xyz`, identical to `../solv_cat_2h2o/job.xyz` (130-atom host
+ 2 cavity waters from ORCA SOLVATOR, `../solv_apo6`).

Why: the unfrozen `../crest_nci_apo` run moved the host 1.4-2.2 A from the seed
in every conformer and left **0 of 58** conformers with a water in the cavity.
Here the host is held to the seed (`$constrain atoms: 1-130`, reference
`coord.ref`, force constant **0.05**) and the bias acts on the waters only
(`$metadyn atoms: 131-136`), so the search permutes water placement.

- `fixhost.inp` is `crest solvcluster.xyz -constrain 1-130` with the force
  constant lowered from its default 0.5. **0.5 aborts the initial
  optimisation** (`Initial geometry optimization failed!`); 0.05 runs.
- Quick test of this setup (`-mquick -mdlen x0.3`, 38 s): 9 conformers, host
  0.003 A from the seed in all, both waters in the cavity in all.
- Launch: `./run.sh` (about 30-60 min on 8 threads for the full search;
  the unfrozen 136-atom run took 2 h 27 min on 4).
- Read it with `computational/crest_io.py` (`load_run`, `host_rmsd`,
  `cavity_occupancy`). A frozen host's energies are not comparable to the
  unfrozen run's: the restraint is in them and the seed is not a GFN-FF minimum.
