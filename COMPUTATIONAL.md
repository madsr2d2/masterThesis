# Computational task register

Pending, running and completed quantum-chemistry work for this thesis, with the
reasoning behind each task and a dated log of what was actually run. Companion
to `MECHANISM.md` (the chemistry), `FITTING.md` (the model fits) and
`DATA_VERIFICATION.md` (the data).

A task earns a place here when a question cannot be settled from the literature
or the existing dataset, and a calculation could settle it. Each entry states
**what it would decide** — a calculation that does not change a conclusion is
not worth running.

## Conventions

- **Status**: `PENDING` (specced, not started) · `RUNNING` · `BLOCKED` ·
  `DONE` · `ABANDONED` (with reason).
- Every task carries a **validation gate** where one exists: a known quantity
  the method must reproduce before its unknown answer is trusted. A result
  without a gate is reported as an estimate, never as a number.
- Log entries are appended newest-first under [Log](#log), with the input files
  and the ORCA version used.
- **File layout, since 2026-09-09.** `computational/` is the curated, permanent
  home for this register's calculations — `orca_stuff/` was only ever the raw
  import of what had been run on the homelab before this layout existed, and
  is not being migrated wholesale (its still-relevant work moves over task by
  task as it's revisited). Layout: `computational/C<n>_<slug>/<reaction>/`,
  where `<reaction>` names the species using `MECHANISM.md`'s own shorthand
  (`K+HOO-_to_KP`, not an ad hoc English or SMILES-like name — `orca_stuff` has
  `catOO+BnOH`, `catHO2(-)+BnOH` and `cat_HO2(-)+Bnal` side by side, close
  enough to confuse). Species geometries reused across tasks (isolated K, Kh,
  KP, HOO⁻, H₂O₂, …) live once in `computational/species/`, not copied per
  task. Each reaction folder carries its own pipeline stages as real
  subdirectories rather than overwriting one `job.inp` in place through scan →
  TS-guess → optTS → IRC — `orca_stuff`'s `catHO2(-)+BnOH/TS/` lost that
  history except for what ORCA's own incidental `.out.v###.xyz` versioning
  happened to preserve:
  - `geometry_<method>/{rc,ts,pc,irc}/` — the QM/XTB geometry/TS/Hessian tier
    (currently `b973c-xtb`, i.e. `B97-3c` QM region on GFN-xTB, both in
    `ddCOSMO(water)`; see the two-tier note below. Replaced `r2scan3c-xtb` on
    2026-09-21).
  - `energy_<method>/` — a single-point energy refinement, when the geometry
    tier's method isn't trusted for the final number (see C7/C8 note below).
  Outputs and scratch stay homelab-local, same policy as `orca_stuff/` and the
  same reasoning (`.gitignore` carries a matching pattern block) — only
  `.inp` and the small `.xyz` chains they need are tracked. `uvvis/` predates
  this convention, tracks its ORCA output in full by a separate prior
  decision, and is not covered by any of the above.
- **Two-tier method for anything using the full 130-atom catalyst**, settled
  2026-09-09 and revised 2026-09-21 (geometry method), worked example
  `K+H2O2_water-relay_to_KP` (C8 item 1). **Geometry, TS search and Hessian:
  `QM/XTB`, `B97-3c` on the reactive QM region, GFN-xTB on the rest.**
  B97-3c replaced r2SCAN-3c on 2026-09-21, after a three-arm test on the reduced
  C8 reactant — same system, QM region, solvent and seed, only the QM1 changed:

  | QM1 | cpu-s/atom | wall | imaginary modes |
  |---|---|---|---|
  | r2SCAN-3c | 32.5 | 1:15:55 | −65.33, −19.67 (a saddle) |
  | PBEh-3c | 36.9 | 1:25:01 | −39.74, −15.21 |
  | **B97-3c** | **25.7** | **0:59:50** | **none** |

  B97-3c is ~21% cheaper per atom than r2SCAN-3c **and the only arm that
  converged to a genuine minimum with no imaginary mode at all** (the r2SCAN-3c
  RC is a *saddle* at B97-3c, −65.33 cm⁻¹ — the surfaces genuinely differ). It
  is also the method the earlier oniom work used for its 23-atom region
  (`orca_stuff/cat/.../ONIOM/ts*/B97-3c_XTB/`). B97-3c is B97-D3(0)/def2-mTZVP +
  gCP + SRB. PBEh-3c was tried first on the manual's geometry recommendation
  ("performs particularly well in the optimization of geometries") but turned
  out *dearer* than r2SCAN-3c — the 42% Fock exchange eats the smaller basis.
  **The saving is ~21%, so it does NOT buy a much larger QM region**: a bond
  which forms across the QM/QM2 boundary needs the region grown, and that needs
  a genuinely cheap QM1 (HF-3c/MINIS, or xtb), not B97-3c.
  **Solvent for this tier is `ddCOSMO(water)`, not `ALPB`** — ALPB only solvates
  the QM2 (xtb) layer, leaving the QM1 reactive region in vacuum electrostatics
  regardless of what the keyword says (ORCA 6.1 manual, Multiscale Simulations:
  "If the ALPB model or CPCM-X are requested [within QM/XTB], the solvation
  effect is just included in the calculation for the large QM2 system"). Plain
  `CPCM`/`SMD` are not an option here at all — ORCA aborts at input-check for
  QM/XTB with "This is not implemented. Provide respective ALPB, ddCOSMO or
  CPCMX keyword instead." ddCOSMO is the one of those three the manual confirms
  gets cavity-charge propagation into QM1 (the "C-PCM/B" scheme it describes
  for real QM1/QM2 pairs).
  Final energies: **do not trust the geometry tier's functional for ΔG/ΔG‡** —
  density functionals (B97-3c included) carry an RMSE of roughly 5 kcal/mol
  (~20 kJ/mol) against CCSD(T)-quality barriers, at C7's own gate width and
  wider than C8's 4 kJ/mol target window. Run the higher-level method **as the
  QM1 of a second QM/XTB single point on the whole system** — the DIRECT route,
  settled 2026-09-20 and the ONLY single-point route since 2026-09-21:

       ! QM/XTB DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C ddCOSMO(Water) TightSCF

  ORCA accepts a correlated method as QM1 (6.1 manual 6.1.1: multiscale single
  points take "all kinds of available electronic structure methods as QM
  method") and embeds it electrostatically (`Embedding Scheme ...
  electrostatic`; it prints `Point charges in QM calc. from MM atoms    120`),
  so the DLPNO density is polarised by the 120 host atoms inside the
  whole-system ddCOSMO cavity. **That embedding is the whole point, and the
  single point is ALWAYS this embedded job.** On C8 item 1 the embedded route
  differs from the old isolated-fragment correction by 6.4 kJ/mol on the RC —
  the fragment simply cannot see the host. The number is
  `orca_io.composite_free_energy_eh(geom_job, T, dlpno_qm1_job, geom_job)` —
  the geometry job passed twice, because the DLPNO job is a single point with
  no thermochemistry: G comes from the geometry tier's frequencies plus an
  electronic level change. One job, ~4.5 min, no `job.QMRegion.xyz` extraction
  and so no geometry-freshness trap.

   **A DLPNO single point on an extracted fragment in `CPCM` is NOT used**
  (2026-09-21). The retired route is recorded only as history: it was
  `E(DLPNO) − E(geometry tier)` on the extracted `job.QMRegion.xyz` in
  `CPCM(water)` (`energy_dlpno-ccsdt/`), blind to the host, and it differs from
  the embedded route by 5–7 kJ/mol on the barrier. **One thing is still open,
  and is not to be assumed:** a species with no host at all (C1's small
  molecules, C9's `H₂O₂ + HPO₄²⁻`) has no QM2 to embed in, so a "QM/XTB ONIOM"
  whose whole system is the QM region reduces to plain DLPNO — whether those
  tasks get a host model, or are the one place plain DLPNO + CPCM survives, is
  undecided.
  **Resources**: benchmarked in place on the homelab (i9-14900K, 8 P-cores /
  16 E-cores, 125 GB RAM) — `nprocs 8` beat `16` and `24` on the DLPNO-CCSD(T)
  fragment job (43.5 s vs 51.5 s vs 50.8 s wall — measured on an 11-atom
  DLPNO-CCSD(T) job before the embedded route replaced it); past the physical
  P-core count, ranks either share a P-core's SMT sibling or land on a slower
  E-core, and MPI syncs to the slowest rank. `maxcore 4000` is already generous
  (the pilot's own peak use was 1178 MB/rank). A single small job cannot
  usefully take more of this machine than that — using the rest of it means
  running further reactions' jobs concurrently, not raising one job's
  `nprocs`.
- **Build the QM2 topology with `AutoFF_QM2_Method GFNFF`, not the default
  `XTB`**, settled 2026-09-17 on the same worked example. ORCA detects the
  QM–QM2 boundary with a distance-based topology builder whose default is
  GFN-xTB; on a folded 130-atom host that builder invented a bond between the
  ketone oxygen O130 and ring carbon C59 2.48 Å away — a 1,5 non-bonded contact
  across the dialkoxy ring, not a bond (a C–O bond is 1.2–1.4 Å) — and capped
  the already proton-accepting carbonyl O with a third link atom, after which
  the link-atom pre-optimisation died (`CANNOT OPEN FILE
  job_S_Link.ORCAFF.prms.tmp`) before any QM1 calculation ran. GFN-FF's
  connectivity sees only the two real cuts, C129–C121 and C129–C125 (scratch
  probe, same input otherwise: default `XTB` → `Created 3 link atoms`, same
  fatal; `GFNFF` → `Created 2 link atoms`, setup clears into the initial
  Hessian). **It changes the topology only**: the QM2 level stays `XTB2` (from
  the `! QM/XTB` line) and the embedding charges stay Hirshfeld charges off
  that GFN2-xTB calculation — ORCA reports `QM2 method … XTB2`, `Method for
  determining QM2 charges … Hirshfeld` and `AutoFF method … GFN-FF` as three
  independent settings — so charge quality is not degraded. The one knock-on is
  the boundary charge-alteration scheme (`ChargeAlteration CS`), which
  redistributes charge using the GFN-FF bond list.
- **Write the geometry and its frequency as ONE job: `! OptTS Freq` (a TS) or
  `! Opt Freq` (a minimum)**, settled 2026-09-25. The frequency runs
  automatically at the converged geometry, in the same job and the same `.out`,
  so there is no geometry file to carry between jobs and no way to point the
  frequency at the wrong structure — on 2026-09-25 a hand-copied `NumFreq`
  launched from the wrong working directory silently re-ran the *previous* step,
  and the geometry it was meant to check sat in a sibling directory untouched.
  ORCA also **refuses to run the frequency when the optimisation did not
  converge** (`As a subsequent Frequencies calculation has been requested /
  ORCA will abort at this point of the run`), a guard the two-job split does not
  have: a separate `NumFreq` will build a Hessian on whatever unconverged
  geometry it is handed. `Freq` is accepted by ORCA 6.1.1 (verified against a
  throwaway job). `Calc_Hess true` still supplies the seed Hessian `OptTS` needs
  to identify the imaginary mode and the trailing `Freq` recomputes it at the
  final geometry — **two Hessians, exactly the two-job split's cost**, so the
  win is provenance and the convergence guard, not time. `orca_io` reads the
  result from the single `.out` as it does any other. **Add `TightOpt` for
  anything a barrier is quoted from.** The default tolerance leaves the
  gradients near 1e-4 Eh/bohr and the ZPE/G are read off the Hessian there, so
  under-convergence shows up in the thermochemistry rather than in `E`.
  `TightOpt` tightens the default about three-fold (`TolRMSG 1.0e-4 → 3.0e-5`,
  `TolMAXG 3.0e-4 → 1.0e-4` Eh/bohr, confirmed in the job's own output). On the
  C8 TS it separated two consistent runs — two `OptTS Freq` jobs (default and
  `TightOpt`) agreeing to 0.09 kJ/mol on G — from the older two-job
  optimisation, which had landed in a different low-frequency well and differed
  by 30 kJ/mol on G through the ZPE alone.
- **Never trust one initial Hessian across a long optimization — use
  `Recalc_Hess`.** `Calc_Hess true` alone computes an exact numerical Hessian
  at cycle 0 and then lets the optimizer's own RFO update approximate it for
  every subsequent step; on this project's floppy 130-atom macrocycle that
  approximation drifts. First seen on a TS search (`K+peroxide+BnOH_bridged_
  UNIDENTIFIED/`, since found to be mis-scoped — see below — but the
  optimizer pathology is real regardless): plateaued 40+ cycles with no
  gradient improvement despite the negative-eigenvalue count looking healthy.
  Seen again, milder, on a plain ground-state minimization
  (`K+H2O2_water-relay_to_KP/geometry_r2scan3c-xtb/rc/`): converged steadily
  through cycle ~111 then started *regressing* (RMS/MAX gradient both got
  worse, not better) by cycle 131. Add `Recalc_Hess N` in `%geom` alongside
  `Calc_Hess true` (ORCA 6.1 manual, "4.3. Transition State Searches":
  recommended when "the PES near the TS can be very far from ideal for a
  Newton-Raphson step") — **`N=10` for a TS search, `N=50` for a plain
  minimization**; minimizations don't drift as fast (BFGS/RFO for a minimum
  is more forgiving than eigenvector-following for a saddle) and each
  recompute costs real wall time (roughly a couple of minutes on this QM
  region size at `nprocs 8`, cheap enough to afford every 10 cycles on a TS,
  wasteful more often than every 50 on a minimization). When restarting a
  run that's drifted, seed from its own latest geometry rather than from
  scratch — killing and resuming with the fix costs one Hessian recompute,
  not the cycles already made.

## Environment

| | |
|---|---|
| ORCA | 6.1.1 — `~/orca_6_1_1/orca` (in `$PATH`); `~/orca_6_0_1/orca` also present, unused |
| verified working | yes, HF/def2-SVP water single point (`computational/hellowater/`) |
| xtb / CREST / Psi4 / NWChem / Gaussian | not installed |
| pyscf / ASE / RDKit / cclib | not installed |

ORCA 6.1 includes the **ESD** module (excited-state dynamics), which computes
vibrationally resolved absorption band shapes with Franck–Condon **and**
Herzberg–Teller terms. That capability is what makes task C1 viable at all —
see the note on vibronic intensity below.

## Task register

| id | task | status | decides |
|---|---|---|---|
| [C1](#c1--extinction-coefficients-at-285-nm) | UV spectra and ε at 285 nm for the absorbing species | **PENDING** | what the absorbance actually measures — the largest open question in `MECHANISM.md` |
| [C2](#c2--criegee-adduct-pka) | pKa of the Criegee adduct at the active-site ketone | backlog | whether the pH-rate inflection can discriminate step 6's branching |
| [C3](#c3--dioxirane-closure-vs-baeyer-villiger) | TS comparison, 3-exo-tet closure vs Baeyer–Villiger | backlog | whether step 6 is energetically sane; no literature comparison exists |
| [C4](#c4--peroxide-cannizzaro-feasibility) | barrier for HOO⁻ vs OH⁻ Cannizzaro hydride transfer | backlog | whether steps 1–2, which have **no literature precedent**, are plausible at all |
| [C5](#c5--why-the-two-substrates-bend-in-opposite-directions) | substituent effect on the hydride-acceptor step, Ph vs 4-MeO-C₆H₄ | **PENDING** | whether the methoxy group is what removes BnOH's autocatalysis in the 4OMe runs |
| [C6](#c6--the-oxidant-attacks-the-product-fifty-times-faster-than-the-substrate) | ArCH₂OH vs ArCHO oxidation barriers, Ph and 4-MeO-C₆H₄ | **PENDING** | what the 4OMe progress curves' slowdown is, and why BnOH's curves do not show it |
| [C7](#c7--what-the-catalyst-is-doing-during-the-induction-period) | hydration equilibrium and dehydration barrier of the active-site ketone in water | **PENDING** | what the induction period is — the one step in the cycle with a measured barrier and no assignment |
| [C8](#c8--is-the-perhydrate-on-the-activation-path-or-off-it) | free-energy profile K + H₂O₂ ⇌ KP → KD, and the KP resting fraction | **PENDING** | whether the peroxide adduct is the catalyst waking up or the state it has to leave |
| [C9](#c9--is-the-buffer-carrying-the-peroxide) | ΔG° of phosphate perhydrate formation, and the same for pyrophosphate | **PENDING** | whether the buffer's kinetic role can be to deliver the oxygen rather than to be a base |
| [C10](#c10--where-the-catalyst-destroys-the-peroxide) | barriers for O₂ release from KP, from KP + H₂O₂, and from KD + H₂O₂ | **PENDING** | where `MECHANISM.md` S4 sits in the cycle — whether the peroxide sink shares step 4's intermediate |

---

## C1 — extinction coefficients at 285 nm

**Status: PENDING.** Specced 2026-08-30, not started.

### What it decides

`MECHANISM.md`'s structural analysis shows the mechanism predicts a
**decelerating** benzaldehyde curve and an **accelerating** total-product curve,
and the data shows the latter (52% of curves reach peak slope more than 15% into
the run). So the mechanism stands only if the absorbance is not tracking
benzaldehyde alone. The observation equation carries

```
signal(t) = [A] + r * [BA]        r = eps_BA / eps_A  at the monitoring wavelength
```

and `r` is currently a fitted parameter. Pinning it independently turns a free
parameter into a constraint, and if the fit converges on a value near the
computed one that is real corroboration of the whole picture.

### What is already known

Settled from the sheets and the literature (see `DATA_VERIFICATION.md`
2026-08-30):

| | |
|---|---|
| monitoring wavelength | **285 nm** for BnOH (43 experiments), **300 nm** for 4OMe-BnOH (**53** — the 46 that declare it, plus exps 2, 4, 5, 7, 8, 9, 10 ruled 300 nm on 2026-08-30, their sheets' 285 being a stale cell copied from the BnOH template) |
| settled 2026-08-30 | the wavelength follows the substrate and `e` is a uniform analysis convention, so all nine sheets declaring 285 nm for a 4OMe run (including exps 57/58, which also carry ε = 1.59) are working notes rather than measurements |
| what C1 could still add | no source has been found for 4-methoxybenzaldehyde at 300 nm, so the **cross-substrate ratio 7.53/1.23** rests on nothing checkable. Within a substrate a wrong ε is one global scale factor absorbed into the fitted constants; between substrates it enters any comparison of rate constants directly. Pinning both aldehydes at their own monitoring wavelengths would close that |
| benzaldehyde, water | **ε ≈ 1400 M⁻¹cm⁻¹ at 278–279 nm** (weak n→π*; the strong π→π* is at 248 nm, ε ≈ 12,000–14,000) |
| the sheets' `e = 1.23 mM⁻¹cm⁻¹` | = 1230 M⁻¹cm⁻¹ at 285 nm, on that band's falling edge — so `e` is benzaldehyde's own ε, **not** a differential coefficient |
| benzoate | strong band 224–230 nm; weak band near 268 nm. **No reliable aqueous ε at 285 nm found** — the two best sources are paywalled (HTTP 403) |

Band-shape bracket only: ε(benzoate, 285 nm) plausibly 100–400 M⁻¹cm⁻¹, i.e.
`r` ≈ 0.08–0.33. That is an estimate, not a result.

### Why the obvious calculation would be wrong

**Do not run a plain TD-DFT single point and broaden it.** Both bands of
interest are weak, symmetry-restricted transitions:

- benzaldehyde 285 nm — carbonyl n→π*, formally allowed but weak;
- benzoate ~268 nm — the benzene-derived **¹L_b** band, **symmetry-forbidden**,
  drawing essentially all its intensity from **vibronic (Herzberg–Teller)
  borrowing**.

A Franck–Condon oscillator strength for a vibronically allowed band is ≈ 0, so
the entire quantity of interest is exactly what a vertical calculation omits.
Compounding it, ε is wanted **on the tail**, 17 nm off the benzoate peak, where
intensity is set by band shape rather than by the vertical transition — and
TD-DFT vertical energies routinely err by 0.2–0.3 eV, about 14 nm here, worth a
factor of 2–3 in ε on a steep edge. A naive calculation returns a confident
number that means nothing.

### Plan

Species: **benzaldehyde**, **benzoate anion**, **benzyl alcohol**,
**perbenzoic acid**. (Benzyl alcohol matters because at 8–10 mM it is 40–80×
the product concentration, so even ε ≈ 20 M⁻¹cm⁻¹ contributes like 1 mM of
aldehyde. At 0.6% conversion it is a constant baseline offset rather than a
kinetic term, but it should be quantified rather than assumed away. Perbenzoic
acid is a state in the reduced model and has never been considered as an
absorber.) The 4-methoxy analogues need the same treatment at **300 nm**.

| stage | work | effort | gate |
|---|---|---|---|
| 1 | ground-state opt + freq, then vertical TD-DFT, CPCM(water), all four species | ~1 h | band positions and state characters match the known assignments |
| 2 | ESD **FC + HT** vibronic spectrum for **benzaldehyde only** | ~half a day | **must reproduce ε ≈ 1400 at 278–279 nm within ~30%, or stop** |
| 3 | same protocol for benzoate; read ε(285 nm); form `r` | ~half a day | only meaningful if stage 2 passed |

The gate at stage 2 is the point of the whole design. Without it the benzoate
number has no error bar; with it the claim becomes "this protocol reproduces a
known ε for the closely related chromophore, so its benzoate value is trusted to
about the same factor" — which is defensible in writing.

Proposed level: `wB97X-D4/def2-TZVP`, `CPCM(water)`, TightSCF. Molecules are
14–15 atoms, so cost is not the constraint. Stage-1 input shape:

```
! wB97X-D4 def2-TZVP TightSCF Opt Freq CPCM(water)
%tddft  nroots 10  end
* xyz <charge> 1
  ...
*
```

(benzoate is charge −1; the others neutral.)

### Alternative that makes this unnecessary

**One cuvette of benzoic acid at working pH, read at 285 nm.** Ten minutes on a
spectrophotometer answers the question directly with no methodological caveats
at all. The calculation is a day plus an argument.

**Recommendation: run stage 1 regardless** — it is an hour and it answers things
the experiment does not (band assignments, and ε for benzyl alcohol and
perbenzoic acid, neither of which anyone is going to measure). Take stages 2–3
only if the spectrophotometer is not available.

---

## C5 — why the two substrates bend in opposite directions

**Status: PENDING.** Specced 2026-09-02, not started.

### What it decides

`product_fate/ANALYSIS.md` establishes that the catalysed 4OMe-BnOH
curves slow **in proportion to the product they have made** — −0.919 ± 0.161 on
the comparison that holds every condition fixed — while the catalysed BnOH
curves do not, at product concentrations nearly twice as high. The two
substrates' progress curves bend in opposite directions, and the difference is
5.8σ.

The conjecture this task tests is that **one substituent turns both signs**:

- steps 1–2 of `MECHANISM.md` make the autocatalysis run through the product
  aldehyde acting as the **hydride acceptor** (`C1 + A → PBA + S`). Accepting a
  hydride at the carbonyl is favoured by electron withdrawal, ρ > 0, so a
  4-methoxy group should **shut that step down** — no autocatalysis;
- attack by an electrophilic oxidant on the ring or the aldehyde is favoured by
  electron donation, ρ < 0 on σ⁺, so the same group should **switch a
  competing consumption of the product on** — the slowdown (C6).

If the calculation gives those two signs with a big enough separation, the
substituent explains the whole contrast and steps 1–2 gain their first
independent support. If step 2's barrier is insensitive to the substituent, the
autocatalysis is not running through the product as hydride acceptor and
`MECHANISM.md`'s step 2 needs rewriting.

### What to compute

Step 2's transition state, `C1 + A → PBA + S`, twice: Ar = phenyl and
Ar = 4-methoxyphenyl on **both** partners, since both the tetrahedral adduct and
the hydride acceptor carry the ring. Same functional, same basis, same implicit
solvation, same conformer search protocol, and report the **difference**, not
the two absolute barriers.

### Validation gate

The classical hydroxide Cannizzaro has a measured Hammett ρ. The method must
reproduce its **sign and rough magnitude** on the OH⁻ reaction before its
answer on the HOO⁻ reaction is quoted. *No source has been read for that ρ yet
— find one before running anything, and record it here.* Without the gate this
is an estimate, per the conventions above.

### Relationship to C4

C4 asks whether steps 1–2 are feasible **at all** (HOO⁻ against OH⁻ as the
initiating nucleophile). C5 asks whether they are **substituent-sensitive in the
direction the kinetics require**. C4's geometries are most of C5's work, so run
C4 first; C5 is then two more substituted analogues on the same protocol.

---

## C6 — the oxidant attacks the product fifty times faster than the substrate

**Status: PENDING.** Specced 2026-09-02, not started.

### What it decides

`product_fate/ANALYSIS.md` §3 shows the 4OMe slowdown has the form `A′ = v − kA`: the rate falls
**linearly** in the accumulated product, on 24 of 29 curves against 0 for the
hyperbolic form that reversible product inhibition would give. That is
production minus a first-order loss of the measured species — the oxidant
attacking the aldehyde it has just made, either destroying it or being diverted
from the alcohol by it. **Absorbance cannot tell those two apart**; they are the
same reaction seen from two ends. A calculation can.

### The gate is already measured, which is what makes this worth running

The stationary level `A∞ = v(S)/k` gives a selectivity directly:

| | |
|---|---|
| k<sub>A</sub>/k<sub>S</sub>, 4OMe-BnOH | **median 54, IQR 42–81** (an upper bound) |
| as a barrier difference at 25 °C | **≈ 9.9 kJ/mol**, i.e. 2.4 kcal/mol |
| the same quantity for BnOH | **not resolvable**: 0.386 mM of benzaldehyde produces no measurable product-driven deceleration, so k<sub>A</sub>/k<sub>S</sub> is small enough to hide |
| independent corroboration | the plateau carries substrate order **+0.610 ± 0.067**, against **+0.577** measured on the rates before this was looked for |

So the calculation has a number to hit and a sign to get right, and the two are
independent of each other.

### What to compute

ΔG‡ for the oxidant attacking

1. the benzylic C–H of ArCH₂OH — the productive step 7, `KD + S → K + A`
2. the aldehyde of ArCHO — Baeyer–Villiger addition at the carbonyl, and
   ring/side-chain attack if the surface offers one

for **Ar = phenyl and Ar = 4-methoxyphenyl**, at one level, in one solvation
model. Four transition states, and what is quoted is the pair of differences.

Whether the oxidant modelled is the dioxirane (`KD`) or the peracid (`PBA`) is a
second axis; start with the dioxirane, which is step 7's own oxidant, and add
the peracid only if the dioxirane answer disagrees with 9.9 kJ/mol.

### Expected result, stated in advance so the run cannot be read backwards

ΔΔG‡(4-MeO) ≈ +10 kJ/mol in favour of attacking the aldehyde, and
ΔΔG‡(H) at least ~10 kJ/mol smaller. If instead the two substrates come out the
same, the slowdown is not a substituent effect on the oxidation and §6's
substrate contrast needs a different explanation — the most likely alternative
being that it belongs to the buffer, since the BnOH set is largely
pyrophosphate and this block is phosphate.

### What it would also give the observation equation

If the 4-methoxy aldehyde is being consumed, the products of that consumption
absorb at 300 nm too — 4-methoxyphenol and its formate, on the
Baeyer–Villiger route. Their ε at 300 nm belongs in **C1**, which currently
scopes only the aldehydes and acids. Fold them in when C1 is run.

---

## C7 — what the catalyst is doing during the induction period

**Status: PENDING.** Specced 2026-09-02, not started.

### What it decides

`induction/ANALYSIS.md` measures the induction of the catalysed 4OMe curves and
gets further than expected without a mechanism, which is exactly the position
that makes a calculation worth running. What is established:

| | |
|---|---|
| it needs the catalyst | 0 of 49 enzyme-free 4OMe curves have one, at matched composition, including a 17934 s run |
| it ends on a clock, not at a product threshold | **−0.025 ± 0.109** on 147 curves, where product control requires −1 |
| its amplitude is a fraction, not a concentration | substrate order **−0.114 ± 0.169** |
| it is 126× faster than turnover | ΔG‡ gap **−11.99 ± 0.62 kJ/mol** at 298 K |
| **its barrier** | **E<sub>a</sub> = 95.0 ± 15.7 kJ/mol**, ΔH‡ = 92.6 ± 15.7 |

The last row is the gate. **The barrier is three to four times too large for
dissolution, diffusion, de-aggregation or a cuvette reaching temperature** —
those run at 15–25 — so the catalyst is making or breaking a bond before it can
turn over. The row above is the one-phase fit's τ over four temperatures;
`induction/ANALYSIS.md` §7f measures the same barrier on the window-free clock
over all six and gets **77 ± 12 kJ/mol**, and the gate is met either way.

**C7 now has two more numbers to hit, and they are signs rather than sizes.**
The induction is *slowed by base* (`d ln τ/d pH = +0.16 to +0.33`, four ladders,
three buffers, both substrates) and *hurried by the substrate* where the buffer
is held fixed (`d ln τ/d ln[S] = −0.16 to −0.62`). A gem-diol hydrate whose
conjugate base cannot dehydrate reproduces the first; a hydrophobic guest in the
cavity shifting the hydration equilibrium reproduces the second. **So compute
the hydrate's pK<sub>a</sub> as well as its ΔG° and barrier**: near 10 puts the
diolate exactly where the boric ladder (pH 8.46–10.34) sits and the scheme
stands; near 13 rules the diolate out and leaves C8's perhydrate trap as the
only survivor.
The concentration orders say the step is unimolecular in everything the cuvette
holds. What is left is a unimolecular change on the catalyst itself, and the
leading candidate has a name.

### The candidate

**The active-site ketone in water is largely its gem-diol hydrate, and only the
free ketone can add H₂O₂.** Hydration is fast and the equilibrium can lie a long
way towards the diol for an electronically activated ketone; **dehydration is
unimolecular**, general acid/base catalysed, and can be slow. That reproduces
every measured feature at once: it needs the catalyst, it is first order in
nothing else, its amplitude is the equilibrium hydrate fraction (a fraction, not
a concentration) and its barrier is a covalent one.

It also makes a prediction the archive cannot test and this calculation can: the
hydrate fraction at equilibrium has to be **large enough to be the measured
amplitude**. The depth of the induction is 0.79 at 15 °C and 0.06 at 40 °C, and
those are *lower bounds* — the rolling window that reads them understates the
amplitude, and understates it most where the induction is fastest. So the
computed hydrate fraction must be at least ~0.8 at 15 °C, and must fall with
temperature.

### What to compute

For the chemzyme's active-site ketone (and, as a calibration, for acetone and
for one α-halo ketone whose K<sub>hyd</sub> is known experimentally):

1. ΔG of hydration, `K=O + H₂O ⇌ K(OH)₂`, in implicit water with explicit
   waters in the first shell. Report K<sub>hyd</sub> at 288 and 313 K, which
   brackets the block.
2. ΔG‡ for the **dehydration** `K(OH)₂ → K=O + H₂O`, water-assisted and
   general-base assisted (one phosphate dianion in the model — the block is
   50–80 mM phosphate at pH 7.00).
3. ΔG‡ for `K=O + H₂O₂ → KP` at the same level, for comparison with (2).

### Validation gate, stated in advance

- ΔG‡ for (2) within about **±20 kJ/mol of 92.6** — the measurement's own error
  is ±15.7 and no solvated barrier is worth more than that.
- K<sub>hyd</sub> at 288 K corresponding to a hydrate fraction **≥ 0.8**, falling
  towards 313 K.
- (3) **smaller** than (2), or the assignment is wrong: if adding peroxide is
  the slow step, the induction should have carried an order in `[H₂O₂]`, and
  §4a's inability to measure one becomes the whole story rather than a caveat.

Failing the first two does not merely weaken this candidate, it eliminates it,
and the next ones in line are the perhydrate collapse (**C8**) and a
conformational change of the cyclodextrin — the last of which is not a
calculation this project can do credibly.

### The experiment that would replace all of this

**Pre-incubate the catalyst with H₂O₂, then add substrate.** If the induction is
catalyst activation it disappears; if it waits for product it does not. One
cuvette, and it settles in an afternoon what C7 and C8 together only make
plausible. Recorded here because the archive is closed and it is the single most
valuable measurement that was never made.

---

## C8 — is the perhydrate on the activation path, or off it?

**Status: PENDING.** Specced 2026-09-02, not started.

### What it decides

`MECHANISM.md` step 4 makes `K + H₂O₂ ⇌ KP` the catalyst's entry into the cycle,
and the natural reading of an induction period is that this is what is being
timed. **The archive's one attempt to test that gets the wrong sign**: in exps
127–131, the only 4OMe cuvettes that move `[H₂O₂]` (3.879 against 195.882 mM),
the induction *time* has an order of **+0.302 ± 0.092** in peroxide. A
relaxation towards `K + H₂O₂ ⇌ KP` goes as `1/τ = k_f[H₂O₂] + k_r` and therefore
has an order between 0 and −1; it cannot be positive.

That result is **not usable as it stands** — the same block's induction
statistic regresses on signal-to-noise at +0.702 ± 0.241 and `[H₂O₂]` is what
sets the signal, so 15 curves cannot separate the two. But the sign it points
at is a coherent mechanism, and it inverts the role of step 4: **the perhydrate
would be a resting state the catalyst has to leave**, not the activation itself.
More peroxide would then mean more catalyst parked as KP and a *longer* wait
before enough dioxirane exists to turn over.

### The gate, added 2026-09-02

`induction/ANALYSIS.md` §4b sharpened this into something a calculation can be
scored against. The scheme constrains **both** peroxide orders at once — taking
log-log slopes in h,

    d ln v / d ln h = 1/(1 + Kh)     and     d ln τ / d ln h = −Kh/(1 + Kh)

so their difference is **1 identically, for every K and every h**. Measured as
one regression on `log(v/t_ind)`, it is +0.502 ± 0.194 on exps 127–131 (2.6σ
short) and +0.304 ± 0.188 on exps 135–151 (3.7σ short), and it stays short at
every log floor from 1 s to 300 s. Separately, **the rate is not first order in
peroxide either**: on the 63-curve two-axis ladder over 2.45–163 mM the free
power is a = 0.654 and a = 1 is rejected at **F = 32**.

**The shortfall is smaller through a clock that carries no window (2026-09-04),
and that weakens this task's premise without removing it.** `t_ind` is a rolling
window a tenth of the run wide, and on exps 135–151 `signal_control` fails, so
the landmark there is partly measuring the spectrophotometer. `tau` and
`tau_slow` come from the progress fit instead and are subject to the same
identity; asked of the O₂-corrected curves they fall **2.0σ and 0.3σ** short of
+1 rather than 3.7, with the substrate control missing by 3.9–8.5σ as it must
(`induction.joint_clocks`). Two consequences for C8. The measured shortfall that
motivates the trap reading is **route-dependent**, so a calculation that lands
on a modest K no longer has a 3.7σ discrepancy to explain, only a 0.3σ one. And
`tau_slow` is resolved on 34 of 110 curves with the estimate moving +0.87 to
+1.26 across cuts, straddling +1 rather than falling short of it, so this is not yet strong enough to overturn the trap reading
either. **Compute the profile; the measurement will not settle it alone.**

Inverting the trap form `d ln τ/d ln h = +Kh/(1 + Kh)` for K, and reading the
same constant off the rates through the saturating fit:

| route | K | ΔG° |
|---|---|---|
| 4OMe induction order, at 40.8 mM | 11 M⁻¹ | **−5.85 kJ/mol** |
| BnOH induction order, at 28.3 mM | 29 M⁻¹ | **−8.31 kJ/mol** |
| BnOH rates, profiled | 29 M⁻¹ (13–54) | −8.36 (−6.43 to −9.90) |

**So C8 now has a number to hit: ΔG° of perhydrate formation between −6 and
−10 kJ/mol**, i.e. K = 11–54 M⁻¹ and 50–80% of catalyst as KP at the archive's
working 82.5 mM. A result in that window corroborates the trap from a direction
with no spectrophotometer in it. A strongly negative ΔG° — a perhydrate that
dominates at millimolar peroxide — contradicts all three routes at once, and
then the positive induction order is signal-to-noise after all and step 4 keeps
its place.

Treat the three-route agreement as a target and not as support: two of the three
come from blocks whose induction statistic tracks its own signal-to-noise, so
their agreement is also what a single shared artefact would produce.

### What to compute

The free-energy profile along `K + H₂O₂ ⇌ KP` and onward, at one level and one
solvation model:

1. ΔG and ΔG‡ for `K + H₂O₂ ⇌ KP` in both directions.
2. ΔG‡ for `KP → KD + H₂O` — direct closure of the perhydrate to the dioxirane,
   the route that would make KP an intermediate.
3. ΔG‡ for `PBA + K → KD + BA` — `MECHANISM.md`'s own step 6, the route that
   makes KP a **dead end** because the acylating agent is the peracid and not
   the peroxide.

### The other discriminating comparison

Compare (2) with (3). If (3) is much the lower, KP is off the path, the
catalyst's activation waits for the first PBA, and step 4 is a trap — which
would also explain why the induction is not first order in peroxide. If (2) is
competitive, the perhydrate is on the path and the positive peroxide order in
§4a has to be an artefact of signal-to-noise after all.

### What it needs from C7

The KP resting fraction at 82.5 mM H₂O₂, which is C7's item (1) and (3)
together. A trap only matters if the catalyst is actually in it: if
K<sub>eq</sub>[H₂O₂] ≪ 1 at 82.5 mM, KP is neither an intermediate nor a trap
and both branches of this question are moot.

---

## C9 — is the buffer carrying the peroxide?

### What it decides

`buffer/` §6 puts a scheme beside general base catalysis that the archive cannot
tell apart from it:

    H₂O₂ + P     ⇌  P–OOH                 the buffer perhydrate
    P–OOH        ⇌  P–OO⁻ + H⁺
    K + P–OO⁻    ⇌  K(O⁻)–OO–P            the Criegee adduct at the ketone
    K(O⁻)–OO–P   →  KD + P–O⁻             the dioxirane, and the buffer back

It is the same shape as `MECHANISM.md`'s autocatalytic step with the buffer
where the peracid stands — now written into the mechanism itself as **step 6b**
(added 2026-09-08) — and it has a reason that is not kinetic: closing a
dioxirane from a Criegee adduct made of plain H₂O₂ means expelling **hydroxide**,
and that is why the same chemistry elsewhere runs on peroxymonosulfate or a
peracid instead. A buffer perhydrate would supply the leaving group the reaction
needs, at the only concentration in the flask large enough to matter.

**The kinetics cannot decide it.** A general base is a term in `[buf]`; a buffer
perhydrate is a term in `[buf][H₂O₂]`. At one pH they differ by that interaction
and by nothing else, and **0 of the archive's 88 runs step both** — all five
titrations sit at 82.5 mM peroxide. The two-pH species test does not separate
them either: nucleophilic attack by the phosphate dianion carries the same pH
signature as general base catalysis by the phosphate dianion, so this scheme is
a fourth survivor of `buffer/` §3 rather than a resolution of it.

So the discrimination has to come from a calculation or from a new experiment,
and the calculation is small.

### What to compute

Aqueous ΔG°, one level and one solvation model, for:

1. `H₂O₂ + HPO₄²⁻ ⇌ HOO–PO₃²⁻ + H₂O` — the perhydrate equilibrium itself.
2. The same for **pyrophosphate**, `H₂O₂ + HP₂O₇³⁻ ⇌ HOO–PO₃²⁻ + HPO₄²⁻`.
   Perhydrolysis needs an electrophilic phosphorus: orthophosphate has none and
   pyrophosphate has a leaving group, so the two should differ by a great deal.
3. pKₐ of the peroxo proton in `HOO–PO₃²⁻`, which decides how much of the
   perhydrate is the reactive anion at pH 6.5–7.5.
4. ΔG‡ for `K + HOO–PO₃²⁻ → KD + HPO₄²⁻`, against C8's item (2) — the same
   closure driven by H₂O₂ alone, where hydroxide is the leaving group. The
   comparison is the whole argument: if the phosphate route is not far lower,
   the scheme has no reason to exist.

### The gate

Item (1) has a known answer in sign and roughly in size: **peroxomonophosphate
hydrolyses**, it is made from P₂O₅ or POCl₃ with H₂O₂ and not from phosphate in
water, so ΔG° must come out **positive**. A method that returns a favourable
perhydrate for orthophosphate has failed its gate and items (2)–(4) are not to
be read.

A positive ΔG° does **not** kill the scheme — catalysis through a reactive minor
species is ordinary chemistry, and that is exactly what a small equilibrium
fraction with a large rate constant looks like. What would kill it is a ΔG° so
positive that no accessible fraction exists at 200 mM buffer: at +30 kJ/mol the
perhydrate is 6 × 10⁻⁶ of the buffer, and item (4) would then have to beat the
direct route by that factor to compensate.

### The experiment that would replace it

Add the peroxide axis to the pH 6.5 titration `buffer/` §3 already asks for: a
**buffer × peroxide block at one pH**, four rungs each. The interaction term is
the question, and no run in this archive has one. The same block would test the
prediction that **pyrophosphate ≫ phosphate** as a buffer catalyst at matched pH
and matched peroxide, which `buffer/` §5 shows the archive cannot ask because
each buffer was used at its own pH range and the two candidate cells share no
peroxide value.

---

## C10 — where the catalyst destroys the peroxide

**Status: PENDING.** Specced 2026-09-04, not started.

### What it decides

`MECHANISM.md` S4 establishes that something in the enzyme-containing cuvette
catalyses `2 H₂O₂ → 2 H₂O + O₂`. The evidence is kinetic: the gas forms in the
sample beam and not the reference, which omits only the enzyme (122 falls
beyond 20σ against 23 rises); its production rate is **+1.477 ± 0.258 in
[H₂O₂] against −0.344 ± 0.093 in [S]** from a fit that never saw a
concentration; and across the whole archive it appears with both substrates in
three buffers, with a hard floor of zero detachments below pH 7.5. It is silent
about *where*, and the archive cannot become less silent: the rate is first
order in peroxide under every branch below, and **no run moves the productive
and unproductive routes against each other**.

**This task assumes the ketone is the catalyst, and that assumption is not
measured.** The cyclodextrin scaffold or a trace transition metal carried in
with the enzyme stock would each give a catalase-like decomposition with the
same peroxide order and the same pH dependence, and the archive holds no run
with cyclodextrin alone, ketone alone, or a chelator added
(`MECHANISM.md`, open questions). If the sink turns out to be metal-catalysed
it is **not** in the catalytic cycle at all and none of the three branches
below applies — so a cheap EDTA control is worth more than this calculation
and should precede it. What the calculation is still worth in the meantime:
branch barriers that all come out prohibitive would themselves be evidence for
the metal reading.

The three branches differ in what they predict for the rest of the cycle:

| branch | consequence if true |
|---|---|
| `KP → K + O₂ + H₂O` | the sink shares step 4's pre-equilibrium and saturation, competes with step 5 for the same KP, and cannot be suppressed without suppressing the catalysis |
| `KP + H₂O₂ → K + O₂ + 2 H₂O` | second order in peroxide overall; the sink is tunable by dilution while the productive route is not |
| `KD + H₂O₂ → K + O₂ + H₂O` | the dioxirane is reduced before it can oxidise; the sink competes with step 7 and scales with turnover, not with the resting state |

The third is the one that would matter most for the thesis's headline: it would
make the peroxide sink a **measure of wasted dioxirane**, and the observed
−0.307 substrate order would follow directly, since substrate intercepts KD
before peroxide can.

### What to compute

At one level and one solvation model, consistent with C7 and C8 so the numbers
compose:

1. ΔG‡ for `KP → K + O₂ + H₂O` (retro-perhydration with O–O scission).
2. ΔG‡ for `KP + H₂O₂ → K + O₂ + 2 H₂O`.
3. ΔG‡ for `KD + H₂O₂ → K + O₂ + H₂O`.
4. For scale, ΔG‡ for step 5 (`KP + S → K + A`) and step 7 (`KD + S → K + A`)
   at the same level — the sink only matters relative to the route it steals
   from, and C6 already needs (4)'s second half.

### How it would be scored

The measured gas rate is **first order in peroxide and weakly negative in
substrate**. Branch 1 predicts an order approaching zero in peroxide once KP
saturates, which the two-axis ladder's own saturation (`a = 0.654`, C8) says is
partly reached — so branch 1 sits least comfortably with +1.389 unless the
perhydrate is far from saturated. Branches 2 and 3 both give first order
cleanly. Branch 3 additionally predicts the negative substrate order without a
further assumption. **A calculation that puts (3) well below (1) and (2) would
be the first mechanistic account of the substrate order in the gas rate.**

### What it does not need

An enzyme-free comparison from the archive. There is one in principle — 16
enzyme-free curves at matched peroxide and pH, none of which bubble — but it is
worth **about one expected event** against its own buffers' catalysed rates
(`scope.gas_enzyme_control`), so it decides nothing and no branch here should
be scored against it. The catalyst-dependence rests on the beam asymmetry
instead, and that is settled.

---

## Backlog

Not yet specced. Each is grounded in an open question in `MECHANISM.md`.

### C2 — Criegee adduct pKa

`MECHANISM.md` argues step 6's dioxirane-vs-Baeyer–Villiger branching is
pH-controlled, but notes that the pH data cannot discriminate the two
mechanisms because HOO⁻ addition (pKa 11.6) produces the same rise across
pH 5.5–11.8. **The discriminator is where the inflection sits**: near 11.6
points to peroxyanion nucleophilicity, well below 10 to adduct ionization. The
document's guess that an electron-poor sugar ketone puts the adduct pKa in the
8–10 range is explicitly flagged as reasoning, not sourced — and no pKa has ever
been measured for the relevant adduct.

A computed pKa would turn that guess into a number and tell us whether the
dataset's pH window can decide the question at all. Needs a thermodynamic cycle
with a proton-affinity reference and careful solvation; absolute pKa accuracy of
1–2 units is realistic, which is enough to distinguish "8–10" from "≈11.6".

### C3 — dioxirane closure vs Baeyer–Villiger

Step 6 (`PBA + K → KD + BA`) is proposed as a 3-exo-tet intramolecular
displacement. `MECHANISM.md` records that the closest thing in existence to a
comparison is the Supporting Information of Porter/Yin/Pratt (JACS 2000,
reference 33), which is paywalled and unretrieved. A relaxed-surface-scan or
NEB-TS comparison of the two collapse channels from the same Criegee adduct
would supply the missing comparison directly, and would also test the
leaving-group penalty argument (benzoate pKa 4.2 vs bisulfate ≈ 2.0, Brønsted
β_lg −0.5 to −1.0 → 10²–10⁴ slower, offsettable by effective molarity).

### C4 — peroxide-Cannizzaro feasibility

Steps 1–2 are the thesis's own hypothesis and a dedicated literature search
found **no precedent anywhere** for an aromatic aldehyde + HOO⁻ Cannizzaro-type
disproportionation. The classical OH⁻ reaction is textbook. Computing both
barriers at the same level — HOO⁻ vs OH⁻ as the initiating nucleophile, same
aldehyde, same solvation — gives a relative number that is far more trustworthy
than either absolute barrier, and directly addresses the objection recorded in
`MECHANISM.md` (that the peroxide version is expected to be slower but is not
thereby ruled out). If the computed penalty is a few kcal/mol the hypothesis
survives; if it is 15+ kcal/mol the step is dead and the fitting strategy has to
change.

This is the highest-value backlog item scientifically, because it is the one
step in the mechanism with no external support of any kind.

---

## Log

Newest first. Record the ORCA version, the input files, the wall time and the
outcome — including failed and abandoned runs, which are the ones most easily
forgotten and most expensive to repeat.

### 2026-09-25 — QM/MM charge goes on the `*xyzfile` line; `Charge_total` is the whole system

The anion model (ketone + HOO⁻, no water) aborted instantly:

    Charge of total system ... -1
    Error : multiplicity (1) is odd and number of electrons (93) is odd -> impossible

An ORCA QM/MM job takes the charge in **two** places and they are not the same
thing:

- `* xyzfile <charge> <mult>` — the **QM region's** charge and multiplicity.
  This is what the SCF uses.
- `%QMMM Charge_total <q>` — the **whole system's** charge.

The input had `Charge_total -1` (correct: the total is −1 because the GFN-FF
host is neutral) but `* xyzfile 0 1`, so the SCF was set up for a **neutral** QM
region — 17 QM atoms + 2 link H = 93 electrons, odd, against a singlet. The fix
is `* xyzfile -1 1`, giving 94 electrons (even), confirmed by
`N(Alpha) : 47.000…`. The split is then −1 (QM) + 0 (XTB) = −1 (total), matching
`Charge_total`. Any charged QM region in a QM/XTB job needs both lines set
consistently; with a neutral QM region (every earlier C8 job) the two coincide
and the distinction is invisible.

### 2026-09-25 — IRC from the 21-atom bridge TS: it IS the perhydrate reaction coordinate

`irc_bridge/` — `! QM/XTB B97-3c ddCOSMO(Water) IRC TightSCF`, 21-atom QM1,
`Direction both`, `MaxIter 50` — ran 1 h 06 m and terminated normally, but hit
MaxIter in both directions, so the endpoints are **directions, not converged
minima**. The connectivity is unambiguous:

- **Forward → the reactant**: C128···O133 stretches 1.607 → 1.918 Å (no bond),
  the carbonyl returns to 1.250 Å, the proton stays on O133 (1.014 Å).
- **Backward → the product**: C128–O133 forms (1.607 → 1.423 Å), the carbonyl
  stretches to 1.335 Å, and the proton relays OFF O133 (1.112 → 2.018 Å) onto
  the water O131 (O131–H132 0.951 Å) — i.e. to hydronium.

So the −490.3 cm⁻¹ mode **is** the reaction coordinate: the saddle connects the
ketone + H₂O₂ pre-complex to the perhydrate, C–O formation coupled to the
proton relay. This is the confirmation a mode tally cannot give, and it resolves
the −313 (17-atom) vs −490 (21-atom) vs −576 (r2SCAN-3c) disagreement in favour
of the enlarged region: the 17-atom saddle was a cut-ring artefact.

Two notes for what follows. (1) The backward endpoint is product-like
(C128–O133 1.423 Å) and is the seed for the PC; the forward endpoint is
reactant-like. (2) Neither reached its minimum inside `MaxIter 50`, so both need
a full relaxation (21-atom QM1, `Opt Freq TightOpt`) before any barrier — the
`−9.81 cm⁻¹` mode under the cutoff on the TS also means the saddle is soft and
should be rechecked after the RC/PC are in hand.

### 2026-09-25 — C8 item 1 at the B97-3c tier: full-system TS, and the composite ΔG‡

**The full-system TS converged.** `ts/optts_full/` — `! OptTS Freq TightOpt`, no
`ActiveAtoms` (all 137 atoms), seeded from the frozen-host TS — converged in
79 cycles, 50 min 43 s, and the Freq ran automatically: `transition state`, one
imaginary mode −313.23 cm⁻¹, E(QM/QM2) = −797.138622 Eh, G(298.15) =
−796.178483 Eh. It converged despite a long step-oscillation (gradients inside
tolerance by cycle 31, steps not shrinking until the end), so the frozen host
was **helpful, not necessary**.

**`ActiveAtoms` was scaffolding and is gone from the final geometry.** It took
the optimiser from a poor full-system start to a frozen-host saddle — enough to
get a good geometry — but a barrier must be one model, and the RC was already a
full 137-atom relaxation. The frozen-host runs (`ts/optts_freq_tight/`, and the
superseded `ts/optts/` + `ts/optts_numfreq/`) are kept as the seed and the
comparison.

**Composite ΔG‡ (embedded DLPNO on the consistent full-system B97-3c
geometries):** `energy_qm1_dlpno_b973c/{ts,rc}/` run
`! QM/XTB DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C ddCOSMO(Water) TightSCF`
(NormalPNO), the 17-atom QM1 electrostatically embedded in the 120 host atoms
(`Embedding Scheme ... electrostatic`, `Point charges in QM calc. from MM
atoms 120`).

- E(QM/QM2): TS −796.589645, RC −796.621965 Eh; B97-3c: −797.138622 / −797.161042 Eh
- level correction DLPNO−B97-3c: +1441.3 (TS) / +1415.3 (RC) kJ/mol — it raises the barrier by **26.0 kJ/mol**
- **ΔG‡ = 97.2 kJ/mol (23.2 kcal/mol)**, ΔE_composite = 84.9 kJ/mol

**This is a large move from the superseded 59.5 kJ/mol** (r2SCAN-3c geometry,
11-atom QM1, embedded). Three things changed at once — geometry tier (r2SCAN-3c
→ B97-3c), QM region (11 → 17 atoms), and full-system vs reduced thermochemistry
— so 97.2 is the current best value, not a settled one. Two caveats: the full
Hessian's softest host modes are 13.5–31 cm⁻¹, so G (and hence ΔG‡) carries real
quasi-RRHO uncertainty; and −313 cm⁻¹ is far from the r2SCAN-3c TS's −576, so an
IRC is owed before this is confirmed as the same saddle.

### 2026-09-25 — the B97-3c C8 TS exists; and write the geometry and its frequency as ONE job

**The B97-3c reduced-system TS is verified.** Seed = frame 0 of
`ts/optts/job_trj.xyz` (the verified r2SCAN-3c TS with its host already relaxed
at B97-3c). `! QM/XTB B97-3c ddCOSMO(Water) OptTS` with the host held by
`%qmmm ActiveAtoms` (17 atoms) converged in **46 cycles, 6 min 41 s** to
E(QM/QM2) = −797.138648 Eh; the active-region Hessian (49 displacements, 4 min
03 s) gives **`transition state` — one imaginary mode at −379.28 cm⁻¹**, the
lowest real mode 48.96. The mode is the relay: the transferring proton
(full-system 132, |v| = 0.63), then the carbonyl C129 (0.40) and the peroxide O
(0.30). **−379 cm⁻¹ is well below the r2SCAN-3c TS's −576 cm⁻¹** — same
character, different point (B97-3c moves the TS to C129···O136 1.736 Å,
O136–H135 1.083 Å) — so an IRC is still owed before this is quoted as *the* C8
TS. The low-level QM/XTB barrier is 58.8 kJ/mol and **not yet a consistent
number**: the reduced RC (`geometry_b973c-xtb/rc`) was a full 137-atom
relaxation while this TS had its host frozen, so the two host geometries differ.

**THE OPTIMISATION AND ITS FREQUENCY GO IN ONE JOB.** `! OptTS Freq` (and
`! Opt Freq` for a minimum) runs the frequency at the converged geometry, in the
same job and the same `.out` — no geometry file to carry across, and ORCA
refuses to run the frequency if the geometry did not converge (`As a subsequent
Frequencies calculation has been requested / ORCA will abort at this point of
the run`). This is not cosmetic. On this date a hand-copied `NumFreq` was
launched from the wrong working directory and silently became a re-run of the
*previous* step; the geometry it was meant to check was never touched and
produced no `.out` at all, and nothing in a two-job split can catch that.
`Calc_Hess true` still supplies the seed Hessian `OptTS` needs and the trailing
`Freq` recomputes it — two Hessians, exactly the two-job split's cost. Recorded
in the Conventions, in the run-orca skill, and in the skill's conventions table.

**`TightOpt` is required for a barrier — and it identified a bad optimisation.**
The default `OptTS` tolerance leaves the gradients near 1e-4 Eh/bohr, and the
ZPE and G are read off the Hessian at that point, so under-convergence hides in
the thermochemistry rather than in `E`. Re-running the TS as
`! OptTS Freq TightOpt` (41 cycles, 9 min 15 s, thresholds `TolRMSG 3.0e-5` /
`TolMAXG 1.0e-4` Eh/bohr from the job's own output) gives
E = −797.138685 Eh, ZPE = 0.129941, G(298.15) = −797.042444, imaginary mode
−369.0 cm⁻¹ — **agreeing with the default-tolerance single-job run to
0.09 kJ/mol on G** (0.02 on ZPE). The two-job result is the outlier: its lowest
real modes are 49.0 and 91.1 cm⁻¹ against 130.4 and 147.2, and its G differs by
30 kJ/mol, all of it in the ZPE. So the "30 kJ/mol systematic" flagged when the
two-job number was read was not a tolerance artefact at all — it was one
optimisation that had converged to a different low-frequency well. The
`TightOpt` run is canonical; `optts_freq/` agrees with it and `optts/` +
`optts_numfreq/` are superseded.

**And the relaxed scan is not the route to a B97-3c TS.** Four constrained
attempts had already failed by the 2026-09-21 entry; the working `ActiveAtoms`
scan then died at the end of step 1 with the same `ERROR (SHARK): Failed to read
input file (job.SHARKINP.tmp)` that killed the first full-host scan (411 cycles)
— disk not full, geometry files written anyway — so the crash is in ORCA's
scan/numerical-Hessian machinery, not in the freeze. The seed was good enough
for `OptTS` without a scan, which is simpler and is what was used.

### 2026-09-21 — to freeze part of a QM/XTB system, use `%qmmm ActiveAtoms`, NOT `%geom Constraints`

Four attempts to optimise the QM region with the 120-atom host frozen failed —
251, 411, then 225 cycles without converging, and the failed PC optimisation
before them (420 cycles, ring broke). The cause was the MECHANISM, and ORCA's
own setup output says so. With `%geom Constraints { C ... C }` on the host:

    optimized atoms = activeRegion ... 137      <- ALL 137 atoms optimised
    Active atoms                   ... All atoms
    activeRegion                   ... NO
    Choice of coordinates  ... (2022) Redundant Internals

The constraints only PIN the host; ORCA still carries every atom through the
step, the Hessian and the convergence test, in a 2022-dimensional redundant-
internal set. So the free QM-region coordinates converge (RMS 2.9e-4, inside
even the tight target) while the pinned host's do not (RMS 0.16, MAX 0.69), and
the convergence criterion can never be met — no tolerance helps, because the
failing coordinates are the ones we deliberately held. Cartesian constraints
belong to a Cartesian optimisation; bolted onto a redundant-internal QM/MM run
they misbehave.

**The native tool is the QM/MM active region.** `%qmmm ActiveAtoms { ... } end`
with `ExtendActiveRegion no` gives, on the same system:

    activeRegion                   ... YES
    optimized atoms = activeRegion ... 17
    Active atoms                   ...  61  64 120 121 122 124 125 126 128 129 ...

ORCA's own QM/MM optimizer then holds everything else. The knock-on is not just
convergence: the numerical Hessian drops from **409 displacements to 49**,
because only the active atoms are differentiated — roughly an order of magnitude
cheaper for the same job. Any future "optimise the QM region with the MM region
frozen" must use `ActiveAtoms`; `%geom` constraints are the wrong layer.

### 2026-09-21 — the project's geometry method becomes B97-3c/XTB; the embedded single point becomes the only route

Two conventions changed, on the three-arm test above and the embedded-route
finding:

1. **Geometry / TS / Hessian: `! QM/XTB B97-3c ddCOSMO(Water)`** (was
   r2SCAN-3c). B97-3c was ~21% cheaper per atom, the only arm that converged to
   a clean minimum, and the method the earlier oniom work used for its 23-atom
   region. The saving is real but modest, and the Conventions now say plainly
   that it does **not** make a much larger QM region affordable.
2. **The DLPNO-CCSD(T) single point is ALWAYS embedded** — run as the QM1 of a
   second QM/XTB job on the whole system,
   `! QM/XTB DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C ddCOSMO(Water)`. A correction
   from an extracted `job.QMRegion.xyz` in `CPCM` is retired: the embedding is
   the point, and the two differ by 6.4 kJ/mol on the C8 RC and 5–7 kJ/mol on
   the barrier.

Updated for this: `COMPUTATIONAL.md` (Conventions and the layout line), the
`run-orca` skill (conventions table and the composite-energy section), and
`PLAN_C7_INDUCTION.md` (§3.1–3.2 and Task 3).

One thing is deliberately left **OPEN** rather than papered over: a host-free
species (C1's small molecules, C9's `H₂O₂ + HPO₄²⁻`) has no QM2 to embed in, so
its single point reduces to plain DLPNO + CPCM. Whether those tasks get a host
model or keep the plain route is undecided, and both documents say so.

**Consequence, NOT yet discharged:** every C8 number now on file is an
r2SCAN-3c-geometry result — the verified TS (−576.3 cm⁻¹), the `displace_minus`
RC, and ΔG‡ = 59.5–60.3 kJ/mol. The three-arm test showed the r2SCAN-3c RC is a
*saddle* at B97-3c (−65.33 cm⁻¹), so the surfaces genuinely differ and those
numbers are **superseded, not merely kept**. Re-running RC, TS and PC at B97-3c
is outstanding work, not a formality.

### 2026-09-21 — the three-arm comparison: B97-3c is the cheapest AND the cleanest, but ~20% is not enough

All three on the reduced RC seed (137 atoms, 17-atom QM region, GFN-FF,
ddCOSMO(Water)); only the QM1 method differs.

| arm | atoms | cpu_seconds | per atom | wall | imaginary modes |
|---|---|---|---|---|---|
| r2SCAN-3c (seed only) | 140 | 4546 (1.26 h) | 32.5 | 1:15:55 | −65.33, −19.67 (saddle) |
| PBEh-3c | 137 | 5051 (1.40 h) | 36.9 | 1:25:01 | −39.74, −15.21 |
| **B97-3c** | 137 | **3520 (0.98 h)** | **25.7** | **0:59:50** | **none** |

**B97-3c is the cheap one**: ~21% less CPU per atom than r2SCAN-3c and ~30% less
than PBEh-3c. It is also the only arm that converged to a genuine **minimum with
no imaginary mode at all** — no "ignoring N under 50 cm⁻¹" clause.

Geometry (old atom names; `new(old) = old − #{131,132,139 < old}`):

| pair | seed | PBEh-3c | B97-3c |
|---|---|---|---|
| C129=O130 (ketone) | 1.225 | 1.212 | 1.232 |
| O136···C129 | 2.920 | 3.832 | 3.098 |
| O136–O137 | 1.464 | 1.411 | 1.459 |
| O130···H133 | 3.125 | 2.069 | 2.052 |

RMSD: PBEh-3c↔seed 0.502 Å, B97-3c↔seed 0.348 Å, B97-3c↔PBEh-3c 0.280 Å. So
B97-3c sits between the seed and PBEh-3c, moving the peroxide off the carbonyl
far less (2.92 → 3.10, against PBEh-3c's → 3.83).

**One real consequence of the water removal, visible in all three:**
O130···H133 collapses from 3.125 to ~2.05 Å. Dropping water1 (whose H131 was the
ketone-O H-bond donor) lets water2's H133 swing in and take that role. The
reduced system is therefore not simply "the same complex minus one water" — it
re-forms the H-bond network around the remaining water.

**And the affordability answer is no.** B97-3c saves ~21% per atom, which is real
but small. Growing the QM region enough to carry the product's ring bridge would
cost a multiple of the present QM1 work, so a 21% saving does not buy it. If the
enlarged region is genuinely needed, the QM1 has to get much cheaper (HF-3c,
MINIS basis; or xtb) or the strategy has to change — B97-3c alone does not
close that gap.

### 2026-09-20 — PBEh-3c: good geometry, but NOT cheaper than r2SCAN-3c; B97-3c launched

The reduced RC (137 atoms) was optimised with `! QM/XTB PBEh-3c ddCOSMO(Water)
opt freq` — the manual's geometry recommendation. It converged (RMS 5.1e-6, MAX
5.1e-5, 70 cycles) and terminated normally in **1 h 25 m**.

**It fails the cost test.** Like-for-like against the r2SCAN-3c RC:

| arm | atoms | cpu_seconds | wall |
|---|---|---|---|
| PBEh-3c | 137 | 5051 (1.40 h) | 1:25:01 |
| r2SCAN-3c | 140 | 4546 (1.26 h) | 1:15:55 |

PBEh-3c costs ~11% MORE CPU on 3 fewer atoms — the 42% Fock exchange eats the
smaller (def2-mSVP) basis, exactly as the manual's ambiguous ordering warned
(HF-3c < B97-3c < PBEh-3c, with no clear placement of PBEh-3c against
r2SCAN-3c). **So PBEh-3c is not the lever for affording a larger QM region.**

**It converged to a genuinely different minimum.** RMSD from the seed 0.502 Å,
and the change is chemical:

| distance | seed | PBEh-3c |
|---|---|---|
| C129=O130 (ketone) | 1.225 | 1.212 |
| O136···C129 (peroxide to carbonyl) | 2.920 | **3.832** |
| O136–O137 (O–O) | 1.464 | 1.411 |
| O130···H133 (water2 to ketone O) | 3.125 | **2.069** |

The peroxide pulls ~0.9 Å off the carbonyl and water2's H133 swings in to
H-bond the ketone O. Several FRAMEWORK atoms move >1 Å (C55 1.05, O61 1.24,
H110 1.53, H117 1.41) — a lot for a rigid host, and a sign of a flat surface
with several nearby basins. Modes: two imaginary, both under the 50 cm⁻¹ cutoff
(−39.74, −15.21), so a "minimum" by the rule; the seed's −65.33 (above cutoff)
is gone, so the one real saddle was resolved.

**B97-3c launched** on the same seed, system and QM region
(`geometry_b973c-xtb/rc/`) — the remaining cheap candidate (plain GGA on
def2-mTZVP), and the method the earlier oniom work used for its 23-atom region.
It gives the third cost point and a third geometry for the comparison.

### 2026-09-20 — the spectator water identified and dropped; PBEh-3c adopted for the geometry tier

**The TS's H-bond chain settles which water is which.** Measured on the verified
TS geometry:

    O136-H135 ... O134(water2) -> O132(water1) -> O130(ketone O)
      1.129            1.307          1.897          1.822

H135 — the peroxide proton — sits 1.307 Å from O134, i.e. delivered to
**water2**; H133 (on O134) H-bonds to O132; H131 (on O132) H-bonds to the ketone
O. The IRC showed H131 never moving in either direction. So the spectator is
**water1 = O132 + H131 + H139**, and the reduced system is the RC with those
three atoms removed: **137 atoms**, QM region of 17 (0-based
`{61 64 120 121 122 124 125 126 128 129 130 131 132 133 134 135 136}`), written
to `geometry_r2scan3c-xtb/rc_reduced/job.xyz`.

**Why the reduction.** The product makes a new bond ACROSS the QM/QM2 boundary,
so the QM region has to grow to carry the complete bridge across the ring — and
at r2SCAN-3c that is not affordable. Dropping the non-participating water buys
three atoms back.

**Method: PBEh-3c**, adopted on the ORCA manual's geometry recommendation —
"highly efficient ... performing particularly well in the optimization of
geometries", and "much better geometries [than HF-3c] ... roughly of
MP2-quality". `geometry_pbeh3c-xtb/rc/` is that job: the reduced system, the
same QM region and GFN-FF topology as the r2SCAN-3c tier, so the only change is
the QM1 method.

**The A/B was started and then abandoned.** `geometry_r2scan3c-xtb/rc_reduced/`
and `geometry_b973c-xtb/rc/` hold the two arms (r2SCAN-3c reference vs B97-3c),
both killed part-way through their initial Hessians when PBEh-3c was chosen.
They are kept for the wall-time comparison if it is wanted. The cost question
they were to answer is now **open**, and the manual's ordering (HF-3c < B97-3c <
PBEh-3c, with no clear placement of PBEh-3c against r2SCAN-3c) means PBEh-3c is
not automatically the cheaper method. The wall time of `geometry_pbeh3c-xtb/rc/`
is the number that settles it.

Also recorded: the two `not implemented` warnings ORCA prints for this job are
benign and method-independent — `QMMM chosen together with analytic frequencies
/ analytical Hessian ... Switching to numerical`.

### 2026-09-20 — the PC optimisation FAILED: it never converged and the ring opened

`geometry_r2scan3c-xtb/pc/`, seeded from the IRC product endpoint
(`irc/job_IRC_B.xyz`), `opt freq` at the 20-atom region. Ran **5 h 12 min**,
terminated normally — and is **unusable**:

- **It did not converge.** It reached the maximum number of optimisation cycles
  (420) and ORCA then refused the frequency:
  `The optimization did not converge but reached the maximum number of
  optimization cycles. As a subsequent Frequencies calculation has been
  requested ORCA will abort at this point of the run.`
  So there is **no G(PC)** from it — C8's ΔG° gate and the reverse barrier stay
  blocked.
- **The geometry collapsed.** The C129–C121 bond broke during the optimisation:
  **1.549 Å in the seed → 5.822 Å** in the result, i.e. the dialkoxy ring
  opened. `job.xyz`'s C129 then has only O130 (1.170), O136 (1.301) and C125
  (1.568) within 1.9 Å. O134 carries three H at 0.78–0.79 Å (a hydronium).
- **The energy fell ~0.6 Eh** (−873.77 at cycle 1 → −874.37 at the end), far
  more than any conformational change — consistent with a rearrangement, not a
  relaxation.

So relaxing the stretched IRC endpoint unconstrained does not give the
perhydrate: the PES prefers to open the ring and move the proton onto a water.
Whether that means KP is not a minimum at this tier, or only that this seed was
too distorted, is open. Re-seeding is the next thing to try — a constrained
relaxation that holds C129–C121, or finishing the IRC properly
(`Direction down`, larger `MaxIter`) and optimising a genuine product minimum —
rather than relaxing a stretched endpoint unconstrained.

### 2026-09-20 — the TS IRC: the coordinate is the C-O bond and H135, not the water relay

`geometry_r2scan3c-xtb/irc/`, from the verified TS geometry (md5
fb47acb010589f690425efa03e53a441). Same tier and solvent as the
geometry/frequency runs. `InitHess` left at ORCA's default, so a fresh numerical
Hessian is built rather than reading the crashed `.hess`. Ran **34 min 54 s**,
terminated normally.

**It did NOT converge.** `MAXIMUM NUMBER OF ITERATIONS REACHED - STOPPING IRC
RUN` appears **twice** — both directions hit `MaxIter 20` without reaching a
minimum. So it does not yet demonstrate that the TS connects reactant to
product; it is a truncated path, not a finished one.

What it does establish, and what it corrects:

- **The coordinate is the C129-O136 bond**, not the proton relay. O136···C129
  runs 1.547 Å (product branch) ↔ 1.883 Å (reactant branch) through the TS's
  1.820 Å. On the product branch C129=O130 weakens 1.265 → 1.319 Å.
- **The proton that moves is H135**, the one on the attacking peroxide oxygen.
  It stretches off O136: 1.129 Å at the TS, 1.501 Å at the product endpoint,
  where it sits nearest O134. On the **verified** Hessian the imaginary mode
  agrees — mode 6 at −576 cm⁻¹ is dominated by **H135 (1.280 Å)**, then C129
  (0.257), H133 (0.225), H139 (0.185), O136 (0.181).
- **H131 is a spectator.** It stays on water1 (O132) at ~0.98 Å in every IRC
  frame, in both directions. **This corrects the earlier reading of this TS as a
  "Grotthuss relay through water1"** (the folder's original label): the earlier
  0.613 Å-atom-131 figure was the RC's **−313 cm⁻¹** mode — a different
  structure and a different mode, not this TS's imaginary mode.
- **The reactant branch stalled.** O136···C129 only moved 1.820 → 1.883 and the
  energy plateaued at −873.742063 Eh, which is **2.25 mEh ABOVE the TS**
  (−873.744316) — so within 20 steps the forward branch did not descend below
  the saddle at all, and the connection to the reactant is unproven. The product
  branch descended properly (23 mEh below the TS at the endpoint).

Follow-up launched: `geometry_r2scan3c-xtb/pc/`, seeded from the product-side
endpoint (`job_IRC_B.xyz`), `opt freq` at the same 20-atom region — to give the
missing G(PC). A restart of the reactant branch (`Direction down`, larger
`MaxIter`) is the documented way to finish the unconverged direction.

### 2026-09-20 — OPI-driven prototype of the same two steps, and why it is the better fit

Same cheap bench (water dimer, QM region = one water; driver
`/tmp/opencode/orca_probe/opi_workflow.py`), two steps driven from Python with
OPI's `Calculator`:

1. `! qm/xtb r2scan-3c opt freq tightscf` with
   `BlockQmmm(qmatoms=IntGroupEnd(values=[0,1,2]), charge_total=0, mult_total=1,
   autoff_qm2_method="gfnff")`.
2. `! qm/xtb dlpno-ccsd(t) def2-tzvpp def2-tzvpp/c tightscf`, the SAME block
   repeated, on step 1's geometry.

Against the Compound probe on the identical system:

| | Compound | OPI-driven |
|---|---|---|
| step 1 E / Eh | −81.500312 | −81.500312 |
| step 1 G(298.15) / Eh | −81.480860 | −81.480861 |
| step 2 E / Eh | −81.418060 | −81.418060 |
| composite G / Eh | −81.398608 | −81.398609 |

They agree to ~1e-6 Eh (0.003 kJ/mol) — two independent ORCA runs of the same
thing, not a difference of method.

**A finding that decides the hand-off:** for a QM/XTB job the property JSON
carries only the **QM-region** geometry — `geometries[-1].geometry` reports
natoms 3 while `calculation_info.numofatoms` says 6 — so OPI cannot hand you the
full system. The full geometry has to come from the `.xyz` ORCA writes for the
job: per-directory, programmatically read, and not a hand-managed file.

**Recommendation: the OPI-driven route.** No new reader is needed — each step
has its own `.out`, so `load_job` works and `composite_free_energy_eh` stays the
single home of the arithmetic — the geometry hand-off is programmatic, and the
inputs are built from OPI's typed keywords. Compound's one real advantage
(geometry inheritance) is matched here, and its costs are worse for this repo:
`%qmmm` silently does not inherit (MM-only, E = 0.0, normal termination), its
steps have no `.out` so `load_job` cannot read them and `terminated_normally()`
returns False, and it puts the workflow in a second place nothing gates.

### 2026-09-20 — Compound prototype: what ORCA's own workflow tool does and does not inherit

Bench probe (`/tmp`, water dimer, QM region = one water, so the whole prototype
is ~100 s) of the shape we actually need — geometry tier `Opt Freq`, then
DLPNO-CCSD(T) as the QM1 single point, then the composite. What it establishes:

- **`%Compound` works in ORCA 6.1.1** and runs both steps from one input. Each
  step writes the input ORCA actually ran, `job_Compound_<n>.inp`, and — with
  `%base "stepN"` and `%output jsonpropfile true end` **inside the step** —
  `stepN.property.json` / `stepN.property.txt`.
- **The GEOMETRY inherits across steps.** Step 2 ran on step 1's optimised
  geometry, to the printed precision, in both a plain and a QM/XTB compound.
  That is the whole provenance gain: there is no geometry file to copy and so
  nothing to poison, and the failure that cost us `energy_qm1_dlpno/ts/job.xyz`
  today cannot happen in this shape.
- **The `%qmmm` BLOCK DOES NOT INHERIT, and failing to repeat it is SILENT.** A
  step with no `%qmmm` block printed `No QM2 atoms in system. Switching to QMMM
  only calculation.` then `No QM atoms in system. Switching to MM only
  calculation.`, ran as **pure MM**, and `get_final_energy()` read **0.0** —
  while still terminating NORMALLY. So every step must repeat its QM region, and
  a compound script needs a per-step assertion that it really is a QM/XTB step.
  This is a new silent-wrong-answer trap and the reason the QM region cannot be
  specified once.
- **The `%qmmm` block is passed through verbatim**, so it must be the working
  form exactly: `QMATOMS {…}` needs its own `END` before `Charge_total`, or
  ORCA dies with `Error in [QMMM] block - Scan QM atoms - Line 6
  (CHARGE_TOTAL)`.
- **OPI reads a Compound step that has no `.out`.** `Output(basename="stepN",
  working_dir=…)` + `parse(read_gbw_json=False)` works, and the thermochemistry
  is at `geometries[-1].thermochemistry_energies[0]` (`freeenergyg`), matching
  ORCA's printed `Final Gibbs free energy`. `orca_io.load_job` cannot load such a
  step — it demands `<basename>.out`.
- **`terminated_normally()` / `scf_converged()` return False without a `.out`**
  — the trap this module exists for — while the property JSON's
  `Calculation_Status.status` reads `NORMAL TERMINATION`. A step loader has to
  validate from the JSON, not from OPI's status calls.

Prototype number (water dimer, deliberately not a project target): step 1
`QM/XTB r2SCAN-3c Opt Freq` E = −81.500312, G(298.15) = −81.480860 Eh; step 2
`QM1 = DLPNO-CCSD(T)` E = −81.418060 Eh; composite
`G = G_step1 + (E2 − E1)` = −81.398608 Eh. The combination is exactly the
formula already in `orca_io.composite_free_energy_eh` — nothing new is needed
for the arithmetic, only for the reading.

**To adopt (superseded the same day by the OPI-driven route above, which is what
we took):** an `orca_io` loader for Compound steps (JSON-only, status from
`Calculation_Status`), a per-step QM-region assertion, and then one compound
input per species (RC, TS, PC).

### 2026-09-20 — DLPNO-CCSD(T) can be the QM1 method directly in QM/XTB, and it moves the barrier

**It can, and it is the more complete route.** `! QM/XTB DLPNO-CCSD(T)
def2-TZVPP def2-TZVPP/C ddCOSMO(Water) TightSCF` runs unchanged (ORCA 6.1 manual
6.1.1: multiscale single points accept "all kinds of available electronic
structure methods as QM method"). ORCA reports `Embedding Scheme ...
electrostatic`, `Charge alteration scheme ... Charge shifting`, `Point charges in
QM calc. from MM atoms ... 120` and 6217 ddCOSMO surface charges — so the
high-level energy sees the cyclodextrin's own embedded QM2 charges, which the
post-hoc fragment correction cannot, because its fragment is isolated (CPCM
only) and ~9.7 Å of host is replaced by continuum. No `job.QMRegion.xyz`
extraction and no geometry-freshness trap; one job, **4 min 25 s** (the fragment
*pair* took ~6 min).

Pilot on C8 item 1, TS and both RCs (`energy_qm1_dlpno/`). The "correction" is
the QM1 level change, `E(QM1=DLPNO) − E(QM1=r2SCAN-3c)`, both QM/QM2 totals:

| species | post-hoc (isolated, CPCM) /kJ | direct (embedded) /kJ | diff |
|---|---|---|---|
| TS | 2068.15 | 2067.12 | −1.0 |
| RC `displace_plus` | 2045.72 | 2039.34 | −6.4 |
| RC `displace_minus` | 2045.68 | 2037.81 | −7.9 |

The host stabilises the **RC more than the TS**, so the barrier rises:

    ΔG‡ post-hoc 54.9 (plus) / 52.7 (minus)  →  direct 60.3 / 59.5 kJ/mol

a **+5.4 to +6.8 kJ/mol** move — larger than C8's own 4 kJ/mol target window.
The H131-rotamer spread grows from 0.04 to 1.5 kJ/mol, still inside it.

**Caveat, not resolved.** The direct route also swaps the continuum model
(ddCOSMO on the whole 140-atom system vs CPCM on the isolated fragment), and its
cavity is the *host's*. The shift therefore mixes "the host's point charges
polarise the QM1" with "a different continuum cavity". Rerunning the fragment in
ddCOSMO does not isolate them, because its cavity would still be the bare
fragment's. **Pick one route and use it for the whole register** — do not mix the
two — and say which when a C7–C10 number is quoted.

**Adopted the same day.** `energy_qm1_dlpno/` is now the energy tier for C8 (and
the pattern for C7/C10); `energy_dlpno-ccsdt/` is kept as the fragment-route
comparison for the isolated-species tasks. Nothing new was needed in the reader:
the two routes are the same two calls, so the pair added for the direct route
(`qm1_level_correction_eh`, `direct_free_energy_eh`) was **withdrawn the same day
as a synonym** — `direct_free_energy_eh(low, high, T)` is exactly
`composite_free_energy_eh(low, T, high, low)` — and
`orca_io.composite_free_energy_eh` now documents both call patterns. The
duplicate was found by asking what the direct route actually computes, not by a
test. `test_orca_io.py` pins both routes and their disagreement. The gate's C8
fixture moved with the RC: it now points at `rc/displace_minus` (total
`-873.7566473144415`), because the re-run replaced the old 12-atom `rc/job.out`
the constants used to describe.

### 2026-09-20 — the re-run RC converged to a saddle; displaced re-minimisations running

The 20-atom RC minimisation converged (61 cycles, gradients inside tolerance,
terminated normally in 1 h 17 m) but the final frequency check makes it a
**saddle**: one imaginary mode past the cutoff at **−313.5 cm⁻¹** (plus the
usual −10.3 cm⁻¹ floppy noise), while the connectivity is still reactant-like
(intact ketone C129=O130 1.226 Å, H₂O₂ not bonded — O136···C129 2.94 Å, H135
still on O136). The BFGS/RFO minimisation slid onto a ridge and stopped there.

`orca_pltvib rc/job.hess 6` gives that mode; its two extremes (frames 5 and 14,
amplitude ±0.70) are the seeds for `rc/displace_plus/` and
`rc/displace_minus/`, each a fresh 20-atom `opt freq`. Whichever returns with
no imaginary mode past the cutoff is the RC minimum. **Both converged normally,
and both are minima by the cutoff — but they are two conformers of the same
pre-complex, not one.** `displace_minus` carries only **−16.3 cm⁻¹** (the floppy
libration `orca_io` documents as noise) with G(298.15) = −872.775959 Eh;
`displace_plus` carries **−49.8 cm⁻¹**, just inside the 50 cm⁻¹ cutoff and on a
heavy-atom/water libration, with G(298.15) = −872.776846 Eh — **2.3 kJ/mol
lower**. Both are reactant-like (ketone 1.225 Å, H₂O₂ intact, O136···C129
2.92–3.05 Å); they differ almost entirely in H131 (O130···H131 2.14 Å in
`displace_plus`, pre-organised towards the TS, against 2.60 Å in
`displace_minus`). The −313 cm⁻¹ mode is gone from both. Against the verified TS
(G(298.15) = −872.764465 Eh) that gives **ΔG‡ = 30.2 kJ/mol** with
`displace_minus` and **32.5 kJ/mol** with `displace_plus`; the 2.3 kJ/mol spread
is the H131-rotamer systematic and sits inside C8's own 4 kJ/mol window.

**The composite barrier, computed the same day** (`energy_dlpno-ccsdt/`: the
DLPNO-CCSD(T)/def2-TZVPP (+def2-TZVPP/C) and r2SCAN-3c single points on the
22-atom `job.QMRegion.xyz` fragments, both CPCM(water), assembled with
`orca_io.composite_free_energy_eh` at 298.15 K; inputs modelled on
`K+peroxide+BnOH_bridged_UNIDENTIFIED/energy_dlpno-ccsdt/`):

| species | G(QM/XTB)/Eh | DLPNO correction/kJ/mol | G(composite)/Eh |
|---|---|---|---|
| TS | −872.764465 | +2068.15 | −871.976748 |
| RC `displace_plus` | −872.776846 | +2045.72 | −871.997675 |
| RC `displace_minus` | −872.775959 | +2045.68 | −871.996803 |

**ΔG‡ = 54.9 kJ/mol (13.1 kcal/mol)** with `displace_plus`, **52.7 kJ/mol
(12.6 kcal/mol)** with `displace_minus`. The DLPNO correction is 22 kJ/mol
larger on the TS than on the RC, so the final tier raises the barrier well above
the 30–33 kJ/mol the geometry tier alone gave — the r2SCAN-3c//GFN-FF barrier is
not the number to quote. The H131-rotamer spread survives essentially unchanged
(−2.3 kJ/mol at both tiers), so it stays inside C8's 4 kJ/mol window and does
not choose between the two RCs.

Two traps recorded while setting this up:

- **`geometry_r2scan3c-xtb/ts/job.xyz` is NOT the verified TS.** A crashed
  OptTS attempt overwrote it (O136–C129 2.596 vs the verified 1.820 Å, RMSD
  0.447 Å), and `ts/job.QMRegion.xyz` belongs to that bad geometry — verified
  atom-by-atom (its O136–C129 is 2.596 Å). The verified geometry is
  `ts/verify_exact_hessian/job.xyz`; its fragment was re-emitted by the
  one-cycle opt in `ts/qmregion_extract/` (heaviest-atom mismatch 0.004 Å,
  O136–C129 1.820 Å). A plain ORCA single point does **not** write
  `job.QMRegion.xyz` — only an optimisation does.
- The TS fragment's DLPNO energy (−647.148010 Eh) is 0.023 Eh above the RC's
  (−647.171059 Eh), while at r2SCAN-3c the gap is only 0.0145 Eh — the whole
  reason the composite barrier is larger than the geometry-tier one.

Still missing for C8: **the product PC**, for which the TS must be followed
along O136–C129 (IRC or relaxed scan) and re-minimised at 20 atoms (`opt freq`),
and then the separated K + H₂O₂ reference if a binding free energy is wanted.

One trap recorded: `orca_pltvib`'s frames carry element + xyz **plus three
displacement columns**, and must be rewritten with element + xyz only. Fed
as-is, ORCA misreads them and aborts with `Zero distance between atoms … in
Cartesian2Internal` (a formatting artefact, not a real clash — the C–C pair it
names is 1.54 Å apart).

### 2026-09-19 — reactant re-run at the enlarged QM region; the C8 pipeline

With the TS verified (one imaginary mode at −576.3 cm⁻¹), the reactant had to
be brought to the same level: `rc/job.inp` still used the **12-atom** QM region
(the carbonyl capped as H2C=O), while the TS uses **20 atoms**. The RC was
re-run at 20 atoms, `AutoFF_QM2_Method GFNFF`, ddCOSMO(Water), `opt freq`,
seeded from the converged 12-atom geometry (`rc/job.xyz`). The 12-atom input
and its only surviving output (`job.property.json`) are kept in
`rc/attempt1_qm12/`.

Two things block a barrier from the existing files, both recorded here so they
are not rediscovered:

- **The old `pc/job.xyz` is not the relay product** — H131 is still on water1
  (0.97 Å to O132) and H133 is mid-transfer, so relaxing it gives a different
  species. The product must be obtained by following the TS (IRC, or a relaxed
  scan along O136–C129) and then `opt freq` at 20 atoms.
- **`rc/job.out`, and every other `job.out`, were deleted** in a home-directory
  cleanup (2026-09-18/19). Only `job.property.json` survives for the old RC.

Remaining pipeline for C8 ΔG‡/ΔG_rxn: (a) this RC re-run; (b) the product as
above; (c) the composite energies — DLPNO-CCSD(T)/def2-TZVPP (+def2-TZVPP/C)
and r2SCAN-3c single points on each species' QM-region fragment, assembled with
`orca_io.composite_free_energy_eh`. The repo rule is explicit: the generic
`reaction-kinetics` script refuses a composite, so the number comes from
`orca_io`, not that skill's `compute_rate.py`.

### 2026-09-19 — the OptTS will not converge on this host; the NumFreq is the TS

Two attempts to wrap the verified geometry in a canonical unconstrained `OptTS`
failed, both for reasons in the optimiser, not the geometry:

- `Calc_Hess true` (forward differences, 418 displacements) reported **4**
  negative eigenvalues where the central-difference `NumFreq` found **1**; the
  gradients were already converged but the P-RFO step chased the spurious
  scaffold modes for 350 cycles.
- `inhess Read` of the NumFreq `job.hess` crashed ORCA in LEANSCF. The `.hess`
  was written with empty active-region records (`$act_atom 0`), and the
  multiscale reader wants the Hessian's active-atom record to match the job's —
  it found "0 corresponding active atoms in Hessian file (should be 140)" and
  aborted parsing the property.
- Adding `ActiveAtoms {0:139} end` (the manual's own remedy) stopped the crash,
  but the warning persisted, so the Hessian was never mapped correctly and the
  optimiser wandered again (3 negative eigenvalues by cycle 22).

**Conclusion: the geometry is already a fully characterised TS** — stationary
gradients plus, from the exact central-difference `NumFreq` in
`verify_exact_hessian/`, exactly one imaginary mode at −576.3 cm⁻¹
(`orca_io.stationary_point` = "transition state") with thermochemistry. The
`OptTS` is redundant. Proceed to use this as the C8 TS: visualise the mode, add
the ddCOSMO `pc`, then run `reaction-kinetics`.

### 2026-09-19 — the TS is verified; the canonical unconstrained OptTS is running

The exact full-system `NumFreq` at attempt4's geometry
(`.../geometry_r2scan3c-xtb/ts/verify_exact_hessian/`) finished in 38 min:
`orca_io.stationary_point` = **transition state**, with exactly **one**
imaginary mode at **−576.3 cm⁻¹** (420 frequencies, full 140 atoms). So the
frozen-host + hybrid-Hessian workaround in attempt4 happened to converge on a
genuine first-order saddle of the *full* PES; the earlier multi-imaginary
counts (−209 to −1150 cm⁻¹) were at other points on the wandering trajectory,
not here.

The canonical **unconstrained** `OptTS` was then launched from that geometry
(`ts/job.inp`, correction 8): full `Calc_Hess true` / `Recalc_Hess 10`, no
`Hybrid_Hess`, no `Constraints`, keeping `TS_Mode {B 135 128}` and
`TS_Active_Atoms`. Gradients at this geometry are already inside tolerance, so
it should converge in ~one cycle and leave a canonical OptTS record with its
own exact Hessian. attempt4's frozen/hybrid input is kept in
`ts/attempt4_frozen_hybrid/`.

### 2026-09-18 — attempt3 stalled too; attempt4 uses a Hybrid Hessian and a frozen host

**attempt3** (`TS_Mode {B 135 128}` + `TS_Active_Atoms`, seeded from attempt2's
cycle-32 frame) also stalled: ~258 cycles, then killed. Gradients were already
inside tolerance by cycle 11 (rms 9.2e-6, max 5.3e-5) but the Hessian's
negative-eigenvalue count flipped 1–5 and the step criteria never held, so it
never declared convergence and wandered ~30 kJ/mol downhill. Forcing the C–O
mode did not fix it. Kept in `.../ts/attempt3_qm20_tsmode_stalled/`.

**attempt4** seeds from attempt3's cycle-20 frame (gradients converged,
neg_eig = 1) and changes the Hessian strategy: `Hybrid_Hess {128 … 139} end`
(numerical second derivatives for the relay atoms only — ORCA 6.0 §6.3.11.2,
for a delocalized concerted mode; the numerical Hessian drops from 418 to 37
displacements), `HESS_Modification EV_Reverse` (§7.26, for the >1-negative
case), and the non-QM host frozen via Cartesian `Constraints` +
`ReduceRedInts true` so the macrocycle's soft modes cannot drive the search.
Setup is clean (`Created 2 link atoms`, `ORCA HYBRID HESSIAN`); running. The
fallback if it still wanders is a relaxed scan / `NEB-TS` from RC→PC.

### 2026-09-18 — the 20-atom `OptTS` stalled; restarted with a forced TS mode

**attempt2 (launched 2026-09-17) ran ~11 h / ~240 cycles without converging and
was stopped.** It reached a stationary point at cycle 32 — rms gradient 2.5e-5
and max 1.25e-4, both inside tolerance — but the full-system numerical Hessian
never settled at one negative eigenvalue (the count flipped 1–4; 3 at cycle 32,
with large imaginary modes at −388, −519 and −1150 cm⁻¹, so it is not the
floppy-scaffold numerical noise the converged `rc` "minimum" shows), and the
search then wandered downhill ~26 kJ/mol for ~200 cycles. Kept in
`C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_r2scan3c-xtb/ts/
attempt2_qm20_fullhess_stalled/`.

**attempt3 restarted from that cycle-32 frame** (frame 32 of the archived
`job_trj.xyz`) with the same full Hessian (`Calc_Hess true`, `Recalc_Hess 10`)
and `AutoFF_QM2_Method GFNFF`, plus the ORCA 6.0 §6.3.11 remedies:
`TS_Mode {B 135 128}` to follow the forming O136–C129 bond rather than whichever
imaginary mode is lowest, and `TS_Active_Atoms {128 129 130 131 132 133 134 135
136}` (factor 1.5) so the relay centres' bonds are in the internal set. Setup is
clean (`Created 2 link atoms`, QM1 = 20 + 2); running. If it wanders again the
next step is a better guess from a relaxed scan / `NEB-TS` along RC→PC.

### 2026-09-17 — QM2 topology builder fixed: `AutoFF_QM2_Method GFNFF`

First attempt at the ddCOSMO `OptTS` for `K+H2O2_water-relay_to_KP`
(`C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_r2scan3c-xtb/ts/`, ORCA
6.1.1, `! QM/XTB r2SCAN-3c ddCOSMO(Water) OptTS Freq TightSCF`, `nprocs 8`,
`maxcore 4000`). The warm start is the previously converged ALPB-level TS
(`ts/alpb_uncorrected/job_converged_alpb.xyz`), inherited from
`orca_stuff/cat/cat+H2O2+H2O+H2O/ts/`.

**Outcome: aborted in multiscale setup, seconds in, before any QM1 energy or
gradient.** The default `AutoFF_QM2_Method XTB` perceived a QM–QM2 bond between
the carbonyl oxygen O130 and ring carbon C59 at 2.48 Å — a 1,5 non-bonded
contact across the folded dialkoxy ring, not a bond — added a third link atom
on top of the two real cuts (C129–C121, C129–C125), and the crude link-atom
pre-optimisation failed (`Optimization not converged after 3 cycles` →
`CANNOT OPEN FILE job_S_Link.ORCAFF.prms.tmp`). The geometry is a valid
stationary point; the fault is in ORCA's distance-based boundary perception,
which is built once at setup.

**Probe (scratch, same input, only `AutoFF_QM2_Method` varies):** default `XTB`
→ `Created 3 link atoms`, same fatal; `GFNFF` → `Created 2 link atoms`, clears
setup into the initial numerical Hessian. The GFNFF output confirms the change
is topology-only: `QM2 method … XTB2`, `Method for determining QM2 charges …
Hirshfeld`, `AutoFF method … GFN-FF`, three independent settings. Adopted as
the tier convention (see Conventions).

**Also 2026-09-17 — the QM region was enlarged before running.** The 12-atom
region left the carbonyl's two α-carbons in QM2, so C129 was capped by two link
H and the QM1 electrophile was H₂C=O, not the real 1,3-dialkoxy ketone (an
aldehyde-like carbonyl is much more electrophilic). The region now also carries
the α-carbons C121/C125, their CH₂ hydrogens, and the β-ether oxygens O65/O62 —
20 atoms, with the boundary now the two O–C(ester) cuts O65–C59 and O62–C56
(capping the β-O with link H models the substituent as –OH rather than –OR, a
residual approximation). A scratch probe confirmed 2 link atoms and a clean
setup with `GFNFF`; the production `OptTS` was then launched from the ALPB warm
start (ORCA 6.1.1). **Still open:** the geometry is an ALPB warm start, and the
run has not finished.

### 2026-09-02 — C9 specced, not started

`buffer/` §6. The buffer question had been asked as "which species" and never as
"acting on what": a general base acts on the catalyst, a buffer perhydrate acts
by carrying the peroxide to it. The two differ by a `[buf] × [H₂O₂]` interaction
and **0 of 88 runs vary both**, so the archive cannot separate them and neither
can the two-pH species test.

Two kinetic results point the same way without settling it. The pre-equilibrium
constraint that `induction/` §4b imposes on H₂O₂ — `d ln v/d ln[X] − d ln τ/d
ln[X] = 1`, exact for every K and every `[X]` — is **met on the buffer axis**
(+1.094 ± 0.150 at the window that passes its signal control) and **missed on
the peroxide axis** (2.6σ and 3.7σ short, on blocks that fail theirs). And the
buffer-free route rises with pH 2.45× faster than [HOO⁻] alone accounts for
(7.90× against 3.23× over 0.50 pH units at matched substrate), which is the
signature of a second base-dependent step. Both are eight-curve, two-run results
against a measured between-day step of 1.80× — pointers, not measurements.

C9 is the calculation, and its gate is that orthophosphate's perhydrate
equilibrium must come out unfavourable, which is known independently.

### 2026-09-02 — C7 and C8 specced, not started

Arose from `induction/ANALYSIS.md`, the same exercise as C5/C6 done at the other
end of the same curves. The archive again settled more than expected — the
induction needs the catalyst, ends on a clock rather than at a product
threshold, has an amplitude that is a fraction rather than a concentration, and
carries a barrier far too large to be physical (95 ± 16 kJ/mol on the one-phase
τ over four temperatures, 77 ± 12 on the window-free clock over six) — and
again stopped exactly where a barrier comparison begins. It also produced a
clean negative worth recording: **the archive holds no usable peroxide lever on
an induction period.** The two blocks that move `[H2O2]` inside a run both have
induction statistics that track their own signal-to-noise, because `[H2O2]` is
what sets the signal, so the one experiment that would test the obvious
candidate cannot be done on these data at all.

No calculation run yet. C7 has three gates and all three are measured; C8's
discriminating comparison is internal to the calculation and needs C7's resting
fraction to be worth running.

*Amended the same day.* C8 gained an external gate after the peroxide question
was pushed further: the adduct scheme fixes the difference of the two peroxide
orders at exactly 1 whatever K is, and both blocks that can test it fall short
by 2.6σ and 3.7σ. Inverting the shortfall as a trap gives ΔG° of perhydrate
formation between −6 and −10 kJ/mol from three routes, which is now what C8 is
scored against. `data/induction.py` carries the machinery
(`joint_peroxide_order`, `peroxide_saturation`, `trap_constant`) and
`test_induction` plants both schemes and checks the test can tell them apart.

### 2026-09-02 — C5 and C6 specced, not started

Arose from the deep dive in `product_fate/ANALYSIS.md` into why the
4OMe progress curves rise to a maximum rate and then fall at 0.3–1.1%
conversion. The archive settled more of it than expected: the fall tracks the
product and not the clock, the same chemistry without the catalyst does the
opposite, the rate is linear in the product rather than hyperbolic in it, and
the stationary level carries a substrate order measured before anyone looked
for it. What the archive cannot do is separate "the oxidant consumes the
product" from "the product scavenges the oxidant", or say why benzaldehyde does
neither — both of which are barrier comparisons, and both of which now have a
measured number to be checked against.

No calculation run yet. C6 is the one with a live validation gate
(9.9 kJ/mol, plus a sign for BnOH); C5's gate is a literature Hammett ρ that has
not been sourced yet.

### 2026-08-30 — C1 specced, not started

Arose from chasing `MECHANISM.md`'s observable question. The literature route
was tried first and partly succeeded: benzaldehyde's aqueous ε at 278–279 nm was
found, and the sheets turned out to declare the monitoring wavelength directly,
which together established that `e = 1.23` is benzaldehyde's own coefficient
rather than a differential one. It failed on benzoate at 285 nm — the two best
sources returned HTTP 403.

No calculation run yet. The Herzberg–Teller point above was the reason for not
simply running a vertical TD-DFT and reporting the number.
