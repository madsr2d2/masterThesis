# Computational task register

Pending, running and completed quantum-chemistry work for this thesis, with the
reasoning behind each task and a dated log of what was actually run. Companion
to `MECHANISM.md` (the chemistry) and `DATA_VERIFICATION.md` (the data).

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

## Environment

| | |
|---|---|
| ORCA | 6.1.0 — `~/orca_6_1_0/orca`, aliased as `orca` |
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

### 2026-08-30 — C1 specced, not started

Arose from chasing `MECHANISM.md`'s observable question. The literature route
was tried first and partly succeeded: benzaldehyde's aqueous ε at 278–279 nm was
found, and the sheets turned out to declare the monitoring wavelength directly,
which together established that `e = 1.23` is benzaldehyde's own coefficient
rather than a differential one. It failed on benzoate at 285 nm — the two best
sources returned HTTP 403.

No calculation run yet. The Herzberg–Teller point above was the reason for not
simply running a vertical TD-DFT and reporting the number.
