# What the buffer does, and to what

The buffer is in this archive three times over and has never been looked at as
one question. It is a **reagent** — the uncatalysed reaction is first order in
it, the catalysed turnover has an order in it that saturates, and the induction
shortens with it. It is a **confound** — substrate was added by volume and
displaced buffer, so `[buf]` falls 80 → 50 mM as `[S]` rises in every 4OMe run
and a substrate order measured there is an order in the pair. And it is the
**leading candidate** for the step `induction/` is trying to name: general
acid/base catalysis of E → E\*.

This folder holds what the archive can say about all three, and one thing it
cannot.

    data/buffer_role.py               the titrations, the species test, the maps
    python data/buffer_role.py        the whole argument, printed
    python data/test_buffer_role.py   the species test and the published numbers
    python buffer/build_figures.py
    python buffer/check_numbers.py

**Figures**: [`index.html`](index.html) is the presentation — five figures, A to
E, one per claim below.

Related: [`../background_reaction/`](../background_reaction/ANALYSIS.md) §5–6b
owns the uncatalysed buffer order and is not restated here;
[`../temperature_series/`](../temperature_series/ANALYSIS.md) §3 the catalysed
one; [`../induction/`](../induction/ANALYSIS.md) §6 the induction's; and
`COMPUTATIONAL.md` C7, whose model already carries a phosphate dianion because
of what is below, and C9, which is §6's calculation.

## 1. A titration at one pH cannot name the species

At a fixed pH the acid, the base and the total buffer are all **proportional to
one another** — the ratio is set by pH and pKa and nothing else. So an order in
`[buf]` measured at one pH is an order in *total* buffer and says nothing about
which species does the work. Every order in this archive, and every order in the
three documents above, is of that kind.

The escape is to titrate the buffer at **more than one pH** and watch the
catalytic coefficient move. Phosphate's second pKa here is **7.20**, so the base
fraction runs 0.387 at pH 7.00 and 0.681 at pH 7.53, and the three hypotheses
separate by a factor of three:

| | predicted ratio of the coefficient, pH 7.53 / pH 7.00 |
|---|---|
| general **base** catalysis | **1.76** |
| general **acid** catalysis | **0.52** |
| buffer a spectator (ionic strength, say) | 1.00 |

**The archive has that design.** It is §3, and it is not strong enough to use.

## 2. There are five buffer titrations, not two

Every run that steps `[buf]` with `[S]` held fixed inside it — 4OMe-BnOH,
40 °C, phosphate, 82.5 mM H₂O₂, **20 curves at three pH values**. Only exps 32
and 34 have been used anywhere in this project before.

| exp | pH | [S] mM | [buf] mM | order in [buf] | R² | coefficient d*v*/d[buf] | intercept |
|---|---|---|---|---|---|---|---|
| 34 | 7.00 | 8.251 | 3.125–25 | **+0.792 ± 0.185** | 0.902 | 2.50 × 10⁻⁶ | 1.50 × 10⁻⁵ |
| 32 | 7.00 | 8.251 | 50–200 | **+0.371 ± 0.018** | 0.995 | 1.86 × 10⁻⁷ | 3.29 × 10⁻⁵ |
| 35 | 7.50 | 8.251 | 50–200 | −0.301 ± 0.459 | 0.176 | −4.43 × 10⁻⁷ | 2.63 × 10⁻⁴ |
| 36 | 7.53 | 57.900 | 50–200 | +0.064 ± 0.082 | 0.232 | 3.19 × 10⁻⁷ | 4.91 × 10⁻⁴ |
| 37 | 7.53 | 12.228 | 50–200 | **+0.122 ± 0.039** | 0.829 | 3.23 × 10⁻⁷ | 2.73 × 10⁻⁴ |

Two things are visible before any fitting.

**The dependence saturates in `[buf]`.** At pH 7.00 the order is +0.792 ± 0.185
over 3.125–25 mM and +0.371 ± 0.018 over 50–200 — about first order at low
buffer, about half order above 50. This is `temperature_series` §3's result, and
it is what makes +0.400 rather than +0.803 the number that applies to a block
sitting at 50–80 mM.

**And it collapses by pH 7.5.** The three runs above the pKa give orders of
−0.301 ± 0.459, +0.064 ± 0.082 and +0.122 ± 0.039. Not zero — exp 37's is three
standard errors from it — but the largest of them is **a third** of the
+0.371 ± 0.018 measured over the same 50–200 mM range at pH 7.00.
Not because the buffer term shrank: its *coefficient* is 1.86 × 10⁻⁷ at pH 7.00
and 3.19–3.23 × 10⁻⁷ at pH 7.53. **The buffer-independent route grew instead**,
from an intercept of 3.29 × 10⁻⁵ to 2.73–4.91 × 10⁻⁴ — roughly eightfold. That
was attributed to [HOO⁻] when this section was written, and it is more than
[HOO⁻] can carry: §6 divides one by the other. The buffer is carrying the same
absolute rate and a much smaller share of a much bigger number.

## 3. Which species — and the archive cannot say

Fit `v/[S]^0.471 = a_run + b(pH group)·[buf]` over the four runs that share the
50–200 mM range: one free level per run, one slope on each side of the pKa. The
substrate normalisation is needed because the two pH groups do not sit at the
same `[S]` (8.251 against 12.228 and 57.900) and a free intercept absorbs that
from the level but not from the slope. It is an assumption — that the
buffer-catalysed route carries the same substrate order as the rest — and the
unnormalised fit is printed beside it so the reader can price it.

| fit | b(pH 7.00) | b(pH ≥ 7.5) | ratio |
|---|---|---|---|
| all four runs | 6.90 × 10⁻⁸ ± 1.4 × 10⁻⁷ | −5.79 × 10⁻⁹ ± 7.9 × 10⁻⁸ | −0.08 ± 1.16 |
| **without exp 35** | 6.90 × 10⁻⁸ ± 4.2 × 10⁻⁸ | 7.32 × 10⁻⁸ ± 3.0 × 10⁻⁸ | **+1.06 ± 0.77** |
| unnormalised | 1.86 × 10⁻⁷ ± 4.2 × 10⁻⁷ | 6.65 × 10⁻⁸ ± 2.4 × 10⁻⁷ | +0.36 ± 1.53 |

Exp 35 is dropped in the second row because its own titration has R² 0.176 and
is non-monotone — 2.46, 2.57, 1.03, 2.24 × 10⁻⁴ across 50 → 200 mM, one cuvette
at less than half its neighbours. That is a judgement about a run, made on the
run's own internal consistency and not on whether it helps.

**Nothing is excluded.** Against the best of the three fits:

| | predicted | distance | |
|---|---|---|---|
| general base | 1.76 | 0.9σ | survives |
| general acid | 0.52 | 0.7σ | survives |
| spectator | 1.00 | 0.1σ | survives |

The reason is in the inputs rather than the output: at pH 7.5 the
buffer-independent route is eight times larger, so the coefficient there is
44–127% uncertain and the ratio inherits it. **The design is right and the data
are too thin.** One titration at pH 6.5, where the base fraction is 0.17 and the
[HOO⁻] route is weakest, would do what four runs at 7.0–7.5 cannot.

## 4. The buffer as a confound, which is where it has done real damage

In every 4OMe run substrate volume displaced buffer volume:

| block | runs with a substrate ladder | median corr(log[S], log[buf]) | d log[buf]/d log[S] |
|---|---|---|---|
| 4OMe catalysed | 28 | **−0.96** | −0.264 |
| 4OMe enzyme-free | 12 | −0.97 | −0.193 |
| the temperature series | 6 | −0.96 | −0.325 |
| **BnOH in scope (135–151)** | 17 | **0.00** | — (17 runs with `[buf]` constant) |

Two orders have already had to be corrected for it, and the corrections went in
opposite directions in their consequences:

- **the rate's substrate order** — `temperature_series` §3, corrected with the
  measured buffer order of the rate, +0.45 → +0.58. It moved the level and not
  the trend.
- **the induction's substrate order** — `induction/` §3, corrected with the
  measured buffer order of the induction, −0.121 ± 0.148 → −0.235 ± 0.158. It
  moved that route from excluding product control to not excluding it.

Exps 135–151 are the one block where `[S]` moves alone, and it is the block the
in-scope fitting is scoped to — which was chosen for its peroxide ladder and
turns out to have been the right choice for this reason too.

## 5. Buffer identity cannot be separated from pH

The three buffers were each chosen for the pH they hold, which is what a buffer
is for and what makes an identity comparison a pH comparison wearing a label.

| substrate | channel | buffer | curves | pH |
|---|---|---|---|---|
| 4OMe | enzyme-free | phosphate | 49 | 6.97–7.50 |
| 4OMe | catalysed | phosphate | 92 | 5.64–8.95 |
| 4OMe | catalysed | boric | 40 | 8.46–10.34 |
| 4OMe | catalysed | pyrophosphate | 15 | 6.94–8.98 |
| BnOH | enzyme-free | phosphate | 22 | 6.71–8.01 |
| BnOH | enzyme-free | boric | 4 | 8.51 |
| BnOH | catalysed | phosphate | 19 | 7.50–8.01 |
| BnOH | catalysed | boric | 27 | 8.51–9.70 |
| BnOH | catalysed | pyrophosphate | 118 | 5.47–9.73 |

The widest pH range any two buffers share inside one substrate/channel cell is
**2.01 units** — phosphate against pyrophosphate on catalysed 4OMe, 6.94–8.95.
**And those two cells share no peroxide.** The phosphate block sits at exactly
82.5 mM; the pyrophosphate one at exactly 3.879 and 195.882 mM. A range overlap
is not a matched pair, and the check is on the values.

`background_reaction` §6b reaches the same wall from the other side: the boric
evidence for buffer catalysis is real and confounded, and boric is unusable for
rates at all (`scope.BORIC_RATE_UNUSABLE`).

## 6. Is the buffer a base, or is it carrying the peroxide?

§3 asks which *species* of buffer does the work and cannot answer. There is a
prior question it does not ask at all: whether the buffer acts on the catalyst
directly, or whether it acts by picking up the peroxide first.

    H₂O₂ + P     ⇌  P–OOH                 the buffer perhydrate
    P–OOH        ⇌  P–OO⁻ + H⁺
    K + P–OO⁻    ⇌  K(O⁻)–OO–P            the Criegee adduct at the ketone
    K(O⁻)–OO–P   →  KD + P–O⁻             the dioxirane, and the buffer back

This is the same shape as `MECHANISM.md`'s autocatalytic step with the buffer
where the peracid stands, and it has a reason to be taken seriously that is not
kinetic: **making a dioxirane out of plain H₂O₂ means expelling hydroxide from
the Criegee adduct**, and hydroxide is a very poor leaving group. That is
exactly why the same chemistry elsewhere is run on peroxymonosulfate or on a
peracid — a peroxide that arrives carrying a leaving group already. At 50–200
mM the buffer is the only species present in bulk that could supply one.

### The archive cannot separate it from general base catalysis

A buffer acting as a general base is a term in `[buf]`. A buffer acting through
its own perhydrate is a term in `[buf][H₂O₂]`, because the active species'
concentration is set by that product. **At a single pH the two schemes differ by
that interaction and by nothing else** — and this archive never varies it:

| | runs |
|---|---|
| runs in the archive | 88 |
| runs that step `[buf]` | 53 |
| runs that step `[H₂O₂]` | 20 |
| **runs that step both** | **0** |

All five titrations sit at exactly **82.5 mM** peroxide. Nor does §3's species
test help: nucleophilic attack by the phosphate dianion has the same pH
signature as general base catalysis *by* the phosphate dianion, so the
perhydrate scheme is a **fourth survivor** of §3 rather than a resolution of it.

### What the archive does say points the same way, twice

**One.** `induction/` §4b's constraint is not a fact about hydrogen peroxide. If
the catalyst is drawn into its active form by any species X held in excess,

    E + X ⇌ E*    1/τ = k_f[X] + k_r,    [E*]/E₀ = K[X]/(1 + K[X])

so `d ln v/d ln[X] − d ln τ/d ln[X] = 1` for **every** K and every `[X]`. It is
a prediction with nothing left to fit, and it transfers to the buffer axis
unchanged. `induction.joint_buffer_order` measures it on the eight curves of
exps 32 and 34:

| axis | joint order | from the required +1 | signal control |
|---|---|---|---|
| `[buf]`, 300 s window | **+0.989 ± 0.159** | 0.1σ | −0.145 ± 0.266 passes |
| `[buf]`, 450 s window | **+1.094 ± 0.150** | 0.6σ | −0.275 ± 0.268 passes |
| `[buf]`, 600 s window | +1.379 ± 0.180 | 2.1σ | −0.660 ± 0.230 **fails** |
| `[buf]`, 900 s window | +1.339 ± 0.207 | 1.6σ | −0.643 ± 0.225 **fails** |
| `[buf]`, 1200 s window | +1.232 ± 0.236 | 1.0σ | −0.565 ± 0.224 **fails** |
| `[H₂O₂]`, exps 127–131 | +0.502 ± 0.194 | 2.6σ | fails |
| `[H₂O₂]`, exps 135–151 | +0.304 ± 0.188 | 3.7σ | fails |

**The buffer axis meets the constraint the peroxide axis misses.** The windows
that fail their own signal control are the windows that overshoot, and they
overshoot in the direction the artefact predicts — more buffer, more signal, an
earlier landmark, an inflated gap — so the two windows that pass are the ones to
read, and they give +0.99 and +1.09.

**Two.** The buffer-*free* route rises with pH faster than one hydroperoxide can
explain. Exps 32 and 35 sit at the same substrate (8.251 mM) over the same
buffer range, so their levels divide with nothing to normalise:

| | pH 7.00 → 7.50 |
|---|---|
| buffer-free level | 2.65 × 10⁻⁵ → 2.09 × 10⁻⁴, **7.90×** |
| [HOO⁻] over the same ΔpH | **3.23×** (order +1.02 in [OH⁻], as it must be) |
| apparent order of the level in pH | **+1.79**, not +1.02 |

**2.45× more than one hydroperoxide accounts for** — the signature of a second
base-dependent step, which is what a scheme that has to deprotonate a
peroxo species *and* a Criegee adduct would show.

### What is wrong with reading either of those too hard

Both are eight-curve, two-run results, and the two runs in each pair are
different days. The archive's own measured between-day step is **1.80×**
(`induction.buffer_join_step`), which is of the same order as the 2.45× excess.
These are pointers. The one thing they are not is silent: the constraint in
route one is parameter-free, and it is met on one axis and missed on the other.

### The experiment, and it is a small change to the one §3 already asks for

Add the peroxide axis to the pH 6.5 titration: a **buffer × peroxide block at
one pH**, four rungs each. The interaction term is the whole question, and no
run in this archive has one. The same design would also test the sharpest
prediction the scheme makes about buffer identity — perhydrolysis needs an
electrophilic phosphorus, orthophosphate has none and **pyrophosphate has a
leaving group**, so the scheme predicts pyrophosphate ≫ phosphate as a buffer
catalyst at matched pH and matched peroxide. §5 is why that comparison cannot be
made from what exists.

`COMPUTATIONAL.md` **C9** is the calculation: ΔG° of `H₂O₂ + HPO₄²⁻ ⇌
HOO–PO₃²⁻ + H₂O`. Peroxomonophosphate hydrolyses readily, so that equilibrium
is expected to lie well to the left; the scheme survives a small equilibrium
fraction, since catalysis through a reactive minor species is ordinary, but it
does not survive one that cannot be reached at 200 mM buffer.

## 7. What this settles, and what it does not

**Settled.**

- Any order in `[buf]` in this project is an order in **total buffer**. The
  species is not identified anywhere and cannot be, from one pH.
- The catalysed turnover's buffer dependence **saturates**: about first order
  below 25 mM, about half order above 50, at pH 7.00.
- It **disappears above the pKa**, not because the buffer term shrinks but
  because the buffer-free route grows eightfold past it — which is **2.45×
  more than [HOO⁻] alone accounts for**, so that route has a second
  base-dependent step in it.
- `[S]` and `[buf]` are collinear at **−0.96** in every 4OMe run and at
  **0.00** in every in-scope run, which is why the two blocks disagree about the
  substrate and why two published orders needed correcting.
- Buffer **identity** is confounded with pH everywhere, and the one 2-unit pH
  overlap has no shared peroxide.

**Not settled.**

- **Acid, base or spectator.** §3: all three survive at 2σ. This is the folder's
  main negative and the archive cannot fix it.
- **Whether the buffer is what carries E → E\*.** `induction/` §6 has the
  direction — more buffer, shorter induction, −0.433 ± 0.201 — and it is eight
  curves at 40 °C, where the induction is nearly over.
- **Whether the buffer acts on the catalyst or on the peroxide.** §6: the two
  differ by a `[buf] × [H₂O₂]` interaction, and **0 of 88 runs** vary both.
  What the archive does have is a parameter-free constraint that the buffer
  axis meets (+1.094 ± 0.150 against +1) and the peroxide axis misses.

**The two measurements that would settle both** are one run and a calculation.
A buffer titration at **pH 6.5 and 25 °C on the catalysed 4OMe system**: the low
pH makes the buffer term the largest share of the rate it can be, the low
temperature makes the induction thousands of seconds long, and one design would
give the species *and* the induction's buffer order at once. And
`COMPUTATIONAL.md` **C7**, whose model already carries a phosphate dianion
alongside the ketone hydrate for this reason — if the computed dehydration
barrier falls when the dianion is present, general base catalysis stops being a
hypothesis with three survivors.
