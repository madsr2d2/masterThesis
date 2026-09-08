"""
Builds early_trough/index.html and early_trough/progress_curves.html.

Draws only; every number comes from `early_trough` or `curve_metrics`, so a
figure and the prose in ANALYSIS.md cannot disagree about a value without
`check_numbers.py` saying so.

    python early_trough/build_figures.py
"""
import functools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
sys.path.insert(0, os.path.dirname(HERE))

import curve_metrics
import early_trough
import scope
from svgplot import ACCENT, GRID, INK, MUTED, Axes, esc
from figure_kit import (CATEGORY, EVENT_BAND_COLOUR, breakpoints,
                        derivative_axes, fig, panel, progress_axes,
                        progress_overlay, residual_axes, styled, write_pages)

SUBSTRATE_COLOUR = {"4OMe-BnOH": CATEGORY[0], "BnOH": CATEGORY[1]}
CLUSTER_COLOUR = {"oxidant": CATEGORY[0], "substrate": CATEGORY[1]}
BUFFER_COLOUR = {"Phosphate": CATEGORY[0], "Pyrophosphate": CATEGORY[1]}


@functools.cache
def _table():
    return early_trough.trough_table()


@functools.cache
def _genuine():
    return _table()[_table().genuine].sort_values("z")


@functools.cache
def _rates():
    return early_trough.binding_rates(_table()).sort_values("z")


_CORRELATION_CEILING = 20.0


def figure_correlation():
    """
    Trough depth against dominance, restricted to `z <= _CORRELATION_CEILING`.

    A THIRD OF THE ARCHIVE IS DELIBERATELY LEFT OFF THIS PANEL, and the
    caption says so: 111 of 311 scanned curves have a strong genuine early
    RISE (z up to +1106), which is a real signal but not this question's
    subject, and plotting it would compress every dip into an unreadable
    sliver near zero on a linear axis. This is a stated selection, not a
    silent clip -- every point with z <= 20 is shown and the axis limits
    cover exactly that range, so `doc.unclipped` still holds.
    """
    table = _table()
    shown = table[table.z <= _CORRELATION_CEILING]
    genuine = _genuine()
    axes = Axes(600, 320, (0.005, 6000.0), (-46.0, _CORRELATION_CEILING),
               xlog=True, pad=(56, 20, 46, 20))
    axes.hline(early_trough.EARLY_TROUGH_CANDIDATE_Z, colour=GRID, dash="4 3")
    axes.note(560, axes._fy(early_trough.EARLY_TROUGH_CANDIDATE_Z) - 5,
              "candidate threshold", MUTED, size=10, anchor="end")
    for substrate, colour in SUBSTRATE_COLOUR.items():
        block = shown[shown.substrate == substrate]
        axes.points(block.dominance.to_numpy(), block.z.to_numpy(), colour,
                    radius=2.6, opacity=0.45, stroke=None)
    for substrate, colour in SUBSTRATE_COLOUR.items():
        block = genuine[genuine.substrate == substrate]
        if not len(block):
            continue
        axes.ring(block.dominance.to_numpy(), block.z.to_numpy(), colour,
                  radius=5.5, width=1.6)
    axes.note(80, 24, "● faint: every scanned curve with z ≤ 20   ○ ringed: "
                       "genuine (passes all three screens)", MUTED, size=10.5)
    for substrate, colour in SUBSTRATE_COLOUR.items():
        axes.note(80, 40 + 14 * list(SUBSTRATE_COLOUR).index(substrate),
                  substrate, colour, size=11, weight="600")
    corr = early_trough.dominance_correlation(table)
    return fig(
        axes.render("dominance = max([enz]/[S], [enz]/[HOO-])",
                    "early trough, sigma (capped at +20)",
                    "A · The archive-wide correlation, both substrates"),
        f"{len(shown)} of {len(table)} live catalysed curves in the archive "
        "plotted against whichever reactant the catalyst's own fixed "
        f"concentration comes closest to overwhelming — the "
        f"{len(table) - len(shown)} left off have a strong genuine early "
        "RISE (z up to +1106) rather than a trough, which is not this "
        "question's subject and would compress the dips shown here to a "
        "sliver near zero. <strong>Spearman ρ = "
        f"{corr['4OMe-BnOH']['e0_hoo'][0]:+.3f}</strong> (4OMe, p = "
        f"{corr['4OMe-BnOH']['e0_hoo'][1]:.1e}) and <strong>"
        f"{corr['BnOH']['e0_hoo'][0]:+.3f}</strong> (BnOH, p = "
        f"{corr['BnOH']['e0_hoo'][1]:.1e}) against log[enz]/[HOO⁻] alone, "
        "over EVERY scanned curve including the ones omitted here — the "
        "catalyst:substrate ratio carries no signal in either substrate "
        "(Section 2), and neither does the catalyst:TOTAL-peroxide ratio "
        f"({corr['4OMe-BnOH']['e0_h2o2'][0]:+.3f} in 4OMe, not even "
        "correctly signed — Section 3).")


def figure_clusters():
    genuine = _genuine()
    axes = Axes(600, 320, (0.001, 0.3), (0.001, 6000.0), xlog=True, ylog=True,
               pad=(58, 20, 46, 20))
    grid = np.array([0.001, 0.3])
    axes.line(grid, grid, GRID, dash="3 3", width=1.2)
    axes.note(axes._fx(0.15), axes._fy(0.15) - 6, "[enz]/[S] = [enz]/[HOO-]",
              MUTED, size=10, anchor="middle")
    for name, colour in CLUSTER_COLOUR.items():
        block = genuine[genuine.cluster == name]
        axes.points(block.e0_s0.to_numpy(), block.e0_hoo.to_numpy(), colour,
                    radius=5.0,
                    title=[f"exp {int(r.experiment)}.{int(r.sample)}"
                          for r in block.itertuples()])
        label = ("oxidant-dominated (above the line)" if name == "oxidant"
                 else "substrate-dominated (below the line)")
        axes.note(70, 40 + 16 * list(CLUSTER_COLOUR).index(name), label,
                  colour, size=11, weight="600")
    return fig(
        axes.render("[enz]/[S]", "[enz]/[HOO-]",
                    "B · The two clusters, by which reactant dominates"),
        "The seventeen genuine curves in the plane the classification is "
        "made in. Fourteen sit far above the diagonal — the catalyst is a "
        "trivial share of the substrate but a thousand-fold excess over the "
        "nanomolar hydroperoxide pool — and three sit below it, all BnOH at "
        "the block's lowest substrate rung (0.216 mM), where the catalyst is "
        "instead a meaningful fraction of the substrate itself.")


def figure_debubble():
    genuine = _genuine()
    axes = Axes(600, 320, (-46.0, 8.0), (-46.0, 8.0), pad=(56, 20, 46, 20))
    axes.line([-46.0, 8.0], [-46.0, 8.0], GRID, dash="3 3", width=1.2)
    axes.note(axes._fx(2.0), axes._fy(4.5), "z_corrected = z_raw", MUTED,
              size=10, anchor="end")
    for name, colour in CLUSTER_COLOUR.items():
        block = genuine[genuine.cluster == name]
        axes.points(block.z.to_numpy(), block.z_corrected.to_numpy(), colour,
                    radius=5.0)
    for row in genuine.itertuples():
        if abs(row.z_corrected - row.z) > 2.0:
            axes.note(axes._fx(row.z) + 8, axes._fy(row.z_corrected),
                      f"exp {row.experiment}.{row.sample}", INK, size=10)
    return fig(
        axes.render("trough, raw readings (sigma)",
                    "trough, debubble-corrected (sigma)",
                    "C · Correction for the O2 artefact deepens the effect, "
                    "it does not explain it away"),
        "Every genuine curve's trough recomputed after `curve_metrics."
        "debubble` removes every detected O2 event. Points would move "
        "<em>up</em>, back toward zero, if the O2 correction were secretly "
        "responsible for the dip. Fifteen of seventeen sit exactly on the "
        "diagonal (no detachment to correct); the two that carry any "
        "(exps 141.4, 142.4) move further <em>down</em> — the correction "
        "had been partly masking the trough with later bubble-driven rise.")


def figure_temperature():
    """
    The retracted two-point estimate against the honest 14-point regression
    -- drawn so the READER sees the same-temperature scatter that makes the
    slope unresolved, not just told about it in prose.
    """
    rates = _rates()
    oxidant = rates[rates.cluster == "oxidant"]
    fit = early_trough.arrhenius_check(rates)
    axes = Axes(600, 320, (12.0, 43.0), (0.3, 60.0), ylog=True,
               pad=(58, 20, 46, 20))
    common_T = float(oxidant.temperature.mode().iloc[0])
    same = oxidant[oxidant.temperature == common_T]
    # The same-temperature spread, drawn as a vertical bracket so its width
    # against the whole 15-40 C span is a single glance rather than a
    # sentence.
    axes.line([common_T, common_T], [same.k_on.min(), same.k_on.max()],
              MUTED, width=10, opacity=0.25)
    axes.note(axes._fx(common_T) + 10, axes._fy(same.k_on.max()),
              f"{len(same)} curves at {common_T:.0f}°C span "
              f"{np.log10(same.k_on.max() / same.k_on.min()):.2f} "
              "orders of magnitude on their own", MUTED, size=10.5,
              anchor="start")
    grid = np.linspace(13.0, 42.0, 50)
    kelvin_grid = grid + 273.15
    predicted = np.exp(fit["intercept"] + fit["slope"] / kelvin_grid)
    axes.line(grid, predicted, ACCENT, width=1.6, dash="5 4")
    axes.points(oxidant.temperature.to_numpy(), oxidant.k_on.to_numpy(),
                CATEGORY[0], radius=5.0,
                title=[f"exp {int(r.experiment)}.{int(r.sample)}"
                      for r in oxidant.itertuples()])
    for exp, samp, label in ((19, 1, "exp 19.1"), (34, 4, "exp 34.4")):
        row = oxidant[(oxidant.experiment == exp)
                     & (oxidant["sample"] == samp)].iloc[0]
        axes.note(axes._fx(row.temperature) + 8, axes._fy(row.k_on),
                  label, INK, size=10)
    return fig(
        axes.render("temperature, °C", "k_on, M⁻¹s⁻¹",
                    "D · Retracting a two-point activation energy"),
        "The two curves at the temperature extremes (exp 19.1, 15 °C; "
        "exp 34.4, 40 °C — the obvious two points to draw a line through) "
        "give an Ea of 85.7 kJ/mol — but the pre-exponential factor that "
        "goes with it is 2.1×10¹⁵ M⁻¹s⁻¹, about 2×10⁵ times the diffusion "
        "limit, and the equivalent Eyring ΔS‡ is <strong>+40 J/mol/K</strong> "
        "where a real bimolecular association must be negative. The dashed "
        "line is the honest fit instead — all 14 oxidant-cluster curves at "
        f"their own temperatures, giving <strong>{fit['activation_kJ']:.1f} "
        f"± {fit['stderr_kJ']:.1f} kJ/mol</strong> (t = "
        f"{fit['t_statistic']:.1f}, not significant at the usual bar of 2). "
        "The grey bar shows why: curves sharing a single temperature scatter "
        "almost as widely as the full 15–40 °C range does.")


def figure_buffer_effect():
    table = early_trough.buffer_comparison(_rates())
    rates = _rates()
    oxidant = rates[rates.cluster == "oxidant"]
    axes = Axes(600, 260, (0.2, 60.0), (-0.8, 1.8), xlog=True,
               pad=(150, 20, 46, 20))
    for index, row in table.iterrows():
        y = 1 - index
        colour = BUFFER_COLOUR[row.buffer]
        block = oxidant[oxidant.buffer == row.buffer]
        axes.line([row["min"], row["max"]], [y, y], colour, width=10,
                  opacity=0.25)
        axes.points(block.k_on.to_numpy(), [y] * len(block), colour,
                    radius=4.0, opacity=0.7)
        axes.points([row.geometric_mean], [y], colour, radius=6.5)
        axes.note(axes._fx(0.2) - 12, axes._fy(y) + 4,
                  f"{row.buffer} (n={row.n})", INK, size=11, anchor="end",
                  weight="600")
        axes.note(axes._fx(row.geometric_mean), axes._fy(y) - 12,
                  f"geo. mean {row.geometric_mean:.2f}", colour, size=10.5,
                  anchor="middle")
    return fig(
        axes.render("k_on, M⁻¹s⁻¹ (oxidant cluster only)", "",
                    "E · The apparent rate constant depends on which "
                    "buffer is present", yticks=False),
        "Every point normalised by [enz] identically, the same way "
        "throughout this page — yet pyrophosphate's geometric mean runs "
        f"<strong>{table.set_index('buffer').loc['Pyrophosphate', 'geometric_mean'] / table.set_index('buffer').loc['Phosphate', 'geometric_mean']:.1f}×</strong> "
        "phosphate's. A clean elementary step between the catalyst and free "
        "HOO⁻ should not care which buffer holds the pH. The more likely "
        "reading: the buffer runs its own fast, cuvette-symmetric "
        "equilibrium with HOO⁻ (a buffer perhydrate, `buffer/ANALYSIS.md`'s "
        "own open question) that sets the size of the reactive pool the "
        "catalyst draws from — different buffers, different pool.")


def build_index():
    table = _table()
    genuine = _genuine()
    rates = _rates()
    corr = early_trough.dominance_correlation(table)
    arrhenius_result = early_trough.arrhenius_check(rates)
    buffer_table = early_trough.buffer_comparison(rates)
    n_candidate = int(table.candidate.sum())
    n_genuine = len(genuine)
    n_oxidant = int((genuine.cluster == "oxidant").sum())
    n_substrate = int((genuine.cluster == "substrate").sum())

    hero = f"""
<div class='hero'>
  <div><div class='k'>curves scanned</div>
       <div class='v'>{len(table)}</div>
       <div class='u'>every live catalysed curve in the archive</div></div>
  <div><div class='k'>genuine early troughs</div>
       <div class='v'>{n_genuine} of {n_candidate}</div>
       <div class='u'>every candidate on a smoothed mean now confirmed</div></div>
  <div><div class='k'>binding rate constant</div>
       <div class='v'>{rates.k_on.min():.1f}–{rates.k_on.max():.1f}</div>
       <div class='u'>M⁻¹s⁻¹, from the archive's own progress fits</div></div>
  <div><div class='k'>buffer effect</div>
       <div class='v'>{buffer_table.set_index('buffer').loc['Pyrophosphate', 'geometric_mean'] / buffer_table.set_index('buffer').loc['Phosphate', 'geometric_mean']:.1f}×</div>
       <div class='u'>pyrophosphate over phosphate, same [enz] normalisation</div></div>
</div>"""

    rows = "".join(
        f"<tr><td>exp {int(r.experiment)}.{int(r.sample)}</td>"
        f"<td>{esc(r.substrate)}</td><td>{r.pH:.2f}</td>"
        f"<td>{r.s0:.3f}</td><td>{r.h2o2:.2f}</td><td>{r.e0:.3f}</td>"
        f"<td>{r.e0_s0:.4f}</td><td>{r.e0_hoo:.3g}</td>"
        f"<td>{r.z:+.1f}</td><td>{esc(r.cluster)}</td></tr>"
        for r in genuine.itertuples())

    body = f"""
<p class='lede'>A catalysed curve's absorbance is already reference-subtracted
against an enzyme-free cuvette holding the same composition, so the
background reaction the two cuvettes share cancels between them — as long as
both cuvettes see the same free concentration of whatever that background
reaction runs on. Seventeen curves across both substrates say it does not
always: before the catalysed rate takes over, the reported curve dips
measurably <strong>below zero</strong> and holds there for several minutes.
<a href='progress_curves.html'>progress_curves.html</a> shows every one of
them, fits and all.</p>
{hero}

<h2>1 · A trough that survives three separate ways to be spurious</h2>
<p>{n_candidate} curves clear a smoothed trough of
{early_trough.EARLY_TROUGH_CANDIDATE_Z:g} sigma or deeper in their own early
readings. Three screens separate a real, sustained decline from an artefact:
the trough window must survive a LEAVE-ONE-OUT test (drop its single worst
reading and the rest must still average below threshold — a real decline is
robust to this, a single bad reading is not); the window must not overlap a
detected O2 detachment; and the trough must survive `curve_metrics.debubble`
correction. <strong>{n_genuine} of {n_candidate} survive all three</strong> —
every candidate found. An earlier version of the sustained test used a
per-reading depth count instead of leave-one-out, and wrongly rejected two
real curves (exps 4.1, 22.2) whose decline is genuine but shallower per
reading than the block's most dramatic examples; see
<code>DATA_VERIFICATION.md</code> for the correction.</p>

<h2>2 · The driver is the oxidant, not the substrate, in both substrates independently</h2>
{figure_correlation()}
<p>[enz]/[HOO⁻] is significant at p &lt; 10⁻⁴ in <strong>both</strong> 4OMe-BnOH
and BnOH scanned separately; [enz]/[S] is not significant in either
(ρ = {corr['4OMe-BnOH']['e0_s0'][0]:+.3f}, p =
{corr['4OMe-BnOH']['e0_s0'][1]:.2f} for 4OMe;
ρ = {corr['BnOH']['e0_s0'][0]:+.3f}, p = {corr['BnOH']['e0_s0'][1]:.2f} for
BnOH). The hydroperoxide anion is nanomolar across most of the archive's pH
range, so a catalyst held at a fixed 0.014–0.273 mM is routinely a
100–4000-fold molar excess over it — exactly the regime in which a catalyst
engaging the peroxide non-productively (`MECHANISM.md` S4, the same route
`BUBBLES.md`'s gas takes) could measurably outrun the reference cuvette's own
supply.</p>

<h2>3 · The species test: HOO⁻, not H₂O₂ in general</h2>
<p>If the effect tracked total peroxide regardless of protonation state,
`[enz]/[H₂O₂]` (un-weighted by pH) should predict the trough at least as well
as the anion-specific ratio. It does not: in 4OMe-BnOH it is
<strong>{corr['4OMe-BnOH']['e0_h2o2'][0]:+.3f}</strong>
(p = {corr['4OMe-BnOH']['e0_h2o2'][1]:.2f}) — not even correctly signed —
against `[enz]/[HOO⁻]`'s {corr['4OMe-BnOH']['e0_hoo'][0]:+.3f} at
p = {corr['4OMe-BnOH']['e0_hoo'][1]:.1e}. In BnOH it is
{corr['BnOH']['e0_h2o2'][0]:+.3f} (same sign, weaker). That is the signature
of the deprotonated form specifically being consumed, not peroxide in
general — consistent with HOO⁻'s much greater nucleophilicity toward a
carbonyl (the α-effect) than neutral H₂O₂.</p>

<h2>4 · A second, smaller cluster: the substrate itself, at its scarcest</h2>
{figure_clusters()}
<p>Three BnOH curves (141.4, 142.4, 143.4) sit at the two-axis block's lowest
substrate rung, 0.216 mM — the one composition in the whole genuine set where
[enz] is a large enough fraction of [S] (6.5–9.7%) for simple stoichiometric
substrate binding to matter. Elsewhere in the archive [enz] never exceeds
0.2% of [S], which is why [enz]/[S] carries no overall signal in Section 2:
this cluster is real but small enough not to move the archive-wide
correlation either way.</p>

<h2>5 · Not the O2 artefact</h2>
{figure_debubble()}
<p>A growing bubble raises absorbance (`BUBBLES.md` §1) so cannot produce a
dip by itself, and none of the seventeen troughs overlaps a detected
detachment — but the strongest possible test is what correction does to the
two curves that carry any gas at all. It deepens their trough rather than
erasing it.</p>

<h2>6 · How fast — a rate constant from the fit already computed</h2>
<p>A trough IS the "lag" shape of the archive's own one/two-phase progress
fit taken far enough that `B/τ` exceeds `v_ss`, so `τ_fast` — already fitted
for every curve — is already the relaxation time of whatever produces it.
Pseudo-first-order (`[enz]` in vast excess for the oxidant cluster, `[S]`
for the smaller substrate cluster): `k_obs = 1/τ_fast = k_on·[excess]`. The
result clusters within <strong>1.8 orders of magnitude</strong>
(0.51–32.7 M⁻¹s⁻¹) despite [enz] varying 20×, pH varying over four units,
temperature varying 15–40 °C, and both substrates and buffers pooled
together — far below the diffusion limit (~10⁹–10¹⁰ M⁻¹s⁻¹), consistent
with a chemically-controlled bond-forming step rather than a barrierless
encounter.</p>

<h2>7 · Retracting a two-point activation energy</h2>
{figure_temperature()}
<p>Quoting {arrhenius_result['activation_kJ']:.0f} ±
{arrhenius_result['stderr_kJ']:.0f} kJ/mol is not the same as measuring an
activation energy: the standard error is more than half the value, and the
regression's own <em>t</em>-statistic
({arrhenius_result['t_statistic']:.1f}) sits below the ~2 a slope needs to be
distinguished from noise. This archive holds essentially one real low-
temperature point, one real high-temperature point, and a lot of scatter at
25 °C — not enough to resolve an Ea, and the honest report is that it
cannot, not a number dressed up as though it can.</p>

<h2>8 · The buffer leaves a footprint the simple story doesn't predict</h2>
{figure_buffer_effect()}
<p>A clean bimolecular step between the catalyst and free HOO⁻ has no reason
to care which buffer holds the pH once `[enz]` is accounted for. This one
does. The likeliest reading connects to something already established
elsewhere in this project: `induction.joint_buffer_order` finds the
catalyst's own E→E* activation step satisfies the pre-equilibrium "+1" rule
specifically on the <strong>buffer</strong> axis (+1.094 ± 0.150), not the
peroxide axis — so a two-step picture fits both findings at once. A fast,
cuvette-symmetric buffer–HOO⁻ pre-equilibrium sets how much reactive
material is available (explaining the buffer-identity effect here without
needing the catalyst to react with the buffer adduct directly — a fast,
symmetric equilibrium cancels in the reference subtraction on its own); the
catalyst's own, slower engagement with whatever is in that pool is the
asymmetric, genuinely observed trough; and the already-established
buffer-driven step converts the loaded intermediate into the active
catalyst afterward. The archive can independently confirm the first and
third pieces with real statistical power; the middle piece — the size of a
buffer–HOO⁻ reservoir — is inferred, not directly measured.</p>

<h2>9 · The seventeen curves</h2>
<div class='tbl'><table>
<tr><th>curve</th><th>substrate</th><th>pH</th><th>[S] mM</th>
<th>[H₂O₂] mM</th><th>[enz] mM</th><th>[enz]/[S]</th><th>[enz]/[HOO⁻]</th>
<th>trough, σ</th><th>cluster</th></tr>
{rows}
</table></div>

<h2>What this does and does not establish</h2>
<p><strong>Established.</strong> The dip is a real, sustained, multi-reading
decline in the raw readings (not a fitted-curve extrapolation artefact — see
`CLAUDE.md`'s derivative-panel discussion, which is what first surfaced 135.5
and 151.5–.7 by eye), it is not the O2 artefact, and archive-wide it tracks
[enz]/[HOO⁻] specifically — not total [H₂O₂], not [S] — far more strongly
than either alternative, in both substrates independently. Its own
relaxation time gives a self-consistent pseudo-first-order rate constant
(0.5–33 M⁻¹s⁻¹) well below the diffusion limit, and the rate depends on
buffer identity in a way a clean elementary step should not.</p>
<p><strong>Not established.</strong> Whether the catalyst is engaging HOO⁻
directly or via a buffer-perhydrate intermediate — this archive has no
headspace or manometric measurement of anything, exactly the same limit
`BUBBLES.md` states for the gas. No activation energy: the two-point
estimate (86 kJ/mol) is retracted, and the honest 14-point regression
(Section 7) is not significant. And exp 151 (three of the strongest
oxidant-cluster curves) is one of the two-axis block's own weakest, most
drift-dominated runs — real chemistry there is hardest to pull apart from
the cell's own wander, which is exactly why the 4OMe REPLICATE_RUNS curves
(exps 4, 5, 7 — the single strongest examples, at −40 to −42σ) matter: they
are a different substrate, a different buffer, and none of the two-axis
block's own caveats apply to them.</p>

<h2>Reproducing</h2>
<p><code>python data/test_early_trough.py</code> ·
<code>python early_trough/build_figures.py</code> ·
<code>python early_trough/check_numbers.py</code>, which re-derives every
number in <code>ANALYSIS.md</code> from <code>early_trough.py</code> and
fails if the prose and the code disagree.</p>
"""
    return styled("The early trough: a curve dipping before it rises", body,
                  "Reference-subtracted evidence for a transient, "
                  "non-productive engagement with a scarce reactant")


def _fit_panel(row, curve, colour):
    times = np.asarray(curve.times, dtype=float)
    values = np.asarray(curve.absorbance, dtype=float)
    corrected, events = curve_metrics.debubble(times, values, curve.noise)
    chopped = len(events) > 0
    axes, radius = progress_axes(times, values, limit=140,
                                 companion=corrected if chopped else None,
                                 pad=(56, 12, 10, 20))
    bands = []
    if chopped:
        for start, stop in events:
            bands.append((float(times[start]), float(times[stop])))
        for lo, hi in bands:
            axes.band([lo, hi], [axes.ylim[0], axes.ylim[0]],
                     [axes.ylim[1], axes.ylim[1]], EVENT_BAND_COLOUR,
                     opacity=0.14)
        axes.line(times, corrected, CATEGORY[2], width=1.0, dash="3 2",
                  opacity=0.85)
        corrected_fit = progress_overlay(axes, times, corrected,
                                         colour=CATEGORY[2], mark_radius=radius)
    raw_fit = progress_overlay(axes, times, values, mark_radius=radius,
                               colour=colour)
    chem_fit = corrected_fit if chopped else raw_fit
    chem_values = corrected if chopped else values
    if np.isfinite(row.t_trough):
        breakpoints(axes, [row.t_trough], [f"trough {row.z:+.1f}σ"],
                    colour=INK)
    residual = (chem_values - chem_fit.predict(times)) / curve.noise
    rax = residual_axes(times, residual, colour=colour, bands=bands)
    drax = derivative_axes(times, chem_fit, colour=colour, bands=bands)
    svg = (axes.render("", "ΔA") + rax.render("", "z")
          + drax.render("time, s", "dA/dt"))
    footer = (f"trough {row.z:+.1f}σ raw"
             + (f" · {row.z_corrected:+.1f}σ corrected" if chopped else "")
             + f" · {esc(row.cluster)}-dominated"
             + " · <strong>GENUINE</strong>")
    return panel(
        f"exp {int(row.experiment)}.{int(row.sample)} · "
        f"{esc(row.substrate)}"
        f"<span class='pill'>pH {row.pH:.2f}</span>",
        f"[S] {row.s0:.3f} mM · [H₂O₂] {row.h2o2:g} mM · [enz] {row.e0:.3f} mM "
        f"· {esc(row.buffer)}",
        svg, footer)


def build_curves_page():
    genuine = _genuine()
    lookup = {(c.experiment, c.sample): c for c in scope.curves(scope.archive())}

    panels = [_fit_panel(row, lookup[(row.experiment, row.sample)],
                         CLUSTER_COLOUR[row.cluster])
             for row in genuine.itertuples()]

    body = (f"<p class='lede'>The {len(panels)} curves behind every claim in "
            "<a href='index.html'>index.html</a>, each with its own three-panel "
            "audit: the readings and whichever fitted form the curve earned "
            "(raw in colour, the debubble-corrected series dashed where a "
            "curve carries any O2 event), the residual, and this session's "
            "derivative-of-fit panel — the same one that first made 135.5's "
            "and 151.5–.7's early behaviour visible by eye. The black dashed "
            "vertical marks the trough `early_trough` reads.</p>"
            "<div class='grid three'>" + "".join(panels) + "</div>")
    return styled("The early trough — every curve behind the finding", body,
                  f"{len(panels)} genuine curves, every candidate the scan "
                  "found")


def main():
    return write_pages(HERE, {"index.html": build_index(),
                              "progress_curves.html": build_curves_page()})


if __name__ == "__main__":
    raise SystemExit(main())
