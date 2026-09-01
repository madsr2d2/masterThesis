"""
Verifies every number quoted in ANALYSIS.md against the code that produces it.

The point is drift. A figure typed into a document is a copy, and copies in this
project have gone stale before -- the whole reason `curve_metrics` exists. So
each claim below is re-derived from the modules AND its rendered form is
required to appear verbatim in ANALYSIS.md. Change the code and this fails;
edit the prose inconsistently and this fails.

    python background_reaction/check_numbers.py
"""
import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))

import scope
from curve_metrics import ACCELERATION_SIGMA
from fit_dataset import source_floor
from summary_kinetics import fit_burst, fit_burst_bounded

DOCUMENT = os.path.join(HERE, "ANALYSIS.md")
FAILURES = []


def _normalise(text):
    """
    Fold the typographic characters prose uses onto the ones %-formatting emits.

    U+2212 MINUS SIGN reads better than a hyphen in a document and is what the
    tables are written with, but `f"{value:+.3f}"` produces U+002D. Comparing
    glyphs would make this checker fail on typography while a genuinely wrong
    number that happened to be typed with a hyphen went through.
    """
    return (text.replace("\u2212", "-")     # minus sign
                .replace("\u2013", "-")     # en dash
                .replace("\u00b1", "+/-"))  # plus-minus


def claim(label, rendered, present=True):
    """Require `rendered` to appear in the document."""
    text = _normalise(io.open(DOCUMENT, encoding="utf-8").read())
    ok = (_normalise(rendered) in text) == present
    print(f"  {'pass' if ok else 'FAIL'}  {label}: {rendered!r}")
    if not ok:
        FAILURES.append(f"{label} -- {rendered!r} not found in ANALYSIS.md")


def main():
    print("checking ANALYSIS.md against the code\n")

    print("the enzyme-free census")
    import pandas as pd
    data = pd.read_csv(os.path.join(os.path.dirname(HERE), "data",
                                    "experiment_data.csv"))
    free = data[data["[enz]"] == 0]
    claim("enzyme-free curves", f"{len(free)} enzyme-free curves over "
                                f"{free.experiment.nunique()} experiments")
    fixed = scope.frame(scope.BUFFER_FIXED)
    titration = scope.frame(scope.BUFFER_CONFOUNDED)
    claim("BUFFER_FIXED size", f"| 65, 67, 69, 70 | {len(fixed)} |")
    claim("BUFFER_CONFOUNDED size",
          f"| 3, 6 | {len(titration)} ({int(titration.live.sum())} live) |")

    print("\nthe separation")
    quad = scope.buffer_dependence(parameter="v0_quad")
    vmax = scope.buffer_dependence(parameter="vmax")
    claim("a, v0_quad", f"**{quad['order_s0_fixed']:+.3f} ± "
                        f"{quad['stderr_s0_fixed']:.3f}** (n = {quad['n_fixed']})")
    claim("a', v0_quad", f"**{quad['order_s0_titration']:+.3f} ± "
                         f"{quad['stderr_s0_titration']:.3f}** "
                         f"(n = {quad['n_titration']})")
    claim("a, vmax", f"{vmax['order_s0_fixed']:+.3f} ± "
                     f"{vmax['stderr_s0_fixed']:.3f} (n = {vmax['n_fixed']})")
    claim("coupling g", f"g = {quad['coupling']:.3f}")
    claim("d, v0_quad", f"**d = {quad['order_buf']:+.2f} ± "
                        f"{quad['stderr_buf']:.2f}** on `v0_quad`")

    print("\nthe result table")
    for estimator, marker in (("v0_quad", "**"), ("vmax", ""),
                              ("v0", ""), ("v0_whole", "")):
        result = scope.buffer_dependence(parameter=estimator)
        cross = scope.background_orders(scope.FREE_4OME_40C,
                                        terms=("s0", "buf"), parameter=estimator)
        claim(f"{estimator} row",
              f"| {marker}{result['order_buf']:+.2f} ± "
              f"{result['stderr_buf']:.2f}{marker} | "
              f"{marker}{cross['order_buf']:+.2f} ± "
              f"{cross['stderr_buf']:.2f}{marker} |")

    print("\nthe joint fit that cannot carry it")
    joint = scope.background_orders(scope.BUFFER_CONFOUNDED,
                                    terms=("s0", "buf"), within=True,
                                    parameter="v0_quad")
    # The live count is derived, not spelled: it went 8 -> 10 on 2026-09-01
    # when exps 3 and 6 moved onto the instrument's own readings and two more
    # of exp 3's cuvettes cleared the live-signal threshold.
    spelled = {8: "eight", 9: "nine", 10: "ten", 11: "eleven", 12: "twelve"}
    claim("joint VIF", f"VIF {joint['vif_s0']:.1f} and {joint['vif_buf']:.1f} "
                       f"on {spelled.get(joint['n'], joint['n'])} live curves; "
                       f"returns {joint['order_buf']:+.2f} ± "
                       f"{joint['stderr_buf']:.2f}")

    print("\nthe burst diagnostics")
    curves = scope.curves(scope.FREE_BNOH_ALL)
    plain = [fit_burst(c.times, c.absorbance) for c in curves]
    bounded_free = [fit_burst_bounded(c.times, c.absorbance, constrain=False,
                                      noise_floor=source_floor(c.source))
                    for c in curves]
    bounded = [fit_burst_bounded(c.times, c.absorbance,
                                 noise_floor=source_floor(c.source))
               for c in curves]
    negative = sum(1 for f in plain if np.isfinite(f.v0) and f.v0 < 0)
    claim("negative v0 count",
          f"**{negative} of {len(curves)} curves return a negative v₀**")
    always = [fit_burst_bounded(c.times, c.absorbance, constrain=True,
                                noise_floor=source_floor(c.source))
              for c in curves]
    claim("bounded count, lag branch shut everywhere",
          f"from **{sum(f.bounded for f in bounded_free)} to "
          f"{sum(f.bounded for f in always)} of {len(curves)}**")
    # `bounded` is the DEFAULT, constrain="auto": the branch is gated on each
    # curve's own acceleration rather than shut for the whole block.
    claim("bounded count, lag branch gated per curve",
          f"That bounds **{sum(f.bounded for f in bounded)} of {len(curves)}**")

    print("\nshape and scope")
    live_fixed = fixed[fixed.live]
    inscope = scope.frame()
    live_inscope = inscope[inscope.live]
    claim("background acceleration",
          f"| **{int((live_fixed.accel_z > ACCELERATION_SIGMA).sum())} of "
          f"{len(live_fixed)}** |")
    claim("in-scope acceleration",
          f"| **{int((live_inscope.accel_z > ACCELERATION_SIGMA).sum())} of "
          f"{len(live_inscope)} live** |")
    buffers = sorted(inscope.buf.round(3).unique())
    claim("in-scope buffer constant",
          f"`[buf]` = {buffers[0]:.3f} mM in all {len(inscope)} in-scope curves")
    everything = scope.frame(scope.FREE_BNOH_ALL)
    curved = int((everything.curvature_t.abs() > 3).sum())
    claim("curvature count", f"{curved}\nof {len(everything)} curves show "
                             f"curvature at |t| > 3")

    print("\nthe in-scope substrate order")
    orders = scope.orders("vmax", scope.PRIMARY_SCOPE, within=True)
    strong = scope.orders("vmax", scope.strong_runs(), within=True)
    claim("in-scope order", f"{orders['order_s0']:+.3f} ± "
                            f"{orders['stderr_s0']:.3f} over all "
                            f"{orders['n']} live in-scope curves")
    claim("strong-run order", f"({strong['order_s0']:+.3f} ± "
                              f"{strong['stderr_s0']:.3f} on the "
                              f"{len(scope.strong_runs())} strong runs)")

    print()
    if FAILURES:
        print(f"{len(FAILURES)} MISMATCH(ES):")
        for item in FAILURES:
            print(f"  - {item}")
        return 1
    print("ANALYSIS.md agrees with the code")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
