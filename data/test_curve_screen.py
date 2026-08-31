"""
Tests for curve_screen.py.

The property that matters here is not that the screen fires, but that it fires
on the right things. Three populations are pinned against each other:

    broken     experiment 25's dead cuvettes -- must be caught
    sound      experiment 26's replicates, and 25's surviving siblings -- must not
    slow       the bottom rungs of the titrations -- must lose their RATE and
               keep their SHAPE, because they are data, not failures

The last is the one that took the longest to get right and is the easiest to
regress: a screen that quietly drops slow curves deletes the K_M information.

    python data/test_curve_screen.py
"""
import sys

import numpy as np
import pandas as pd

from curve_screen import (LADDER_FACTOR, MINIMUM_PEERS, MINIMUM_POINTS,
                          MINIMUM_SHAPE_SIGMA, MINIMUM_WINDOW_QUANTA, PEER_Z,
                          POWER_SHIFT_LIMIT, add_eligibility, eligibility,
                          ladder_anomalies, peer_scores, peer_sets, screen,
                          unsafe_blocks, validate_power_cut)
from fit_dataset import build_curves
from summary_kinetics import buffer_concentrations, summarise, to_frame

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


def curve_row(**overrides):
    row = dict(experiment=1, sample=1, points=100, v0=5e-5, window_quanta=20.0,
               amplitude=0.05, noise=6e-4, conversion=0.01, s0=4.0, h2o2=82.5,
               hoo=2.5e-3, buf=75.0, e0=0.0, group="g")
    row.update(overrides)
    return row


def test_eligibility():
    print("per-use eligibility")
    good = eligibility(curve_row())
    check("a healthy curve is eligible for both", good.rate and good.shape)
    check("and says so", good.summary == "rate+shape", good.summary)

    slow = eligibility(curve_row(window_quanta=0.8, v0=4e-6))
    check("a curve below the power floor loses its RATE", not slow.rate)
    check("...but keeps its SHAPE -- it is a slow rung, not a failure", slow.shape)
    check("...and the reason names the resolution, not a fault",
          any("resolution" in text for _, text in slow.reasons),
          f"{slow.reasons}")

    faint = eligibility(curve_row(amplitude=3e-3, noise=6e-4))
    check("a curve that barely moved loses its SHAPE", not faint.shape)
    check("...and keeps its rate if the window still rises", faint.rate)

    short = eligibility(curve_row(points=MINIMUM_POINTS - 1))
    check("too few points loses both", not short.rate and not short.shape)
    check("ten-point curves are NOT too few -- experiment 26 is sound",
          eligibility(curve_row(points=10)).rate)

    check("conversion above 100% loses both",
          not eligibility(curve_row(conversion=1.5)).rate)
    check("absorbance past the detector range loses both",
          not eligibility(curve_row(amplitude=2.5)).shape)
    check("a curve with no fitted slope loses its rate",
          not eligibility(curve_row(v0=np.nan)).rate)
    check("the power floor is stated in quanta", MINIMUM_WINDOW_QUANTA > 0)
    check("the shape floor is stated in noise units", MINIMUM_SHAPE_SIGMA > 1)


def test_ladder():
    print("ladder dips and spikes")
    base = [curve_row(experiment=9, sample=i + 1, s0=s, v0=v)
            for i, (s, v) in enumerate(zip([1.0, 2.0, 4.0, 8.0],
                                           [1e-4, 1e-4, 1e-4, 1e-4]))]
    check("a flat ladder raises nothing",
          ladder_anomalies(pd.DataFrame(base)) == [])

    dipped = pd.DataFrame(base).copy()
    dipped.loc[1, "v0"] = 5e-6                      # 20x below both neighbours
    hits = ladder_anomalies(dipped)
    check("a rung far below both neighbours is a dip",
          len(hits) == 1 and hits[0]["kind"] == "dip", f"{hits}")

    spiked = pd.DataFrame(base).copy()
    spiked.loc[1, "v0"] = 2e-3                      # 20x above both
    hits = ladder_anomalies(spiked)
    check("a rung far above both neighbours is a spike",
          len(hits) == 1 and hits[0]["kind"] == "spike", f"{hits}")

    # the alternating case: two dead rungs make the SOUND one between them
    # look like a spike. This is why layer 2 may only nominate.
    alt = pd.DataFrame(base).copy()
    alt.loc[1, "v0"] = 4e-6
    alt.loc[3, "v0"] = 8e-6
    kinds = {h["sample"]: h["kind"] for h in ladder_anomalies(alt)}
    check("with two dead rungs the sound one between them reads as a spike",
          kinds.get(3) == "spike", f"{kinds}")
    check("...which is exactly why the ladder layer cannot convict",
          "spike" in kinds.values() and "dip" in kinds.values(), f"{kinds}")

    monotone = pd.DataFrame([curve_row(experiment=9, sample=i + 1, s0=s, v0=v)
                             for i, (s, v) in enumerate(zip([1., 2., 4., 8.],
                                                            [1e-5, 3e-5, 9e-5, 2.7e-4]))])
    check("a steeply rising but smooth ladder raises nothing",
          ladder_anomalies(monotone) == [],
          "a real trend must not look like a spike")
    check("the factor is above 1", LADDER_FACTOR > 1)


def test_peers():
    print("condition-matched peers")
    rows = [curve_row(experiment=e, sample=1, s0=4.0, v0=8e-5) for e in range(1, 6)]
    rows.append(curve_row(experiment=6, sample=1, s0=4.0, v0=4e-6))   # 20x low
    rows.append(curve_row(experiment=7, sample=1, s0=99.0, v0=8e-5))  # no peers
    frame = pd.DataFrame(rows)
    peers = peer_sets(frame)
    check("cuvettes at matched conditions find each other",
          len(peers[0]) == 5, f"got {len(peers[0])}")
    check("a curve at a unique condition finds none",
          len(peers[6]) == 0, f"got {len(peers[6])}")

    z, counts, centres = peer_scores(frame)
    check("the odd one out scores strongly negative", z[5] < -PEER_Z,
          f"z={z[5]:.2f}")
    check("its healthy peers do not", abs(z[0]) < 2.0, f"z={z[0]:.2f}")
    check("a curve without peers is NaN, not a pass",
          np.isnan(z[6]) and counts[6] == 0)
    check("the peer centre is reported for context",
          abs(centres[5] - 8e-5) / 8e-5 < 0.05, f"{centres[5]:.2e}")
    check("at least three peers are required", MINIMUM_PEERS >= 3)


def test_screen_combines():
    print("layer bookkeeping")
    rows = [curve_row(experiment=e, sample=1, s0=4.0, v0=8e-5) for e in range(1, 6)]
    rows.append(curve_row(experiment=6, sample=1, s0=4.0, v0=4e-6, window_quanta=0.5,
                          amplitude=3e-3))
    out = screen(pd.DataFrame(rows))
    bad = out[out.experiment == 6].iloc[0]
    check("a curve failing both convicting layers reaches layers=2",
          bad.layers == 2, f"layers={bad.layers}, l1={bad.layer1}, l3={bad.layer3}")
    check("healthy curves reach layers=0",
          (out[out.experiment != 6].layers == 0).all())
    check("the ladder column exists even when nothing fires", "ladder" in out.columns)

    # layer 2 must NOT count towards a conviction
    ladder = pd.DataFrame([curve_row(experiment=9, sample=i + 1, s0=s, v0=v)
                           for i, (s, v) in enumerate(zip([1., 2., 4., 8.],
                                                          [1e-4, 5e-6, 1e-4, 1e-4]))])
    out = screen(ladder)
    flagged = out[out.ladder != ""]
    check("the ladder layer fires", len(flagged) >= 1)
    check("...but does not on its own count as a conviction",
          (flagged.layers < 2).all(), f"{flagged[['sample','ladder','layers']].to_dict()}")


def test_on_real_data():
    print("the real dataset")
    curves, _ = build_curves()
    frame = screen(to_frame(summarise(curves, buffer_concentrations())))

    survivors = frame[frame.experiment == 25]
    check("experiment 25's surviving cuvettes are eligible for a rate",
          survivors.eligible_rate.all(),
          f"{survivors[['sample','eligible_rate','window_quanta']].to_dict('records')}")
    replicates = frame[frame.experiment == 26]
    check("experiment 26's ten-point replicates are eligible for both",
          replicates.eligible_rate.all() and replicates.eligible_shape.all())
    check("no surviving curve is convicted by two layers",
          (frame.layers < 2).all(),
          f"{frame[frame.layers >= 2][['experiment','sample']].to_dict('records')}")

    only_shape = frame[~frame.eligible_rate & frame.eligible_shape]
    check("a substantial set keeps its shape after losing its rate",
          len(only_shape) > 20, f"got {len(only_shape)}")
    check("those are the slow ones, not the broken ones",
          only_shape.v0.median() < frame[frame.eligible_rate].v0.median(),
          "expected the rate-ineligible curves to be slower")

    # the safety gate must still be catching the block we know it should
    unsafe, report = unsafe_blocks(frame)
    check("the power cut is judged unsafe where it truncates a ladder",
          "('BnOH', 25.0, 'Pyrophosphate')" in unsafe, f"unsafe: {unsafe}")
    check("...and safe on the best-designed block",
          "('4OMe-BnOH', 40.0, 'Phosphate')" not in unsafe)
    row = report[report.group == "('4OMe-BnOH', 40.0, 'Phosphate')"].iloc[0]
    check("on that block the cut moves the order by essentially nothing",
          abs(row.shift_sigma) < 0.1, f"shift {row.shift_sigma:+.2f} sigma")
    check("the safety limit is one standard error", POWER_SHIFT_LIMIT == 1.0)


if __name__ == "__main__":
    test_eligibility()
    test_ladder()
    test_peers()
    test_screen_combines()
    test_on_real_data()
    print(f"\n{len(FAILURES)} failure(s)" + (": " + ", ".join(FAILURES) if FAILURES else ""))
    sys.exit(1 if FAILURES else 0)
