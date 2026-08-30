"""
Self-test for data/validate_dataset.py.

A validator that has never failed is indistinguishable from one that cannot
fail. This injects each class of corruption the validator claims to catch into
a copy of the dataset and asserts that it is caught -- so the clean run on the
real data means something.

Usage:
    python data/test_validator.py
"""
import os
import tempfile

import pandas as pd

from validate_dataset import validate, DATASET_PATH, MANIFEST_PATH


def _run(corrupt, tmpdir):
    """Applies `corrupt` to a copy of the dataset and returns the findings."""
    data = pd.read_csv(DATASET_PATH)
    corrupt(data)
    path = os.path.join(tmpdir, "corrupted.csv")
    data.to_csv(path, index=False)
    return validate(path, MANIFEST_PATH)


def _errors_for(findings, check):
    return [e for e in findings.errors if e["check"] == check]


def _set(data, experiment, column, value, sample=None):
    """
    Overwrites one column for one experiment (or one sample of it).

    Returns the dataframe so corruptions can be chained; the CASES above rely
    on the in-place edit and ignore the return value.
    """
    mask = data.experiment == experiment
    if sample is not None:
        mask &= data["sample"] == sample
    data.loc[mask, column] = value
    return data


# These fire on the clean dataset by design -- their enzyme block describes the
# planned with-enzyme cuvettes, which were never run -- so a deep case must not
# count them as having caught its injected fault.
ACCEPTED_ENZYME_DEVIATIONS = {32, 34, 35, 36, 37}

# These fire on the clean dataset by design -- their enzyme block describes the
# planned with-enzyme cuvettes, which were never run -- so a deep case must not
# count them as having caught its injected fault.
ACCEPTED_ENZYME_DEVIATIONS = {32, 34, 35, 36, 37}

CASES = [
    # (name, corruption applied to a copy of the dataset, check that must fire)
    ("silently wrong pH",
     lambda d: _set(d, 2, "pH", 9.99), "metadata"),
    ("silently wrong buffer",
     lambda d: _set(d, 13, "buffer", "Carbonate"), "metadata"),
    # e and the wavelength moved from "invariant" to "optics" on 2026-08-30,
    # when they were reframed as substrate conventions rather than per-sheet
    # measurements. Relabelling a substrate leaves its optics behind, which is
    # now an optics finding.
    ("substrate swapped without updating e",
     lambda d: _set(d, 2, "substrate", "BnOH"), "optics"),
    ("enzyme concentration lost",
     lambda d: _set(d, 2, "[enz]", 0.0), "invariant"),
    ("negative concentration",
     lambda d: _set(d, 2, "[sub]", -1.0), "invariant"),
    ("null concentration",
     lambda d: _set(d, 2, "[buf]", float("nan")), "invariant"),
    ("a sample dropped",
     lambda d: d.drop(d.index[(d.experiment == 2) & (d["sample"] == 4)], inplace=True),
     "structure"),
    ("pH varying within one experiment",
     lambda d: _set(d, 2, "pH", 3.0, sample=1), "structure"),
    ("an experiment appearing from nowhere",
     lambda d: _set(d, 2, "experiment", 9999), "coverage"),
    # The wavelength was declared in the manifest but compared against nothing
    # until 2026-08-30, so a drift in the SUBSTRATE_PROPERTIES hardcode would
    # have passed silently. Exp 14 is a 300 nm run with no open question.
    ("monitoring wavelength drifting from the sheet",
     lambda d: _set(d, 14, "abs", 285), "optics"),
    # Uniformity is what the dataset's comparability rests on: within a
    # substrate a wrong e is one global scale factor absorbed into the fitted
    # constants, but a per-experiment e silently rescales runs against each
    # other. Nothing asserted it until 2026-08-30.
    ("a per-experiment e sneaking in",
     lambda d: _set(d, 14, "e", 9.9), "optics"),
]


# The deep checks run outside validate(), so they are exercised directly. Each
# case corrupts the dataset and names the check that must fire.
DEEP_CASES = [
    # The case this check was built for: exps 79 and 80 compiled with [enz] = 0
    # because their "[Enz] mmol/l" cell holds 0.000001, and nothing noticed
    # until it was spotted by eye. Restoring that state must now raise.
    ("enzyme concentration lost to a broken cell",
     lambda d: _set(_set(d, 79, "[enz]", 0.0), 80, "[enz]", 0.0), "enzuse"),
    # A plausible-looking but wrong enzyme concentration -- the failure mode a
    # scanner cannot catch, since the value is neither missing nor absurd.
    ("enzyme concentration silently halved",
     lambda d: _set(d, 14, "[enz]", float(d.loc[d.experiment == 14, "[enz]"].iloc[0]) / 2),
     "enzuse"),
]


def _run_deep(corrupt, tmpdir):
    """Applies a corruption and returns verify_enzyme's findings."""
    from verify_enzyme import analyse as analyse_enzyme
    data = pd.read_csv(DATASET_PATH)
    corrupted = corrupt(data)
    path = os.path.join(tmpdir, "deep.csv")
    (corrupted if corrupted is not None else data).to_csv(path, index=False)
    findings, _ = analyse_enzyme(path)
    return findings


def main():
    print(f"baseline: validating the real dataset")
    baseline = validate()
    print(f"  {len(baseline.errors)} error(s), {len(baseline.warnings)} warning(s)")
    assert not baseline.errors, (
        "the real dataset must validate cleanly before the fault tests mean anything")

    passed = failed = 0
    with tempfile.TemporaryDirectory() as tmpdir:
        for name, corrupt, expected_check in CASES:
            findings = _run(corrupt, tmpdir)
            caught = _errors_for(findings, expected_check)
            if caught:
                passed += 1
                print(f"  PASS  {name:38s} -> {expected_check}: {caught[0]['message'][:56]}")
            else:
                failed += 1
                print(f"  FAIL  {name:38s} -> nothing raised under '{expected_check}'")

    with tempfile.TemporaryDirectory() as tmpdir:
        for name, corrupt, expected_check in DEEP_CASES:
            findings = _run_deep(corrupt, tmpdir)
            caught = [f for f in findings if f[1] == expected_check
                      and f[0] not in ACCEPTED_ENZYME_DEVIATIONS]
            if caught:
                passed += 1
                print(f"  PASS  {name:38s} -> {expected_check}: {caught[0][2][:56]}")
            else:
                failed += 1
                print(f"  FAIL  {name:38s} -> nothing raised under '{expected_check}'")

    print(f"\n{passed} passed, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    print("fault-injection tests\n")
    raise SystemExit(main())
