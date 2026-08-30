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
    """Overwrites one column for one experiment (or one sample of it)."""
    mask = data.experiment == experiment
    if sample is not None:
        mask &= data["sample"] == sample
    data.loc[mask, column] = value


CASES = [
    # (name, corruption applied to a copy of the dataset, check that must fire)
    ("silently wrong pH",
     lambda d: _set(d, 2, "pH", 9.99), "metadata"),
    ("silently wrong buffer",
     lambda d: _set(d, 13, "buffer", "Carbonate"), "metadata"),
    ("substrate swapped without updating e",
     lambda d: _set(d, 2, "substrate", "BnOH"), "invariant"),
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
]


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

    print(f"\n{passed} passed, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    print("fault-injection tests\n")
    raise SystemExit(main())
