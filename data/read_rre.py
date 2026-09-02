"""
Reads the instrument's own .rre files, the layer below the .txt exports.

Every progress curve in the dataset arrives as a tab-separated .txt export, and
`data/data/dataNN.txt` is what the pipeline parses. Those exports were made by
hand, one per run, and the archive shows that the hand stopped before the end:
`data/Mads` holds 43 instrument runs with no corresponding .txt, and two of
them -- exps 33 and 133 -- have no .xls sheet either, so nothing in the
repository knew they existed.

The .rre is a VisPro binary. Its layout, worked out by hand and confirmed
against the exports:

    "Sample00N"                  block header, one per cuvette
    ... float64 t0, dt, %T ...   little-endian doubles, %T in percent
    "BestFit1"                   end of the numeric block

Absorbance is the export's own transformation of that, A = -log10(%T/100), and
reading it back reproduces `data34.txt` to the last of its three decimals in all
four samples -- which is what licenses using this module on a run that has no
export to check against.

Since 2026-08-31 this module also supplies the READINGS for runs that were
exported, through `read_all`. The .txt export is rounded to three decimals of
absorbance, and on 67 of the 119 in-scope curves that rounding erases the
point-to-point scatter entirely -- their measured noise is exactly zero, which
is why every noise estimate in this package is floored at QUANTISATION_SIGMA
(2.89e-4 AU). The binary stores %T at about 2.1e-4 %, or 9.3e-7 AU, and from it
the same curves have a real noise of 1.15e-4 to 7.5e-4 AU, median 1.8e-4: the
floor had been overstating the instrument's noise by about 1.6x.
`fit_dataset.read_all_curves` therefore prefers the .rre wherever one exists
AND agrees with the export, and records which source each curve came from in
`Curve.source`. Since 2026-09-02 that is **all 402 fittable curves**; the
fallback path is kept because a disagreement must never be silently
substituted, not because anything currently takes it.

Two regex defects got it there, and they are the same defect twice. Until
2026-09-01, 125 curves were on the export because `experiment_number` matched
only `rate<n>.rre` and never the 32 `mads_t<n>.rre` files covering exps 2-32.
Until 2026-09-02, 28 more were, because the instrument wrote sample 3's label
as `sample003` -- lowercase -- in 31 files across exps 1-32 and the block
pattern was case-sensitive. Neither announced itself: the export exists, so a
fallback produces a curve rather than an error, and the second one had even
been written into `test_read_rre` as an expected count.

What this module still does NOT do is add CURVES to the dataset. A .rre
carries no conditions: no pH, no temperature, no concentrations, not even which
cuvette a sample sat in. For exp 33 those would all have to be inferred from a
neighbouring sheet, and an inferred condition record is exactly what
DATA_VERIFICATION.md exists to prevent. The curves are recoverable; whether they
are usable is a separate ruling that has not been made.

Usage:
    python data/read_rre.py                    # list the unexported runs
    python data/read_rre.py data/Mads/rate033.rre
    python data/read_rre.py --verify           # check the reader against a .txt
"""
import argparse
import glob
import os
import re

import numpy as np

ARCHIVE_DIR = "data/Mads"
EXPORT_DIR = "data/data"

# A block's numeric run starts within this many bytes of the "Sample00N" label;
# the exact offset varies with the string encoding, so it is searched for.
SEARCH_WINDOW = 80
# %T for a real reading. Wide enough for an over-range cuvette, tight enough
# that the first byte pattern that is not transmittance ends the run.
MIN_PERCENT, MAX_PERCENT = 0.5, 130.0
MIN_INTERVAL, MAX_INTERVAL = 1.0, 3600.0


def read_rre(path):
    """
    Reads every sample block of a .rre file.

    Args:
        path (str): Path to the .rre file.

    Returns:
        list: (sample_number, times, absorbance) tuples, in file order.
            times is in seconds and absorbance is -log10(%T/100).
    """
    with open(path, "rb") as handle:
        raw = handle.read()

    samples = []
    # [Ss]: the instrument wrote sample 3's label in LOWERCASE in 31 files of
    # the early campaign -- exps 1-32 -- and a case-sensitive pattern dropped
    # that cuvette from every one of them. Those 28 curves were being read from
    # the .txt export instead, at 1096x the noise floor, silently: the export
    # exists, so nothing reported a missing curve. Found 2026-09-02 from a plot,
    # because one panel per run looked coarser than its neighbours. Same class
    # of defect as the regex that hid 32 whole files until 2026-09-01.
    for match in re.finditer(rb"[Ss]ample00(\d)", raw):
        stop = raw.find(b"BestFit1", match.end())
        if stop < 0:
            continue
        found = None
        for offset in range(match.end(), match.end() + SEARCH_WINDOW):
            length = (stop - offset) - (stop - offset) % 8
            if length < 80:
                break
            values = np.frombuffer(raw[offset:offset + length], dtype="<f8")
            head = values[2:12]
            if (values[0] == 0.0 and MIN_INTERVAL <= values[1] <= MAX_INTERVAL
                    and np.all(np.isfinite(head))
                    and np.all((head > MIN_PERCENT) & (head < MAX_PERCENT))):
                found = values
                break
        if found is None:
            continue

        start, interval, rest = found[0], found[1], found[2:]
        good = np.isfinite(rest) & (rest > MIN_PERCENT) & (rest < MAX_PERCENT)
        count = len(rest) if good.all() else int(np.argmin(good))
        percent = rest[:count]
        times = start + interval * np.arange(count)
        samples.append((int(match.group(1)), times, -np.log10(percent / 100.0)))
    return samples


def read_export(path):
    """
    Reads a dataNN.txt export into the same shape read_rre returns.

    Kept local rather than imported from kinetics_io so that --verify compares
    the binary against the file on disk, not against a shared parser that could
    be wrong in both directions at once.
    """
    with open(path) as handle:
        lines = [line.rstrip("\n").split("\t") for line in handle]
    header = next(i for i, row in enumerate(lines) if row and row[0] == "ss.sss")
    columns = len(lines[header]) // 2
    samples = []
    for column in range(columns):
        times, values = [], []
        for row in lines[header + 1:]:
            if len(row) <= 2 * column + 1 or not row[2 * column].strip():
                continue
            times.append(float(row[2 * column]))
            values.append(float(row[2 * column + 1]))
        samples.append((column + 1, np.array(times), np.array(values)))
    return samples


def unexported_runs(archive=ARCHIVE_DIR, exports=EXPORT_DIR):
    """
    Instrument runs that were never exported to .txt.

    Returns:
        list: (number, rre_path, has_sheet) sorted by experiment number.
    """
    exported = {int(re.search(r"data(\d+)\.txt$", p).group(1))
                for p in glob.glob(os.path.join(exports, "data*.txt"))}
    sheets = set()
    for path in glob.glob(os.path.join(archive, "**", "*.xls"), recursive=True):
        match = re.search(r"mads_t(\d+)", os.path.basename(path))
        if match:
            sheets.add(int(match.group(1)))

    missing = []
    for path in sorted(glob.glob(os.path.join(archive, "*.rre"))):
        number = int(re.search(r"(\d+)\.rre$", path).group(1))
        if number not in exported:
            missing.append((number, path, number in sheets))
    return missing


def verify(number=34, archive=ARCHIVE_DIR, exports=EXPORT_DIR):
    """
    Checks the reader against an experiment that has both a .rre and a .txt.

    Returns:
        tuple: (matched_samples, total_samples, worst_absolute_difference)
    """
    candidates = [p for p in glob.glob(os.path.join(archive, "*.rre"))
                  if int(re.search(r"(\d+)\.rre$", os.path.basename(p)).group(1)) == number]
    binary = read_rre(candidates[0])
    export = read_export(os.path.join(exports, f"data{number}.txt"))

    matched, worst = 0, 0.0
    for (_, _, from_rre), (_, _, from_txt) in zip(binary, export):
        length = min(len(from_rre), len(from_txt))
        gap = float(np.max(np.abs(np.round(from_rre[:length], 3) - from_txt[:length])))
        worst = max(worst, gap)
        matched += gap == 0.0
    return matched, len(export), worst


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("path", nargs="?", help="a .rre file to summarise")
    parser.add_argument("--verify", action="store_true",
                        help="check the reader against exp 34's .txt export")
    arguments = parser.parse_args()

    if arguments.verify:
        matched, total, worst = verify()
        print(f"exp 34: {matched}/{total} samples reproduce the .txt export "
              f"exactly at 3 dp (worst difference {worst:g} AU)")
        return 0 if matched == total else 1

    if arguments.path:
        for sample, times, absorbance in read_rre(arguments.path):
            print(f"sample {sample}: {len(times)} points, "
                  f"{times[0]:.0f}-{times[-1]:.0f} s at dt = {times[1] - times[0]:.0f} s, "
                  f"A {absorbance[0]:+.4f} -> {absorbance[-1]:+.4f}")
        return 0

    missing = unexported_runs()
    print(f"{len(missing)} instrument run(s) with no .txt export")
    orphans = [m for m in missing if not m[2]]
    print(f"{len(orphans)} of those have no .xls sheet either, so nothing in the "
          f"repository records them:")
    for number, path, _ in orphans:
        samples = read_rre(path)
        span = max((t[-1] for _, t, _ in samples), default=0)
        print(f"  exp {number:<4d} {os.path.basename(path):<16s} "
              f"{len(samples)} cuvettes, {span:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# The smallest step seen between distinct %T values in these files, and the
# absorbance it is worth near 100 %T. This is the .rre's resolution, and it
# stands in for ABSORBANCE_QUANTUM wherever a curve came from one: the .txt
# export rounds to 0.001 AU, which is about 1000x coarser.
TRANSMITTANCE_QUANTUM = 2.1e-4
RRE_QUANTUM = TRANSMITTANCE_QUANTUM / 100.0 / np.log(10.0)
RRE_SIGMA = RRE_QUANTUM / np.sqrt(12)


# The archive names instrument runs two ways: `rate<n>.rre` for most, and
# `mads_t<n>.rre` for the early ones. Until 2026-09-01 only the first was
# matched, so 32 files -- every run of exps 2-32 -- were skipped and those
# curves kept the .txt export's floor for no reason. See DATA_VERIFICATION.md.
RRE_NAME = re.compile(r"(?:rate|mads_t)0*(\d+)\.rre")


def experiment_number(filename):
    """
    The experiment a `rate<n>.rre` or `mads_t<n>.rre` belongs to, or None.

    THE NUMBER IN THE NAME IS NOT TAKEN ON TRUST. It only proposes a pairing;
    `fit_dataset._prefer_rre` then requires the binary to have the same number
    of points as the export and to track it to within the export's own
    rounding, per sample, before a single reading is substituted. That test is
    what identifies the `mads_t` files: `mads_t003.rre` carries seven blocks of
    227 points that match all seven of exp 3's exported cuvettes, which no
    filename convention could establish and no coincidence explains.

    Filenames get copied forward between runs, which is why the standing rule
    is sheet over filename -- the same reason this function's output is a
    proposal checked against the data rather than an answer.
    """
    match = RRE_NAME.fullmatch(os.path.basename(filename))
    return int(match.group(1)) if match else None


def read_all(archive=ARCHIVE_DIR):
    """
    {experiment: {sample: absorbance array}} for every rate<n>.rre in `archive`.

    Only the top level is scanned. `good data BnOH/` and `done/` hold copies of
    the same runs under the same names, and a copy quietly overwriting the
    original is the class of error DATA_VERIFICATION.md exists to prevent.

    A sample label appears twice in the binary, once for the data block and
    once in the trailer; the first with a decodable block wins.
    """
    found = {}
    if not os.path.isdir(archive):
        return found
    for filename in sorted(os.listdir(archive)):
        number = experiment_number(filename)
        if number is None:
            continue
        samples = {}
        for sample, _, absorbance in read_rre(os.path.join(archive, filename)):
            samples.setdefault(sample, absorbance)
        if samples:
            found[number] = samples
    return found


def covered(archive=ARCHIVE_DIR):
    """The experiments a readable .rre exists for."""
    return set(read_all(archive))
