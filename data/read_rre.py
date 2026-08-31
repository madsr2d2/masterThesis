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

What this module deliberately does NOT do is add anything to the dataset. A .rre
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
    for match in re.finditer(rb"Sample00(\d)", raw):
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
