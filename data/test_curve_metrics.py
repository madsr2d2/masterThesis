"""
Tests for curve_metrics.py, and the guard that keeps it the only definition.

test_no_duplicate_definitions is the one that earns its place. Six functions
were defined in two modules each and four had diverged; the lag statistic's two
copies disagreed on 96 of 402 curves and would have put 21% in a thesis where
the evidence says 34%. That class of bug is invisible in review and expensive
in print, so it is now a test failure.

    python data/test_curve_metrics.py
"""
import ast
import glob
import os
import sys

import numpy as np

from curve_metrics import (ACCELERATION_SIGMA, BUBBLE_DROP_SIGMA,
                           DETACHMENT_SNR_FLOOR,
                           EXCURSION_RECOVERY_CEILING,
                           EXCURSION_RECOVERY_DEPTH,
                           INITIAL_WINDOW, LAG_THRESHOLD,
                           QUANTISATION_SIGMA, acceleration, curve_noise,
                           early_trough,
                           initial_rate, line_fit, line_slope, peak_position,
                           OUTLIER_SIGMA, apply_gains, bubble_drops,
                           bubble_gains, bubble_load,
                           bubble_profile, bubble_rate, bubble_shortfall,
                           local_outlier_z, OUTLIER_SIGMA,
                           debubble, detachments, isolated_outliers,
                           monotone_bound, tail_excess,
                           terminal_gas, _is_excursion,
                           local_outlier_z, model_residual, quadratic_rate,
                           segmented_fit, segment_breaks,
                           segment_selection, _segment_errors,
                           SEGMENT_RATIO_STEEP,
                           whole_slope, window_size)
from fit_dataset import source_floor
from read_rre import RRE_SIGMA
import scope

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  pass  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILURES.append(name)


# Names every module is allowed to define for itself: per-module entry points
# and test harness helpers, not shared measurements. `build_index` and
# `build_curves_page` are the analysis folders' two page entry points -- one of
# each per folder, called only by that folder's `main`, so they can no more
# drift into each other than `main` can.
PERMITTED_DUPLICATES = {
    "main", "build", "analyse", "report", "check", "close",
    "test_regressions", "build_index", "build_curves_page",
}

# Constants every module is allowed to define for itself. Each is a module's
# own LOCATION or its own HARNESS STATE, not a shared quantity: they carry the
# same value in every module that has one, and a module that stopped agreeing
# would fail loudly on the next read rather than drift quietly. Nothing here
# is a measurement, a threshold or a colour -- those must be imported.
PERMITTED_DUPLICATE_CONSTANTS = {
    "FAILURES",                       # each test module's own failure list
    "HERE", "REPOSITORY",             # each module's own path bootstrap
    "DOCUMENT", "MECHANISM_DOC",      # each folder's own document
    "DATASET_PATH", "MANIFEST_PATH",  # where the dataset is, stated per script
    "SHEET_DIR", "EXPORT_DIR", "CURVE_DIRECTORY",
    "TOLERANCE",                      # each verify script's own agreement bar
}

REPOSITORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _guarded_files():
    """
    Every module the duplicate rule covers: `data/`, the root, and the folders.

    IT COVERED `data/` ALONE UNTIL 2026-09-02, and the drift moved to where it
    could not see. Five copies of the document-comparison contract lived in the
    five `check_numbers.py` and no two were the same; the palettes and the
    figure wrapper were declared five times in the `build_figures.py`; and
    `_table` meant a memoised data frame in two folders and an HTML renderer in
    a third. None of it was visible to a guard globbing one directory.
    """
    return sorted(
        glob.glob(os.path.join(REPOSITORY, "data", "*.py"))
        + glob.glob(os.path.join(REPOSITORY, "*.py"))
        + glob.glob(os.path.join(REPOSITORY, "*", "build_figures.py"))
        + glob.glob(os.path.join(REPOSITORY, "*", "check_numbers.py")))


def _defined_names(path):
    """
    Top-level functions, classes AND CONSTANTS a module defines, by name.

    Constants were not covered until 2026-09-04, which is the whole category
    the guard exists for: the five palettes and the `TEMPERATURES` ramp that
    differed between folders are all constants, and one of those shadowings
    was still live in `temperature_series/build_figures.py` on the day this
    was widened. 27 constants were defined in more than one of the 58 modules
    and five had genuinely diverged.

    Three shapes are skipped, structurally rather than by name:

      - anything starting with `_`, which is private to its module;
      - `X = X`, which is a RE-EXPORT of an imported name and not a copy;
      - names that are not constant-shaped (lowercase module state).
    """
    tree = ast.parse(open(path).read())
    found = {node.name: node.lineno for node in tree.body
             if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef,
                                  ast.ClassDef))}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 \
                and isinstance(node.targets[0], ast.Name):
            name, value = node.targets[0].id, node.value
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            name, value = node.target.id, node.value
        else:
            continue
        if name.startswith("_") or not (name.isupper() or name[0].isupper()):
            continue
        if isinstance(value, ast.Name) and value.id == name:
            continue
        if name in PERMITTED_DUPLICATE_CONSTANTS:
            continue
        found.setdefault(name, node.lineno)
    return found


def _duplicate_names(files, root=None):
    """The names defined at top level in more than one of `files`."""
    root = root or REPOSITORY
    seen, duplicates = {}, {}
    for path in files:
        module = os.path.relpath(path, root)[:-3]
        for name, line in _defined_names(path).items():
            if name in PERMITTED_DUPLICATES:
                continue
            if name in seen:
                duplicates.setdefault(name, [seen[name]]).append(
                    f"{module}:{line}")
            else:
                seen[name] = f"{module}:{line}"
    return duplicates


def test_the_duplicate_guard_catches_a_planted_duplicate():
    """
    Fault injection, because a guard that has never failed is not a guard.

    This one was widened on 2026-09-02 from `data/` to the whole repository,
    and a widened glob that silently matched nothing would pass exactly as
    loudly as one that works.
    """
    print("\nthe duplicate guard, planted")
    import tempfile
    with tempfile.TemporaryDirectory() as root:
        first = os.path.join(root, "one.py")
        second = os.path.join(root, "two.py")
        open(first, "w").write("def lag_time():\n    return 1\n"
                                  "def main():\n    return 0\n")
        open(second, "w").write("def lag_time():\n    return 2\n"
                                   "def main():\n    return 0\n")
        found = _duplicate_names([first, second], root=root)
        check("a name defined in two modules is caught",
              "lag_time" in found, f"{found}")
        check("and the report names both places",
              len(found.get("lag_time", [])) == 2, f"{found}")
        check("a permitted name is not caught", "main" not in found)
        check("one module alone is not a duplicate",
              not _duplicate_names([first], root=root))

        # CONSTANTS, which the guard did not see at all until 2026-09-04 --
        # and which are the whole category it exists for. A diverged palette
        # is a constant; so was the `TEMPERATURES` ramp.
        third = os.path.join(root, "three.py")
        fourth = os.path.join(root, "four.py")
        open(third, "w").write('BURST_COLOUR = "#12856a"\n'
                               'DATASET_PATH = "data/experiment_data.csv"\n'
                               '_PRIVATE = 1\n'
                               'lowercase_state = 2\n')
        open(fourth, "w").write('BURST_COLOUR = "#7a4bb8"\n'
                                'DATASET_PATH = "data/experiment_data.csv"\n'
                                '_PRIVATE = 1\n'
                                'lowercase_state = 2\n'
                                'from x import BURST_COLOUR as _\n')
        found = _duplicate_names([third, fourth], root=root)
        check("a constant defined in two modules is caught",
              "BURST_COLOUR" in found, f"{found}")
        check("a permitted per-module path is not",
              "DATASET_PATH" not in found)
        check("nor is a private name", "_PRIVATE" not in found)
        check("nor is lowercase module state", "lowercase_state" not in found)

        # `X = X` is a re-export of an imported name, not a second copy.
        fifth = os.path.join(root, "five.py")
        open(fifth, "w").write('from curve_metrics import QUANTISATION_SIGMA\n'
                               'QUANTISATION_SIGMA = QUANTISATION_SIGMA\n')
        sixth = os.path.join(root, "six.py")
        open(sixth, "w").write('QUANTISATION_SIGMA = 0.000288\n')
        check("a re-export is not a duplicate of the thing it re-exports",
              "QUANTISATION_SIGMA" not in _duplicate_names([fifth, sixth],
                                                           root=root))


def test_every_test_is_actually_run():
    """
    A test that is defined but never called is worse than no test at all.

    It reads as coverage, it survives review, and it reports nothing. Six of
    them were found on 2026-09-04, five of them written in the days before
    while building the O2 correction: the suite said "0 failures" and had
    skipped every check on the repair it was there to guard. The sixth,
    `test_the_joint_buffer_order_reads_back_its_own_scheme`, had never run at
    all -- it was the planted-recovery test behind `+1.094 +/- 0.150`, and when
    it was finally called it failed.

    So the runner is checked rather than trusted. Every `test_*` defined at top
    level in a test module must be CALLED somewhere in that same module, which
    is what a `__main__` block does. This cannot be done by importing and
    running them -- that is the suite itself -- so it reads the source.
    """
    print("\nevery test is actually run")
    modules = sorted(glob.glob(os.path.join(REPOSITORY, "data", "test_*.py"))
                     + glob.glob(os.path.join(REPOSITORY, "test_*.py")))
    check("there are test modules to check", len(modules) >= 6,
          f"{len(modules)} modules")
    dead = []
    total = 0
    for path in modules:
        tree = ast.parse(open(path, encoding="utf-8").read())
        defined = [node.name for node in tree.body
                   if isinstance(node, ast.FunctionDef)
                   and node.name.startswith("test_")]
        called = {node.func.id for node in ast.walk(tree)
                  if isinstance(node, ast.Call)
                  and isinstance(node.func, ast.Name)}
        total += len(defined)
        dead += [f"{os.path.basename(path)}::{name}"
                 for name in defined if name not in called]
    check(f"every one of the {total} tests defined is called by its runner",
          not dead, "; ".join(dead))


def test_the_runner_finds_every_gate():
    """
    `run_gates.py` must discover every gate the repository actually has.

    It globs rather than listing, for the reason `test_every_test_is_run`
    exists one level down: a hardcoded list is a list that drifts, and the
    gate left off it is never run by anybody. CLAUDE.md's list named 9 and the
    repository has 20 -- `test_curve_flags`, `test_curve_screen`,
    `test_kinetic_model`, `test_read_rre`, `test_solution_chemistry` and
    `test_summary_kinetics` were in the tree and in no documented suite.

    A glob that silently stopped matching would otherwise pass by finding
    nothing to complain about, so this counts the files independently and
    requires the runner to have found all of them.
    """
    print("\nthe gate runner finds every gate")
    sys.path.insert(0, REPOSITORY)
    import run_gates

    found = set(run_gates.gate_paths(include_slow=True))
    expected = {"data/validate_dataset.py"}
    for pattern in ("test_*.py", os.path.join("data", "test_*.py"),
                    os.path.join("*", "check_numbers.py")):
        expected |= {os.path.relpath(p, REPOSITORY)
                     for p in glob.glob(os.path.join(REPOSITORY, pattern))}
    check(f"the runner finds all {len(expected)} gates",
          found == expected, f"missing {sorted(expected - found)}; "
                             f"extra {sorted(found - expected)}")
    check("and there are enough of them for the glob to be working",
          len(found) >= 20, f"{len(found)} gates")
    check("this module is one of them",
          "data/test_curve_metrics.py" in found)
    check("every folder's check_numbers is one of them",
          sum(1 for g in found if g.endswith("check_numbers.py")) == 7)

    # The slow suite is excluded from the routine run but must still EXIST --
    # `gate_paths` raises if it does not, so an optimiser suite cannot go
    # missing the way an unlisted gate would.
    routine = set(run_gates.gate_paths(include_slow=False))
    check("the slow suite is held back from the routine run",
          set(run_gates.SLOW_GATES) and not (routine & set(run_gates.SLOW_GATES)))
    check("but it is still required to exist",
          set(run_gates.SLOW_GATES) <= found)


def test_no_duplicate_definitions():
    """
    No name may be defined at top level in two modules, anywhere in the repo.

    A shared measurement belongs in `curve_metrics` and is imported; a shared
    drawing belongs in `figure_kit`, a shared document check in `doc_check`. If
    this fails, do not rename one of the copies -- delete one and import the
    other, or the two will drift apart exactly as they did before. Renaming is
    right only when the two are genuinely different things that happened to
    share a name, and then the new name has to say which one it is.
    """
    print("\nno duplicate definitions")
    files = _guarded_files()
    duplicates = _duplicate_names(files)
    check(f"no shared name is defined in two of the {len(files)} modules",
          not duplicates,
          "; ".join(f"{n} in {', '.join(w)}" for n, w in sorted(duplicates.items())))
    # The guard is worth nothing if it stopped covering something. These are
    # the three trees it exists for, and each has to be in the list.
    covered = {os.path.relpath(p, REPOSITORY) for p in files}
    for required in ("data/curve_metrics.py", "figure_kit.py", "doc_check.py",
                     "buffer/build_figures.py", "buffer/check_numbers.py"):
        check(f"the guard covers {required}", required in covered)
    # And that it still reads CONSTANTS, not only functions and classes. It
    # read only the latter until 2026-09-04, so a narrowing back would be
    # invisible: every folder would pass while a palette drifted.
    names = _defined_names(os.path.join(REPOSITORY, "figure_kit.py"))
    check("and it reads constants, not only functions and classes",
          {"RUNGS", "TEMPERATURES", "CATEGORY", "SURFACE"} <= set(names),
          f"missing {{'RUNGS','TEMPERATURES','CATEGORY','SURFACE'}} - {set(names)}")


def test_lag_statistic():
    """
    The published 160/402 has to come out of the canonical implementation.

    It has climbed twice, both times for the same reason and neither time
    because the statistic changed: the export rounds to 0.001 AU and that
    rounding flattens real lags below the threshold.

      136/402 (34%)    until 2026-08-31, all readings from the .txt exports
      151/402 (37.6%)  when the rate<n>.rre files were adopted
      158/402 (39.3%)  when mads_t<n>.rre was too -- 97 further curves
      160/402 (39.8%)  when the first reading of every run was dropped

    THAT LAST STEP IS SMALL AND THE STATISTIC IS NOT. Removing one reading
    flips this verdict on 46 of 402 curves; 24 gain a lag and 22 lose one, and
    the net +2 hides both. `acceleration` flips on 10. Quote this figure as
    "about 40%", never to three digits.

    See read_rre.py.
    """
    print("\nthe lag statistic")
    from fit_dataset import build_curves
    curves, _ = build_curves()
    positions = np.array([peak_position(c.absorbance, c.times) for c in curves],
                         dtype=float)
    lagging = int(np.nansum(positions > LAG_THRESHOLD))
    check("402 fittable curves", len(curves) == 402, f"got {len(curves)}")
    check("160 of them lag, as MECHANISM.md and FITTING.md report",
          lagging == 160, f"got {lagging}")

    # The shape the statistic is supposed to detect, and the shape it is not.
    times = np.linspace(0, 3000, 300)
    straight = 1e-5 * times
    sigmoid = 0.1 / (1 + np.exp(-(times - 2000) / 200))
    check("a straight line does not lag", peak_position(straight, times) == 0.0)
    check("a sigmoid whose inflection is late does lag",
          peak_position(sigmoid, times) > LAG_THRESHOLD)


def test_noise_and_rate():
    print("\nnoise and initial rate")
    rng = np.random.default_rng(0)
    times = np.linspace(0, 3000, 400)
    check("noise floors at the quantisation sigma on a smooth curve",
          curve_noise(1e-5 * times) == QUANTISATION_SIGMA)
    check("noise floors rather than dividing by nothing on a short curve",
          curve_noise([0.1, 0.2, 0.3]) == QUANTISATION_SIGMA)
    sigma = 0.004
    estimate = curve_noise(1e-5 * times + rng.normal(0, sigma, len(times)))
    check("noise recovers a known sigma through a strong trend",
          abs(estimate - sigma) < 0.25 * sigma, f"{estimate:.5f} vs {sigma}")

    slope, stderr, rms = initial_rate(times, 2e-5 * times)
    check("initial_rate recovers a known slope",
          abs(slope - 2e-5) < 1e-9, f"got {slope}")
    check("initial_rate returns (slope, stderr, rms)",
          np.isfinite([slope, stderr, rms]).all())
    check("initial_rate reads only the leading window",
          window_size(len(times), INITIAL_WINDOW) == int(len(times) * INITIAL_WINDOW))
    check("a floor keeps stderr finite on a perfectly straight window",
          stderr >= 0 and np.isfinite(stderr))
    check("line_slope is line_fit without the intercept",
          line_slope(times, 2e-5 * times) == line_fit(times, 2e-5 * times)[1:])


def test_acceleration():
    """
    The autocatalysis statistic, and the case peak_position gets wrong.

    A lag phase starts flat, and a flat start makes the first point-wise
    gradient a coin flip about zero -- which trips peak_position's
    `slope[0] <= 0` guard and scores the curve as having no lag at all. On the
    two-axis block that guard silences 31 of 96 live curves, including all six
    live curves of exp 142, whose lag is visible by eye. `acceleration` fits
    slopes over blocks instead, so it survives a flat start.
    """
    print("\nthe acceleration statistic")
    times = np.linspace(0, 1000, 80)

    z, _ = acceleration(times, 1e-4 * times)
    check("a straight line does not accelerate", abs(z) < ACCELERATION_SIGMA,
          f"z={z:.2f}")
    z, _ = acceleration(times, 0.05 * (1 - np.exp(-times / 200)))
    check("a saturating curve does not accelerate", z < 0, f"z={z:.2f}")
    z, where = acceleration(times, 0.05 / (1 + np.exp(-(times - 600) / 80)))
    check("a sigmoid accelerates", z > ACCELERATION_SIGMA, f"z={z:.2f}")
    check("and its steepest block sits near the inflection",
          0.4 < where < 0.8, f"where={where:.2f}")

    # A lag phase read at three decimals: exactly flat, then a ramp.
    flat_then_ramp = np.round(np.concatenate(
        [np.zeros(40), 4e-5 * (times[40:] - times[40])]), 3)
    z, _ = acceleration(times, flat_then_ramp)
    check("a flat start does not defeat it", z > ACCELERATION_SIGMA,
          f"z={z:.2f}")
    check("where peak_position's first-point guard does defeat it",
          peak_position(flat_then_ramp, times) == 0.0,
          "guard no longer fires -- peak_position changed, revisit this test")

    check("too few points for two blocks returns nan",
          not np.isfinite(acceleration(times[:6], times[:6] * 1e-4)[0]))


def test_floor_belongs_to_the_source():
    """
    The variance floor is an argument, and it changes the verdict.

    `curve_noise` has taken its floor as an argument since 2026-08-31, but
    `line_fit` hardcoded QUANTISATION_SIGMA until 2026-09-01 -- so every
    standard error in the package, and the acceleration z-score that divides
    by two of them, was floored at the .txt export's 0.001 AU rounding even on
    .rre curves read a thousand times finer. It bound on 52 of the 110 live
    two-axis curves and cost 3 of them their acceleration verdict (48/110 read
    where the instrument says 51/110).

    These checks fail if the floor is ever hardcoded again.
    """
    print("\nthe floor is a property of the source")
    times = np.linspace(0, 1000, 80)
    fine = RRE_SIGMA

    # A signal below the export's quantisation but far above the .rre's: real
    # to the instrument, invisible to the export. This is the regime the whole
    # 2026-08-31 .rre swap was about.
    straight = 1e-7 * times

    _, _, coarse_stderr, _ = line_fit(times, straight)
    _, _, fine_stderr, _ = line_fit(times, straight, fine)
    check("line_fit takes a floor and a smaller one gives a smaller stderr",
          fine_stderr < coarse_stderr,
          f"{fine_stderr:.3e} vs {coarse_stderr:.3e}")
    check("and the ratio is the ratio of the floors, since both are floored",
          abs(coarse_stderr / fine_stderr
              - QUANTISATION_SIGMA / RRE_SIGMA) < 1e-6 * (QUANTISATION_SIGMA / RRE_SIGMA),
          f"ratio {coarse_stderr / fine_stderr:.1f}")

    check("line_fit's default is still the export's floor",
          line_fit(times, straight) == line_fit(times, straight,
                                                QUANTISATION_SIGMA))

    # The verdict itself moves: a rise of a few .rre quanta accelerates when
    # judged against the instrument's own floor and does not when judged
    # against the export's.
    # Total rise 1e-4 AU: a tenth of one export quantum, ~380 .rre quanta.
    ramp = np.concatenate([np.zeros(40), 2e-7 * (times[40:] - times[40])])
    check("a sub-quantum acceleration is invisible at the export's floor",
          acceleration(times, ramp)[0] < ACCELERATION_SIGMA,
          f"z={acceleration(times, ramp)[0]:.2f}")
    check("and visible at the instrument's",
          acceleration(times, ramp, floor=fine)[0] > ACCELERATION_SIGMA,
          f"z={acceleration(times, ramp, floor=fine)[0]:.2f}")

    check("source_floor maps the sources to the two constants",
          source_floor("rre") == RRE_SIGMA
          and source_floor("txt") == QUANTISATION_SIGMA)


def test_whole_curve_estimators_return_rates():
    """
    `whole_slope` and `quadratic_rate` must return a RATE, not an intercept.

    This test exists because they did not. `whole_slope` was written as
    `line_fit(...)[:2]`, and `line_fit` returns (intercept, slope, stderr,
    rms), so it handed back the intercept -- an absorbance -- as the rate, and
    the slope as its standard error. It shipped in c41f459 and produced a
    plausible-looking number that was written up as a finding before a figure
    made it obvious. See DATA_VERIFICATION.md, 2026-09-01.

    The check is that a line of KNOWN slope and a large NON-ZERO offset comes
    back as the slope: an intercept-returning implementation fails only if the
    offset is not zero, which is exactly why the bug was invisible on
    baseline-subtracted curves that start near zero.
    """
    print("\nwhole-curve estimators return rates")
    times = np.linspace(0.0, 1000.0, 200)
    slope, offset = 4.0e-5, 0.35
    values = offset + slope * times

    measured, stderr = whole_slope(times, values)
    check("whole_slope returns the slope of a straight line",
          abs(measured - slope) < 1e-12, f"got {measured:.4e}, want {slope:.4e}")
    check("...and not its intercept",
          abs(measured - offset) > abs(measured) * 0.5,
          f"got {measured:.4e}, which is the offset {offset}")
    check("whole_slope's second return is an error, not the slope",
          stderr < abs(slope), f"stderr {stderr:.3e} vs slope {slope:.3e}")

    v0, v0_stderr, curvature_t = quadratic_rate(times, values)
    check("quadratic_rate returns the slope on a straight line",
          abs(v0 - slope) < 1e-10, f"got {v0:.4e}, want {slope:.4e}")
    check("...and reports no curvature there",
          abs(curvature_t) < 3.0, f"curvature t {curvature_t:.2f}")

    # A decelerating curve: the initial rate must exceed the average rate.
    bend = slope * times - 1.2e-8 * times ** 2
    initial, *_ = quadratic_rate(times, offset + bend)
    average, _ = whole_slope(times, offset + bend)
    check("on a decelerating curve the initial rate exceeds the average",
          initial > average > 0, f"initial {initial:.3e}, average {average:.3e}")
    check("and the curvature is flagged negative",
          quadratic_rate(times, offset + bend)[2] < -3.0,
          f"curvature t {quadratic_rate(times, offset + bend)[2]:.1f}")


def test_outlier_flagging():
    """
    A single spike is found; a real kinetic feature is not.

    The distinction is a timescale argument. At 30-60 s sampling nothing
    chemical moves in one interval and reverts in the next, so ONE reading out
    of line with both neighbours is an artefact. Two or more consecutive ones
    are not separated from chemistry, and `curve_screen` is explicit that curve
    shape is never a defect -- so `isolated_outliers` reports the two
    separately and neither is ever removed here.
    """
    print("\noutlier flagging")
    noise = 3e-4
    times = np.linspace(0.0, 1200.0, 40)
    clean = 4e-5 * times
    rng = np.random.default_rng(11)
    noisy = clean + rng.normal(0, noise, len(times))

    isolated, in_runs = isolated_outliers(times, noisy, noise)
    check("a clean noisy line flags nothing",
          len(isolated) == 0 and len(in_runs) == 0,
          f"isolated {list(isolated)}, runs {list(in_runs)}")

    spiked = noisy.copy()
    spiked[17] += 12 * noise
    isolated, in_runs = isolated_outliers(times, spiked, noise)
    check("an interior spike is found and called isolated",
          list(isolated) == [17] and len(in_runs) == 0,
          f"isolated {list(isolated)}, runs {list(in_runs)}")

    # The endpoint case, which is the whole reason the window is one-sided:
    # the first reading is the archive's worst-behaved and carries the most
    # leverage on v0, so it has to be SCORABLE at all.
    first = noisy.copy()
    first[0] -= 14 * noise
    scores = local_outlier_z(times, first, noise)
    check("a bad FIRST reading is scored, not skipped",
          np.isfinite(scores[0]) and abs(scores[0]) > OUTLIER_SIGMA,
          f"z[0] = {scores[0]:.1f}")
    # ...but it drags its neighbour past the threshold, so the pair reads as a
    # RUN and not as an isolated spike. This is exactly why scope.frame takes
    # `first_point_flagged` from z[0] rather than from `isolated`.
    isolated, in_runs = isolated_outliers(times, first, noise)
    check("a bad first reading can drag its neighbour into a run",
          list(in_runs) == [0, 1] and len(isolated) == 0,
          f"isolated {list(isolated)}, runs {list(in_runs)}")

    # A GRADUAL induction period must survive, which is the real case: exp 65
    # sample 3's flat-to-rise spans several readings and scores -2.7.
    ramp = np.clip((times - times[18]) / (times[24] - times[18]), 0, 1)
    gradual = 8e-5 * ramp * (times - times[18]) + rng.normal(0, noise, len(times))
    isolated, _ = isolated_outliers(times, gradual, noise)
    check("a gradual flat-then-rise transition is not called an artefact",
          len(isolated) == 0, f"isolated {list(isolated)}")

    # A genuinely INSTANTANEOUS kink is flagged. Stated rather than asserted
    # away: it is a real limitation, and one more reason nothing is removed
    # automatically.
    # Its own generator, so this does not depend on how many draws the checks
    # above happened to consume.
    kink_noise = np.random.default_rng(5).normal(0, noise, len(times))
    kinked = np.concatenate([np.zeros(18),
                             1.6e-4 * (times[18:] - times[18])]) + kink_noise
    scores = local_outlier_z(times, kinked, noise)
    check("an instantaneous kink IS flagged -- a known limitation",
          abs(scores[18]) > OUTLIER_SIGMA, f"z[18] = {scores[18]:.1f}")

    # Degree 2 is load bearing: a local LINE reads a GRADUAL transition as an
    # outlier, which is how an automatic filter deletes real chemistry.
    straight_fit = local_outlier_z(times, gradual, noise, degree=1)
    curved_fit = local_outlier_z(times, gradual, noise, degree=2)
    check("a local line scores a gradual transition worse than a quadratic",
          np.nanmax(np.abs(straight_fit)) > np.nanmax(np.abs(curved_fit)),
          f"line {np.nanmax(np.abs(straight_fit)):.1f} vs "
          f"quadratic {np.nanmax(np.abs(curved_fit)):.1f}")

    # Masking: adjacent spikes sit in each other's fitting windows and pull
    # them toward themselves, so the second often falls under the threshold.
    two = noisy.copy()
    two[20] += 9 * noise
    two[21] += 9 * noise
    scores = local_outlier_z(times, two, noise)
    check("adjacent spikes mask each other, so both scores shrink",
          abs(scores[20]) < 9 and abs(scores[21]) < 9,
          f"z[20] = {scores[20]:.1f}, z[21] = {scores[21]:.1f}")
    check("...and at least one of them is still caught",
          max(abs(scores[20]), abs(scores[21])) > OUTLIER_SIGMA,
          f"z[20] = {scores[20]:.1f}, z[21] = {scores[21]:.1f}")


def test_model_residual():
    """
    One residual definition for every form, and it says what `bounded` cannot.

    The point of the function is that it separates "is this parameter pinned"
    from "does this form fit", which came apart on exp 65 when the burst/lag v0
    was proposed as the headline: bounded v0 on a fit sitting 7-8x above noise.
    """
    print("model residual")
    rng = np.random.default_rng(11)
    times = np.arange(0.0, 600.0, 10.0)
    noise = 1e-4
    truth = 2.0e-5 * times

    check("a perfect fit scores 0", model_residual(truth, truth, 3, noise) == 0.0)

    scattered = truth + rng.normal(0.0, noise, len(times))
    at_noise = model_residual(scattered, truth, 3, noise)
    check("scatter at the noise scores about 1", 0.8 < at_noise < 1.2,
          f"{at_noise:.3f}")

    # The parameter count matters: the same residuals over fewer degrees of
    # freedom score higher, which is what makes a 3-parameter and a
    # 4-parameter form comparable at all.
    check("more parameters, larger residual",
          model_residual(scattered, truth, 4, noise) >
          model_residual(scattered, truth, 3, noise))

    # A form that misses by a constant offset is caught in units of noise,
    # which is the property the panels rely on.
    offset = model_residual(scattered + 5 * noise, truth, 3, noise)
    check("a 5-sigma offset scores about 5", 4.5 < offset < 5.6, f"{offset:.2f}")

    check("zero noise is nan, not a divide", np.isnan(
        model_residual(scattered, truth, 3, 0.0)))

    # And the real case the function was added for.
    frame = scope.frame((65,))
    check("exp 65's burst fits are worse than its quadratic ones",
          bool((frame.v0_burst_resid > frame.v0_quad_resid).all()),
          f"burst {frame.v0_burst_resid.round(1).tolist()} vs "
          f"quad {frame.v0_quad_resid.round(1).tolist()}")
    check("and every one of them is beyond 3x noise",
          bool((frame.v0_quad_resid > 3).all()))
    # bounded says the parameter is pinned; the residual says the form is
    # wrong. If these ever agree on exp 65 the ANALYSIS.md argument for
    # keeping v0_quad as the headline needs rewriting, not just re-running.
    check("yet all four report a bounded v0",
          bool(frame.v0_burst_bounded.all()))


def test_segmented_fit():
    """
    The two-line split, and the shape every other statistic here steps over.

    `late_over_early`, `acceleration` and `peak_position` all compare a curve's
    START to its END. A curve that breaks upward in the MIDDLE and plateaus
    looks ordinary to all three -- exp 65 ranked mid-pack on `late_over_early`
    while carrying the most distinctive shape in its block. This test pins both
    halves: that the split recovers a known break, and that it separates exp 65
    from the run it is compared against.
    """
    print("segmented fit")
    times = np.arange(0, 40) * 10.0

    # An exact two-slope curve: the break and both slopes come back exactly.
    values = np.where(times < 200, 1e-4 * times,
                      1e-4 * 200 + 5e-4 * (times - 200))
    where, before, after, ratio = segmented_fit(times, values)
    check("recovers the break time", abs(where - 200) <= 10)
    check("recovers the slope ratio", abs(ratio - 5.0) < 0.01)

    # A straight line must not invent a break: some split always minimises the
    # residual, so the guard is the RATIO, not the existence of a breakpoint.
    _, _, _, straight = segmented_fit(times, 2e-4 * times)
    check("a straight line gives ratio 1", abs(straight - 1.0) < 0.01)

    # Too short to split is nan, not a coincidence of five points.
    check("too few points is nan",
          not np.isfinite(segmented_fit(times[:6], values[:6])[3]))

    # The real case. Exp 65's four cuvettes steepen across a shared break;
    # exp 67, matched to it in substrate and peroxide, decelerates. If this
    # ever stops separating them, ANALYSIS.md section 6b is out of date.
    frame = scope.frame((65, 67))
    frame = frame[frame.live]
    boric = frame[frame.experiment == 65]
    phosphate = frame[frame.experiment == 67]
    check("all four boric cuvettes steepen",
          bool((boric.break_ratio > SEGMENT_RATIO_STEEP).all()))
    check("no phosphate cuvette does",
          bool((phosphate.break_ratio <= SEGMENT_RATIO_STEEP).all()))
    check("and the boric breaks are synchronised to two sampling intervals",
          float(boric.break_time.max() - boric.break_time.min()) <= 56.0)


def test_two_breakpoints():
    """
    Two breakpoints are found when there are two, and refused when there is one.

    `segmented_fit` searches for ONE break, so on a curve whose rate rises and
    then falls it lands on whichever change is stronger and never reports the
    other. That is not a tuning problem, it is the shape of the search, and it
    is why the early break on the 40 C curves was invisible until it was seen
    by eye.
    """
    print("two breakpoints")
    times = np.arange(0, 180) * 40.0

    # rise then fall: slopes 1e-5, 4e-5, 2e-5, breaking at 2000 and 4800 s
    def piecewise(edges, slopes):
        out = np.zeros(len(times))
        level = 0.0
        previous = 0.0
        for edge, slope in zip(list(edges) + [times[-1] + 1], slopes):
            inside = (times >= previous) & (times < edge)
            out[inside] = level + slope * (times[inside] - previous)
            level += slope * (min(edge, times[-1] + 1) - previous)
            previous = edge
        return out

    generator = np.random.default_rng(1)
    curve = piecewise((2000.0, 4800.0), (1e-5, 4e-5, 2e-5))
    noisy = curve + generator.normal(0, 3e-4, len(times))
    result = segment_selection(times, noisy)
    check("a rise-then-fall curve takes two breakpoints",
          result["breaks"] == 2, f"{result['breaks']}, F {result['f_statistic']:.1f}")
    check("and places them near the truth",
          len(result["times"]) == 2
          and abs(result["times"][0] - 2000) < 600
          and abs(result["times"][1] - 4800) < 600,
          f"{[round(v) for v in result['times']]}")
    check("and names the pattern", result["pattern"] == "rise then fall",
          result["pattern"])

    # one break, monotone: the pattern must not claim a maximum
    single = piecewise((2500.0,), (1e-5, 3e-5))
    result = segment_selection(times, single + generator.normal(0, 3e-4, len(times)))
    check("a single bend is called rising", result["pattern"] == "rising",
          result["pattern"])

    # The one-break search is a special case and must not beat the two-break
    # search: more freedom cannot fit worse.
    for values in (noisy, single):
        _, _, one = segment_breaks(times, values, 1)
        _, _, two = segment_breaks(times, values, 2)
        check("two breakpoints never fit worse than one",
              two <= one * 1.0001, f"{two:.3e} against {one:.3e}")

    # And the prefix-sum errors must equal an honest least squares.
    errors = _segment_errors(times, noisy, 4)
    for start, stop in ((0, 40), (17, 96), (100, len(times))):
        design = np.column_stack([np.ones(stop - start), times[start:stop]])
        beta, *_ = np.linalg.lstsq(design, noisy[start:stop], rcond=None)
        residual = noisy[start:stop] - design @ beta
        check("the prefix-sum stretch error matches a direct fit",
              abs(errors[start, stop] - float(residual @ residual))
              <= 1e-9 * float(residual @ residual),
              f"{errors[start, stop]:.6e} against {float(residual @ residual):.6e}")


def _sawtooth(times, chemistry, edges, rate, ends_holding=True):
    """
    A curve plus a bubble that grows at `rate` and sheds at each edge.

    `ends_holding` is what the run is doing when the recording stops: True
    leaves production running past the last release, so the curve ends
    carrying a bubble that never detached, and False stops it there. The
    reconstruction subtracts only gas it watched leave, so the two endings are
    the two sides of that rule and both are planted below.
    """
    artefact = np.zeros(len(times))
    start = 0
    for edge in edges:
        artefact[start:edge + 1] = rate * (times[start:edge + 1] - times[start])
        start = edge + 1
    if ends_holding:
        artefact[start:] = rate * (times[start:] - times[start])
    return chemistry + artefact, artefact


def _stitch(values, drops):
    """The repair this module refuses: add each step back to everything after
    it. Defined HERE and not in curve_metrics, because it is the wrong answer
    and the only thing it is wanted for is to fail against the right one."""
    out = np.asarray(values, dtype=float).copy()
    steps = np.diff(values)
    for index in drops:
        out[index + 1:] -= steps[index]
    return out


def test_the_bubble_correction():
    print("\nthe bubble correction")
    times = np.arange(0, 3600, 60.0)
    noise = 1e-4
    # A decelerating reaction: the chemistry the repairs have to give back.
    chemistry = 0.08 * (1 - np.exp(-times / 1200.0))
    edges = np.array([12, 26, 41, 52])

    spoilt, artefact = _sawtooth(times, chemistry, edges,
                                 rate=1.0 * chemistry[-1] / times[-1])
    found = bubble_drops(spoilt, noise)
    check("every planted detachment is found",
          set(found) == set(edges), f"{sorted(found)} against {list(edges)}")
    check("a curve with no bubble has no detachment",
          len(bubble_drops(chemistry, noise)) == 0)

    events = detachments(spoilt, noise)
    check("each planted detachment is one event, not several",
          [start for start, _ in events] == list(edges),
          f"{events}")

    corrected, _ = debubble(times, spoilt, noise)
    before = np.abs(spoilt - chemistry).max()
    after = np.abs(corrected - chemistry).max()

    # ONLY GAS THAT WAS WATCHED TO LEAVE IS SUBTRACTED, so this planting --
    # which is still making gas when the recording stops -- keeps the bubble
    # it never shed, and what is left over is exactly that bubble and nothing
    # else. Asserting the residual rather than bounding it is the point: an
    # over-correction of the same size would pass a bound and fails this.
    unshed = np.zeros(len(times))
    unshed[edges[-1] + 1:] = artefact[edges[-1] + 1:]
    # To a tenth of a bubble, not to the noise: the rate is pinned to the LEAST
    # that pays for every detachment, so it comes out a few percent under the
    # planted one and the reconstruction removes slightly less than was put in.
    # That is the direction the whole model is built to err in.
    check("what the correction leaves is the bubble that never detached",
          np.abs(corrected - chemistry - unshed).max() < 0.1 * artefact.max(),
          f"worst {np.abs(corrected - chemistry - unshed).max():.2e} "
          f"against a bubble of {artefact.max():.4f}")
    check("so it is closer than the readings, but not by everything",
          0.05 * before < after < 0.75 * before,
          f"{after:.2e} against {before:.2e}")

    # The same planting stopped at its last release: nothing is left holding
    # gas, and there the repair is exact.
    shed, _ = _sawtooth(times, chemistry, edges,
                        rate=1.0 * chemistry[-1] / times[-1],
                        ends_holding=False)
    shed_fixed, _ = debubble(times, shed, noise)
    check("a run that stops making gas is corrected in full",
          np.abs(shed_fixed - chemistry).max() < 0.1
          * np.abs(shed - chemistry).max(),
          f"{np.abs(shed_fixed - chemistry).max():.2e} against "
          f"{np.abs(shed - chemistry).max():.2e}")
    check("and the tail of one that does not is left ON the readings",
          abs(corrected[-1] - spoilt[-1]) < 1e-12,
          f"{corrected[-1] - spoilt[-1]:+.2e}")

    # THE PROPERTY THE OLD SEGMENT RAMP DID NOT HAVE, and the reason this
    # module was rewritten: the repaired curve may not fall. `bubble_profile`
    # lets the gas grow by at most what the reading itself gained, so
    # `A_obs - b` is non-decreasing at every ordinary step, and pays for each
    # detachment in full, so it is non-decreasing across those too.
    check("the reconstruction is non-decreasing everywhere",
          bool((np.diff(corrected) >= -1e-15).all()),
          f"worst step {np.diff(corrected).min():+.2e}")
    check("and the readings it was built from are not",
          np.diff(spoilt).min() < -8 * noise,
          f"worst step {np.diff(spoilt).min():+.2e}")

    # THE CASE AGAINST STITCHING, as an assertion rather than a caveat.
    stitched = _stitch(spoilt, found)
    check("stitching ends ABOVE the truth by the sum of the drops",
          stitched[-1] - stitched[0] > 1.5 * (chemistry[-1] - chemistry[0]),
          f"{stitched[-1] - stitched[0]:.4f} against "
          f"{chemistry[-1] - chemistry[0]:.4f}")
    check("stitching is further from the truth than doing nothing",
          np.abs(stitched - chemistry).max()
          > np.abs(spoilt - chemistry).max())
    check("the reconstruction is closer than either", after < before
          and after < np.abs(stitched - chemistry).max())
    check("stitching IS this model at rate zero",
          np.allclose(spoilt - bubble_profile(times, spoilt, events, 0.0)
                      + np.cumsum(np.concatenate(
                          [[0.0], -np.diff(spoilt) * np.isin(
                              np.arange(len(times) - 1), found)])),
                      stitched))

    # A bubble starts at nothing and ends at nothing or more, so the repair can
    # only ever LOWER the net rise. This is the mass balance stitching breaks.
    check("the correction never raises the net rise",
          (corrected[-1] - corrected[0]) <= (spoilt[-1] - spoilt[0]) + 1e-12,
          f"{corrected[-1] - corrected[0]:.4f} against "
          f"{spoilt[-1] - spoilt[0]:.4f}")

    held = bubble_profile(times, spoilt, events, bubble_rate(
        times, spoilt, events))
    check("the gas is never negative", held.min() >= 0.0)
    check("the gas never outruns the curve it rides on",
          bool((np.diff(held) <= np.maximum(np.diff(spoilt), 0.0)
                + 1e-15).all()))
    check("the rate pays for every detachment and no more",
          abs(bubble_shortfall(times, spoilt, events,
                               bubble_rate(times, spoilt, events))) < 1e-9)
    check("a rate below it cannot pay",
          bubble_shortfall(times, spoilt, events,
                           0.5 * bubble_rate(times, spoilt, events)) > 0)

    # A bubble that grew before the first reading leaves no rise to date it
    # from, so no rate explains it and the curve is returned untouched. The
    # fall has to be well clear of the excursion test -- on this chemistry the
    # first interval alone rises 0.0039, so a 0.01 fall is 61% recovered by
    # the very next reading and is rejected as a spike, correctly.
    early = chemistry.copy()
    early[1:] -= 0.03
    untouched, _ = debubble(times, early, noise)
    check("a detachment in the first interval has no affordable rate",
          not np.isfinite(bubble_rate(times, early, detachments(early, noise))))
    # THE FALLS MODEL leaves this curve exactly alone -- there is no rate to
    # apply. `bubble_gains` is not quite as clean beside a fall this extreme:
    # a 0.03 AU, ~120 sigma drop in the very first interval is far outside
    # anything the real archive carries, and it distorts `local_outlier_z`'s
    # local fit for a few readings after it the same way `isolated_outliers`
    # documents for two adjacent real spikes ("masking") -- z climbs past
    # +30 sigma at reading 4 from the fit trying to bend around the drop, not
    # from anything arriving. The whole-archive sweep in this file finds no
    # such case on a real curve; this is what a fall four orders of magnitude
    # past anything real does to a test built for something else.
    check("and the readings the falls model leaves alone are moved by only "
          "a masking artefact of the planted 120 sigma fall, not a real "
          "gain",
          np.abs(untouched - early).max() < 0.001,
          f"{np.abs(untouched - early).max():.6f}")

    # ONE BUBBLE, TWO READINGS. A fall spread over consecutive readings is one
    # detachment; counting it as two gave the second a growth window of zero
    # seconds and left the whole of it uncorrected.
    slow = spoilt.copy()
    slow[28:] -= 0.004
    pair = detachments(slow, noise)
    check("a fall over two readings is one detachment",
          (26, 28) in pair, f"{pair}")
    rebuilt, _ = debubble(times, slow, noise)
    check("and it is corrected in full",
          bool((np.diff(rebuilt) >= -1e-15).all()),
          f"worst step {np.diff(rebuilt).min():+.2e}")

    load = bubble_load(spoilt, found)
    check("the load is the absorbance lost over the net rise",
          abs(load - (-np.diff(spoilt)[found].sum()
                      / (spoilt[-1] - spoilt[0]))) < 1e-12,
          f"{load:.3f}")
    check("a clean curve carries no load",
          bubble_load(chemistry, bubble_drops(chemistry, noise)) == 0.0)
    check("a curve that went nowhere has no load to report",
          not np.isfinite(bubble_load(np.zeros(20), np.array([], dtype=int))))
    check("a clean curve is returned untouched, not merely close",
          np.array_equal(debubble(times, chemistry, noise)[0], chemistry))


def test_the_excursion_test_on_planted_spikes():
    """
    The excursion test on synthetic spikes of known kind.

    Gas that leaves the beam does not come back, and a bubble does not grow
    half its size in one 60 s reading. So a fall flanked by a single reading
    that climbs a comparable amount is an instrument spike, and correcting it
    as gas removes real chemistry: on exp 149 cuvette 5 that cost 0.0097 AU
    off a curve that rose 0.0262.
    """
    print("\na fall that comes straight back is not gas")
    times = np.arange(0, 3600, 60.0)
    noise = 1e-4
    chemistry = 0.02 + 1e-5 * times

    # A SPIKE DOWN: one reading dips and the next is back on the line.
    spike = chemistry.copy()
    spike[30] -= 0.004
    check("the detector sees the fall",
          30 - 1 in set(bubble_drops(spike, noise)) or 29 in
          set(bubble_drops(spike, noise)),
          f"{list(bubble_drops(spike, noise))}")
    check("but it is not a detachment", detachments(spike, noise) == [],
          f"{detachments(spike, noise)}")
    check("and the curve is returned untouched",
          np.array_equal(debubble(times, spike, noise)[0], spike))

    # A SPIKE UP: one reading jumps and the fall is the return off it.
    perched = chemistry.copy()
    perched[30] += 0.004
    check("a fall off a spike is not a detachment either",
          detachments(perched, noise) == [], f"{detachments(perched, noise)}")

    # A REAL STEP: the level goes down and STAYS down.
    step = chemistry.copy()
    step[30:] -= 0.004
    check("a persistent step IS a detachment",
          detachments(step, noise) == [(29, 30)],
          f"{detachments(step, noise)}")
    rebuilt, _ = debubble(times, step, noise)
    check("and it is corrected in full",
          rebuilt[30] - rebuilt[29] >= -1e-12,
          f"{rebuilt[30] - rebuilt[29]:+.2e}")

    # THE TEST MUST NOT LOOK ACROSS THE FALL. `local_outlier_z` does, so a
    # genuine step flags itself -- which is why this is a separate test and
    # not a call to that one.
    z = local_outlier_z(times, step, noise)
    check("local_outlier_z flags the step change itself",
          abs(z[29]) > OUTLIER_SIGMA or abs(z[30]) > OUTLIER_SIGMA,
          f"z[29]={z[29]:+.1f} z[30]={z[30]:+.1f}")


def test_the_monotone_bound():
    print("\nthe monotone bound")
    times = np.arange(0, 3600, 60.0)
    chemistry = 0.08 * (1 - np.exp(-times / 1200.0))
    spoilt, _ = _sawtooth(times, chemistry, np.array([12, 26, 41, 52]),
                          rate=1.0 * chemistry[-1] / times[-1])
    bound = monotone_bound(spoilt)
    check("the bound never decreases", (np.diff(bound) >= -1e-15).all())
    check("the bound lies under the readings", (bound <= spoilt + 1e-15).all())
    check("and it is the GREATEST such function -- it touches the readings",
          np.isclose(bound, spoilt).any())
    # The bound is an upper bound on the chemistry: A_chem <= min of every
    # later reading, because product accumulates and gas only ever adds.
    check("the bound stays above the chemistry it brackets",
          (bound >= chemistry - 1e-12).all())
    check("a clean rising curve is its own bound",
          np.allclose(monotone_bound(chemistry), chemistry))


def test_the_bubble_the_run_never_shed():
    print("\nthe bubble the run never shed")
    times = np.arange(0, 3600, 60.0)
    noise = 1e-4
    chemistry = 0.08 * (1 - np.exp(-times / 1200.0))
    edges = np.array([12, 26, 41, 44])
    rate = 1.0 * chemistry[-1] / times[-1]

    # STILL MAKING GAS WHEN THE RECORDING STOPS. `debubble` hands the tail
    # back to the readings, so what it leaves behind is exactly the bubble
    # that never detached -- and this bounds it.
    spoilt, artefact = _sawtooth(times, chemistry, edges, rate=rate)
    corrected, events = debubble(times, spoilt, noise)
    planted = float(artefact[-1])
    fitted_rate = bubble_rate(times, spoilt, events)
    held, stripped = terminal_gas(times, corrected, events, fitted_rate)
    # NOT BOUNDED FROM ABOVE IN GENERAL -- that was this check's original
    # name and claim, and it held here only because edges 41 and 44 (3
    # readings apart) used to collapse to 3 detected events instead of 4:
    # the OLD `_is_excursion`, comparing a recovering step only to the SIZE
    # OF ITS OWN DROP, wrongly read edge 44's recovery as a spike, the same
    # bug exp 130 cuvette 2 was caught on (`_local_step_scale`'s docstring).
    # Fixing that gives the TRUE 4-event detection and a shorter, correct
    # final span (840 s from event 44's stop, not the wrong 900 s+ from 41's).
    # `held = fitted_rate * span` over that span, uncapped by the tail-rise
    # clause here, so `held` inherits `bubble_rate`'s own shortfall from the
    # PLANTED rate EXACTLY -- 7.4% under, both here -- rather than being
    # bounded above it. `terminal_gas` was never entitled to more than that:
    # it is built on the same "least rate that pays" `bubble_rate` is, and
    # can only ever be as accurate as that rate is, not more. This asserts
    # the relationship that is actually true, not the stronger one the old
    # name claimed.
    check("the terminal bound tracks the fitted rate's own shortfall, not more",
          abs(held - planted * (fitted_rate / rate)) < 1e-9,
          f"held={held:.6f} planted*(fitted/true rate)="
          f"{planted * (fitted_rate / rate):.6f}")
    check("and not by more than half again",
          held < 1.5 * planted, f"{held:.4f} against {planted:.4f}")
    check("taking it off brings the tail back towards the chemistry",
          np.abs(stripped[edges[-1] + 1:] - chemistry[edges[-1] + 1:]).max()
          < 0.5 * np.abs(corrected[edges[-1] + 1:]
                         - chemistry[edges[-1] + 1:]).max(),
          f"{np.abs(stripped[edges[-1] + 1:] - chemistry[edges[-1] + 1:]).max():.2e}"
          f" against {np.abs(corrected[edges[-1] + 1:] - chemistry[edges[-1] + 1:]).max():.2e}")
    check("the strip never lifts the curve",
          bool((stripped <= corrected + 1e-15).all()))
    check("and never takes off more than the tail rose",
          held <= corrected[-1] - corrected[events[-1][1]] + 1e-15)
    check("a curve with no detachment is untouched",
          terminal_gas(times, chemistry, [], 0.0)[0] == 0.0)
    check("and so is one whose rate could not be pinned",
          terminal_gas(times, corrected, events, np.inf)[0] == 0.0)

    # THE BOUND CANNOT TELL THE TWO ENDINGS APART -- it asks the fitted rate,
    # not the readings -- so it charges the run that stopped making gas too.
    # That is why it is a bracket and not a repair, and why `tail_excess`
    # exists.
    shed, _ = _sawtooth(times, chemistry, edges, rate=rate,
                        ends_holding=False)
    shed_fixed, shed_events = debubble(times, shed, noise)
    quiet, _ = terminal_gas(times, shed_fixed, shed_events,
                            bubble_rate(times, shed, shed_events))
    check("a run that stopped making gas is charged all the same",
          quiet > 0.0, f"{quiet:.4f} where the truth is 0")

    # AND THIS IS WHAT SEPARATES THEM. The two plantings are the same
    # chemistry and differ only by the bubble that never left, so the gap
    # between their tail excesses IS the gas rate -- to 10% of it here.
    holding = tail_excess(times, corrected, events)
    stopped = tail_excess(times, shed_fixed, shed_events)
    check("the gap between the two endings measures the gas rate",
          0.7 * rate < holding - stopped < 1.3 * rate,
          f"{holding - stopped:.3e} against a planted {rate:.3e}")

    # A ONE-SIDED TEST, AND THE PLANTING SAYS WHICH SIDE. This chemistry
    # DECELERATES, so its tail is flatter than its body whether or not a
    # bubble is growing in it, which pulls the excess of a run that is STILL
    # making gas back down towards zero rather than leaving it positive.
    # Before the `_is_excursion` fix this landed comfortably negative
    # (-2e-5-ish, on the wrong 3-event span); on the TRUE 4-event span it
    # lands at +1.9e-6 -- under 9% of the planted gas rate, and the sign is
    # not the point. What matters, and is asserted here, is that the
    # deceleration confound pulls the reading down NEAR zero and far below
    # the clean, unambiguous positive signal an accelerating chemistry gives
    # below (`tail_excess > 0` with no gas in the tail at all) -- so a small
    # or negative excess is not evidence the run stopped making gas, it is
    # the absence of evidence that it did not, and `terminal_gas`'s bound
    # stays standing and uncredited either way.
    check("a decelerating curve hides the gas signal near zero, not evidence against it",
          abs(holding) < 0.2 * abs(rate),
          f"{holding:+.2e} against a planted gas rate of {rate:.2e}")
    speeding = 0.02 * (np.exp(times / 1800.0) - 1.0)
    fast, _ = _sawtooth(times, speeding, edges, rate=rate, ends_holding=False)
    fast_fixed, fast_events = debubble(times, fast, noise)
    check("and an accelerating one shows a bubble that is not there",
          tail_excess(times, fast_fixed, fast_events) > 0,
          f"{tail_excess(times, fast_fixed, fast_events):+.2e} with no gas "
          f"in the tail")
    check("neither statistic exists without a detachment",
          np.isnan(tail_excess(times, chemistry, [])))


def test_bubble_drop_sigma_enrichment():
    """
    Why BUBBLE_DROP_SIGMA is 6, not 8: the archive-wide sweep behind the
    2026-09-07 change, kept as a check rather than a one-off calculation so it
    cannot silently stop being true of the data.

    8 was a hard cut through a smooth tail: sweeping every step in the block
    that falls short of the OLD cutoff and would survive `_is_excursion`
    unchanged, the count thins gradually from 5 sigma to 8 with no gap in it
    anywhere -- not the signature of a clean threshold. What separates real
    gas from noise here is not the fall's own size, since gas and noise share
    a size distribution in this band; it is whether the CURVE it is on already
    carries a confirmed (>=8 sigma) detachment. Noise would not know that;
    the artefact would.
    """
    print("\nwhy BUBBLE_DROP_SIGMA is 6, not 8")
    OLD_CUTOFF = 8.0
    curve_list = scope.curves(scope.TWO_AXIS_BLOCK)
    bubbling, clean = set(), set()
    per_curve_sigmas = {}
    for curve in curve_list:
        values = np.asarray(curve.absorbance, dtype=float)
        noise = curve.noise
        key = (curve.experiment, curve.sample)
        # floor=0: this sweep is about the SIGMA cutoff alone, from before
        # DETACHMENT_SNR_FLOOR existed -- `test_the_detachment_snr_floor` is
        # where exp 150 cuvette 1's own exclusion is checked.
        events = detachments(values, noise, sigma=OLD_CUTOFF, floor=0)
        (bubbling if events else clean).add(key)
        covered = set()
        for start, stop in events:
            covered.update(range(start, stop))
        sigmas = []
        for index, step in enumerate(np.diff(values)):
            sigma = -step / noise
            if sigma < OLD_CUTOFF and index not in covered and not _is_excursion(
                    values, (index, index + 1)):
                sigmas.append(sigma)
        per_curve_sigmas[key] = sigmas

    def rate_at(cutoff, pool):
        hit = sum(1 for key in pool
                  if any(s >= cutoff for s in per_curve_sigmas[key]))
        return hit / len(pool)

    check("every curve in the block falls in one bin or the other",
          len(bubbling) + len(clean) == len(curve_list),
          f"{len(bubbling)} + {len(clean)} against {len(curve_list)}")

    at_new_cutoff = rate_at(BUBBLE_DROP_SIGMA, bubbling)
    off_new_cutoff = rate_at(BUBBLE_DROP_SIGMA, clean)
    check("curves that already bubble carry a near-threshold fall far more "
          "often than curves that never do",
          at_new_cutoff > 8 * off_new_cutoff,
          f"{at_new_cutoff:.2f} against {off_new_cutoff:.2f}")

    # 6.0 is the lowest cutoff this split supports without picking up a
    # second curve outside the confirmed-bubbling set -- exp 149 cuvette 1's
    # 7.97 sigma fall is the sole exception at every cutoff from 7.5 down to
    # 6.0, and by eye it is a sustained level drop held for four readings
    # after, not a spike, so it reads as the same near-miss population as
    # 143.3 rather than a false positive.
    touched_clean = {key for key in clean
                     if any(s >= BUBBLE_DROP_SIGMA for s in per_curve_sigmas[key])}
    check("only one curve outside the confirmed-bubbling set is touched at "
          "the new cutoff",
          touched_clean == {(149, 1)}, f"{sorted(touched_clean)}")
    touched_below = {key for key in clean
                     if any(s >= BUBBLE_DROP_SIGMA - 0.5 for s in per_curve_sigmas[key])}
    check("half a sigma lower already picks up more, which is why 6.0 and "
          "not lower",
          len(touched_below) > len(touched_clean),
          f"{len(touched_below)} against {len(touched_clean)}")


def test_the_recovery_depth_extension():
    """
    Why `_is_excursion` reaches past the one adjacent reading: the sweep
    behind the 2026-09-07 depth extension, kept as a check for the same
    reason `test_bubble_drop_sigma_enrichment` is -- so the archive-wide
    claim cannot silently stop being true.

    Exp 150 cuvette 1 and exp 151 cuvette 6 are the block's two weakest,
    most drift-dominated curves (net signal 14-21x noise, against 27-143x
    for their sibling cuvettes). Both carry falls that recover only 13-48%
    of themselves in the single adjacent reading the old test looked at, and
    the rest of the way one or two readings later -- a shape the old test
    could not see, and one no real detachment in the archive shares.
    """
    print("\nwhy the recovery test reaches past one reading")
    archive_curves = {(c.experiment, c.sample): c
                      for c in scope.curves(scope.archive())}

    def get(experiment, sample):
        return archive_curves[(experiment, sample)]

    # The two curves the extension was built for. Isolated from
    # DETACHMENT_SNR_FLOOR (floor=0) so this checks the depth extension on
    # its own -- the floor's own effect on exp 150 cuvette 1's remaining
    # four is `test_the_detachment_snr_floor`'s.
    weak_before, weak_after = {}, {}
    for experiment, sample in ((151, 6), (150, 1)):
        curve = get(experiment, sample)
        weak_before[(experiment, sample)] = detachments(
            curve.absorbance, curve.noise, sigma=BUBBLE_DROP_SIGMA, floor=0)
    check("exp 151 cuvette 6 no longer carries any detachment",
          weak_before[(151, 6)] == [], f"{weak_before[(151, 6)]}")
    check("exp 150 cuvette 1 keeps four of its eight",
          len(weak_before[(150, 1)]) == 4, f"{weak_before[(150, 1)]}")

    # Every real detachment and confirmed excursion pinned elsewhere in the
    # package is unmoved by the extension -- this is the regression guard.
    pinned_real = {
        (143, 3): 3, (149, 1): 1, (135, 2): 15, (144, 2): 4, (140, 4): 7,
        (135, 1): 19, (139, 2): 3, (130, 2): 6,
    }
    for (experiment, sample), count in pinned_real.items():
        curve = get(experiment, sample)
        events = detachments(curve.absorbance, curve.noise,
                             sigma=BUBBLE_DROP_SIGMA, floor=0)
        check(f"exp {experiment} cuvette {sample} keeps its {count} "
              f"real detachments",
              len(events) == count, f"{len(events)}: {events}")
    excursion_curve = get(149, 5)
    check("exp 149 cuvette 5 still has none",
          detachments(excursion_curve.absorbance, excursion_curve.noise,
                     sigma=BUBBLE_DROP_SIGMA) == [], "")

    # The extension is one-directional: reaching backward from `start` the
    # same way conflates genuine pre-fall acceleration with a spike, and
    # exp 135 cuvette 1's largest detachment (41.3 sigma) is the case that
    # would be lost. It sits right after four readings of real, fast rise.
    curve = get(135, 1)
    event = (222, 223)
    values = np.asarray(curve.absorbance, dtype=float)
    check("the pre-fall rise into exp 135 cuvette 1's largest detachment is "
          "real acceleration, not noise",
          float(values[222] - values[218]) > 20 * curve.noise,
          f"{(values[222] - values[218]) / curve.noise:.1f} sigma of climb "
          f"over the four readings before it")
    check("and that detachment is not read as an excursion",
          not _is_excursion(values, event), "")

    # The ceiling is what keeps a genuine acceleration right after a fall
    # from being read as the fall's own recovery. Exp 135 cuvette 1's
    # 6.2 sigma detachment at (272, 273) is real and sits right before the
    # curve accelerates hard; uncapped, that acceleration alone would cross
    # the anomaly threshold within the extended window.
    event = (272, 273)
    check("this detachment survives only because recovery is capped",
          not _is_excursion(values, event), "")
    check("  uncapped, the same reach would have rejected it",
          _is_excursion(values, event, ceiling=1e9), "")


def test_the_detachment_snr_floor():
    """
    Why `DETACHMENT_SNR_FLOOR` exists: exp 150 cuvette 1's remaining four
    "detachments" survived every per-event test tried against them --
    recovery depth, capped or not, local noise computed with every other
    candidate fall excluded -- because a curve this weak (net/noise 20.7,
    barely over `live`'s own 20) puts real gas and noise excursions at the
    same size relative to what its own noise estimate can resolve. No
    per-event statistic can tell them apart there; the curve itself has to
    be excluded from the presence question.

    Kept as a check, not a one-off calculation, so the archive-wide gap the
    floor sits in cannot silently close.
    """
    print("\nwhy DETACHMENT_SNR_FLOOR is 30, not fit to one curve")
    archive_curves = {(c.experiment, c.sample): c
                      for c in scope.curves(scope.archive())}

    def get(experiment, sample):
        return archive_curves[(experiment, sample)]

    # The curve the floor was built for loses all four of its remaining
    # candidates; the real detachments already validated against every other
    # test in this file are untouched.
    weak = get(150, 1)
    check("exp 150 cuvette 1 carries none",
          detachments(weak.absorbance, weak.noise, sigma=BUBBLE_DROP_SIGMA)
          == [], "")
    weak_snr = float(weak.absorbance[-1] - weak.absorbance[0]) / weak.noise
    check("  and it sits below the floor",
          weak_snr < DETACHMENT_SNR_FLOOR, f"{weak_snr:.1f}")

    pinned_real = {
        (143, 3): 3, (149, 1): 1, (135, 2): 15, (144, 2): 4, (140, 4): 7,
        (135, 1): 19, (139, 2): 3, (130, 2): 6,
    }
    for (experiment, sample), count in pinned_real.items():
        curve = get(experiment, sample)
        events = detachments(curve.absorbance, curve.noise,
                             sigma=BUBBLE_DROP_SIGMA)
        check(f"exp {experiment} cuvette {sample} still keeps its "
              f"{count} real detachments", len(events) == count,
              f"{len(events)}: {events}")

    # THE GAP THE FLOOR SITS IN. Two curves that are genuinely heavy, dense
    # bubblers -- exp 131 cuvettes 1 and 2, whose own bubble_load (6.5 and
    # 8.3) is as high as exp 150 cuvette 1's (5.4), so load alone cannot
    # separate them -- sit at net/noise 36.8-44.6 and keep every one of
    # their 18 and 19 detachments. Nothing in the archive sits between the
    # weak curve's 20.7 and the heavy bubblers' 36.8.
    for experiment, sample, count in ((131, 1, 18), (131, 2, 19)):
        curve = get(experiment, sample)
        events = detachments(curve.absorbance, curve.noise,
                             sigma=BUBBLE_DROP_SIGMA)
        check(f"exp {experiment} cuvette {sample}, a genuine heavy "
              f"bubbler, keeps all {count}",
              len(events) == count, f"{len(events)}")
        snr = float(curve.absorbance[-1] - curve.absorbance[0]) / curve.noise
        check("  well above the floor", snr > DETACHMENT_SNR_FLOOR,
              f"{snr:.1f}")

    # THE GAP ITSELF, over every curve in the archive that carries even a
    # CANDIDATE fall (bubble_drops, before the excursion test or the floor
    # are applied) -- the widest population the floor could possibly matter
    # to, not just the ones that end up confirmed.
    candidate_snrs = []
    for curve in archive_curves.values():
        if not len(bubble_drops(curve.absorbance, curve.noise,
                                sigma=BUBBLE_DROP_SIGMA)):
            continue
        net = float(curve.absorbance[-1] - curve.absorbance[0])
        if curve.noise > 0 and np.isfinite(net):
            candidate_snrs.append(net / curve.noise)
    candidate_snrs.sort()
    below = [s for s in candidate_snrs if s < DETACHMENT_SNR_FLOOR]
    above = [s for s in candidate_snrs if s >= DETACHMENT_SNR_FLOOR]
    check("the floor sits in a real gap, not a hand-picked line",
          bool(below) and bool(above) and max(below) < DETACHMENT_SNR_FLOOR
          <= min(above),
          f"{max(below):.1f} < {DETACHMENT_SNR_FLOOR} <= {min(above):.1f}")


def test_bubble_gains():
    """
    Gas arriving in the beam, not leaving it -- the rare mirror of a
    detachment, and why it cannot be found the way one is.

    A fall past `BUBBLE_DROP_SIGMA` needs no further test to be suspect: real
    chemistry never falls. A rise past the same threshold is not suspect on
    its own -- most large rises in the two-axis block are the reaction, 809
    against 303 falls -- so `bubble_gains` needs a rise to pass two tests a
    fall does not: it must not reverse (recovery, reused from `_is_excursion`
    on the negated curve) and it must be a KINK against the curve's own local
    trend (`local_outlier_z`), never merged across readings the way a fall
    is, because a genuine multi-reading acceleration would merge into one
    giant false jump if it were.
    """
    print("\ngas arriving, not leaving")
    times = np.arange(0, 3600, 60.0)
    noise = 1e-4
    # An ordinary step here is ~0.6 sigma -- comfortably under every
    # threshold, so a planted event is unambiguous against the background.
    chemistry = 0.02 + 1e-6 * times

    # A SPIKE UP THAT REVERTS: one reading jumps and the next is back on the
    # line. The recovery test must reject it -- this is real chemistry (or
    # noise), not gas that arrived and stayed.
    perched = chemistry.copy()
    perched[30] += 0.004
    check("a spike that reverts is not a gain",
          bubble_gains(times, perched, noise) == [],
          f"{bubble_gains(times, perched, noise)}")

    # A PERSISTENT STEP: the level jumps and stays. This is what a gain is
    # for, and its size should read off almost exactly, net of the ordinary
    # step the curve was already taking there.
    step = chemistry.copy()
    step[30:] += 0.004
    found = bubble_gains(times, step, noise)
    check("a persistent step is a gain",
          len(found) == 1 and found[0][0] == 30, f"{found}")
    check("  and its size is the jump, not the ordinary step under it",
          found and abs(found[0][1] - 0.004) < 1e-9,
          f"{found[0][1]:.6f}" if found else "none")

    rng = np.random.default_rng(0)
    noisy = step + rng.normal(0, noise, len(times))
    noisy_found = bubble_gains(times, noisy, noise)
    check("the same step survives realistic noise",
          len(noisy_found) == 1 and noisy_found[0][0] == 30
          and abs(noisy_found[0][1] - 0.004) < 5 * noise,
          f"{noisy_found}")

    # A SMOOTH, FAST ACCELERATION: real kinetics can rise by many sigma a
    # step for many consecutive readings. NOT ONE of those steps may score as
    # a gain -- this is exactly the failure mode `bubble_gains` exists to
    # avoid, and it is why events are never merged across readings the way a
    # fall's are.
    fast = 0.02 + 0.06 * (1 - np.exp(-times / 300.0))
    check("a smooth acceleration has no gain, at any step",
          bubble_gains(times, fast, noise) == [],
          f"{bubble_gains(times, fast, noise)}")

    check("a clean curve has no gain",
          bubble_gains(times, chemistry, noise) == [])

    # apply_gains is the level shift alone, checked independent of detection.
    shifted = apply_gains(step, [(30, 0.004)])
    check("apply_gains lowers everything from its index on, and nothing "
          "before it",
          np.allclose(shifted[:30], step[:30])
          and np.allclose(shifted[30:], step[30:] - 0.004))

    # THE REAL CASE THE KINK TEST WAS BUILT FOR: exp 135 cuvette 5's jump at
    # 9780 s, and exp 146 cuvette 4's jump near the end of its run -- a
    # bubble that arrived and never left before the recording stopped, with
    # no detachment anywhere on the curve.
    archive_curves = {(c.experiment, c.sample): c
                      for c in scope.curves(scope.archive())}
    jump = archive_curves[(135, 5)]
    found = bubble_gains(np.asarray(jump.times, dtype=float),
                         np.asarray(jump.absorbance, dtype=float), jump.noise)
    check("exp 135 cuvette 5's jump at 9780 s is a gain",
          len(found) == 1 and jump.times[found[0][0]] == 9780.0, f"{found}")

    holding = archive_curves[(146, 4)]
    check("exp 146 cuvette 4 carries no detachment at all",
          detachments(holding.absorbance, holding.noise) == [])
    found = bubble_gains(np.asarray(holding.times, dtype=float),
                         np.asarray(holding.absorbance, dtype=float),
                         holding.noise)
    check("but it does carry a gain, never watched to leave",
          len(found) == 1, f"{found}")

    # THE NEGATIVE CASE THE MERGE BUG PRODUCED: exp 144 cuvette 2 climbs
    # 20-30 sigma a step for readings 29-42, real and smooth. Merged into one
    # span the way a fall's consecutive candidates are, this scored as a
    # single ~0.034 AU jump -- larger than any real gain found anywhere else
    # in the block. Unmerged, no gain may fall inside that stretch.
    fourteen = archive_curves[(144, 2)]
    found = bubble_gains(np.asarray(fourteen.times, dtype=float),
                         np.asarray(fourteen.absorbance, dtype=float),
                         fourteen.noise)
    check("exp 144 cuvette 2's real 14-reading acceleration is not a gain",
          not any(29 <= index <= 43 for index, _ in found), f"{found}")
    check("  and nothing on that curve is anywhere near that size",
          all(gain < 0.01 for _, gain in found), f"{found}")

    # THE SAME CURVE-LEVEL GATE AS `detachments`. exp 150 cuvette 1 sits
    # below DETACHMENT_SNR_FLOOR, and the gate excludes it here for the
    # identical reason -- a curve this weak cannot license a per-event call,
    # rise or fall alike.
    weak = archive_curves[(150, 1)]
    check("exp 150 cuvette 1 carries no gain either, gated by the same "
          "floor",
          bubble_gains(np.asarray(weak.times, dtype=float),
                      np.asarray(weak.absorbance, dtype=float),
                      weak.noise) == [])


def test_debubble_with_gains():
    """
    `debubble` folds `bubble_gains` on top of the falls model, and neither
    guarantee the falls model already had may be weaker for it.

    `worst_at_event` (every detachment corrected in full) and `gas_at_end`
    being exactly the falls-component's own zero -- not the readings' or a
    gain's -- are both PROVEN by construction of `unreleased_gas`, not merely
    observed; this checks they still hold with gains folded in, over the
    whole two-axis block and not just the curves that carry one.
    """
    print("\ndebubble with gains folded in")
    worst_at_event = []
    rebuilt_worst = np.inf
    for curve in scope.curves(scope.TWO_AXIS_BLOCK):
        times = np.asarray(curve.times, dtype=float)
        values = np.asarray(curve.absorbance, dtype=float)
        rebuilt, events = debubble(times, values, curve.noise)
        for start, stop in events:
            worst_at_event.append(
                float(rebuilt[stop] - rebuilt[start]) / curve.noise)
        if len(rebuilt) > 1:
            rebuilt_worst = min(rebuilt_worst,
                                float(np.diff(rebuilt).min()) / curve.noise)
        gains = bubble_gains(times, values, curve.noise)
        gain_total = sum(gain for _, gain in gains)
        held_at_end = float(values[-1] - rebuilt[-1])
        check(f"exp {curve.experiment} cuvette {curve.sample}: gas held at "
              "the end is exactly its own gains, nothing from the falls "
              "model",
              abs(held_at_end - gain_total) < 1e-9,
              f"{held_at_end:.6f} against {gain_total:.6f}")

    # THE ONE DOCUMENTED EXCEPTION, unmoved: exp 135 cuvette 6's fall is in
    # the first interval, and `bubble_rate` returns `inf` for it, so
    # `debubble` leaves that curve untouched. Gains do not touch it either --
    # there is no gain on that curve -- so the exception is exactly what it
    # was before this function existed.
    check("every detachment is corrected in full, except the one first-"
          "interval case",
          min(worst_at_event) > -9.62 and min(worst_at_event) < -9.60,
          f"{min(worst_at_event):.4f}")

    # THE MERGE BUG'S SIGNATURE: it turned exp 144 cuvette 2's real
    # acceleration into a false gain that landed in the MIDDLE of a real
    # detachment's span, producing a step far worse than any real excursion
    # in the block. Unmerged, the worst step in any reconstruction is back to
    # the excursion `_is_excursion` deliberately leaves alone.
    check("no gain corrupts a real detachment's own step",
          rebuilt_worst > -61.2 and rebuilt_worst < -61.0,
          f"{rebuilt_worst:.4f}")


def test_early_trough():
    """
    `early_trough` finds a sustained decline and rejects a single deep
    outlier that a smoothed mean alone cannot tell apart from one, via a
    LEAVE-ONE-OUT test: drop the trough window's single worst reading and
    require the rest to still average below threshold.

    The real-curve case is calibrated to the one that motivated the design:
    an earlier version required a fixed COUNT of the window's readings to
    individually clear a per-reading depth bar, and it wrongly rejected
    exp 4.1 -- a real, sustained decline (9 of 9 readings in its trough
    window negative) that is simply shallower per reading than exp 4.2's
    more dramatic one. Leave-one-out passes both: removing the single worst
    reading barely moves either curve's mean.
    """
    print("\nthe early trough: sustained decline against a lone outlier")
    times = np.arange(0, 3600, 60.0)
    noise = 2e-4

    # A SUSTAINED DECLINE: every reading in the window is genuinely low, not
    # just the smoothed mean.
    declining = np.zeros_like(times)
    declining[:15] = -0.0015 * (1 - np.exp(-times[:15] / 600.0))
    declining[15:] = declining[14] + 3e-6 * (times[15:] - times[14])
    z, t_trough, start, sustained = early_trough(times, declining, noise)
    check("a sustained decline is found", z < -4.0, f"z={z:.2f}")
    check("and marked sustained", sustained, f"z={z:.2f}")
    check("start is within the declining region", 0 <= start < 20,
          f"start={start}")

    # ONE DEEP OUTLIER, otherwise flat: the smoothed window it sits in can
    # still average out negative, but removing that single reading collapses
    # the mean back toward zero -- the leave-one-out test's whole point.
    rng = np.random.default_rng(0)
    spiky = rng.normal(0.0, noise * 0.3, size=len(times))
    spiky[10] -= 0.006
    z_spike, _, _, sustained_spike = early_trough(times, spiky, noise)
    check("a lone spike can still smooth to a deep z",
          z_spike < -4.0, f"z={z_spike:.2f}")
    check("but is not sustained (leave-one-out collapses it)",
          not sustained_spike, f"z={z_spike:.2f}")

    # CLEAN NOISE: no decline at all.
    clean = rng.normal(0.0, noise * 0.5, size=len(times))
    z_clean, _, _, sustained_clean = early_trough(times, clean, noise)
    check("clean noise is not sustained", not sustained_clean,
          f"z={z_clean:.2f}")

    # TOO SHORT: fewer readings than one window needs.
    z_short, t_short, start_short, sustained_short = early_trough(
        times[:5], declining[:5], noise)
    check("too few readings returns nan, not a spurious trough",
          np.isnan(z_short) and start_short == -1 and not sustained_short)

    print("\nthe early trough against two real curves (exp 4.2, exp 4.1)")
    lookup = {(c.experiment, c.sample): c for c in scope.curves((4,))}
    real_dip = lookup[(4, 2)]
    z_real, _, _, sustained_real = early_trough(
        np.asarray(real_dip.times, dtype=float),
        np.asarray(real_dip.absorbance, dtype=float), real_dip.noise)
    check("exp 4.2's real, sustained dip is found",
          z_real < -30.0 and sustained_real, f"z={z_real:.2f}")
    sibling = lookup[(4, 1)]
    z_sibling, _, _, sustained_sibling = early_trough(
        np.asarray(sibling.times, dtype=float),
        np.asarray(sibling.absorbance, dtype=float), sibling.noise)
    check("exp 4.1's shallower, genuine dip clears -4 sigma",
          z_sibling < -4.0, f"z={z_sibling:.2f}")
    check("and IS sustained under leave-one-out (the corrected verdict)",
          sustained_sibling, f"z={z_sibling:.2f}")


if __name__ == "__main__":
    test_every_test_is_actually_run()
    test_the_runner_finds_every_gate()
    test_no_duplicate_definitions()
    test_the_duplicate_guard_catches_a_planted_duplicate()
    test_lag_statistic()
    test_noise_and_rate()
    test_acceleration()
    test_floor_belongs_to_the_source()
    test_whole_curve_estimators_return_rates()
    test_outlier_flagging()
    test_model_residual()
    test_segmented_fit()
    test_two_breakpoints()
    test_the_bubble_correction()
    test_the_excursion_test_on_planted_spikes()
    test_the_monotone_bound()
    test_the_bubble_the_run_never_shed()
    test_bubble_drop_sigma_enrichment()
    test_the_recovery_depth_extension()
    test_the_detachment_snr_floor()
    test_bubble_gains()
    test_debubble_with_gains()
    test_early_trough()
    print(f"\n{len(FAILURES)} failure(s)" + (": " + ", ".join(FAILURES) if FAILURES else ""))
    sys.exit(1 if FAILURES else 0)
