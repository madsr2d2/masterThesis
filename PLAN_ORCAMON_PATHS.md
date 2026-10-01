# orcamon — reaction paths: IRC and NEB (PA1–PB2)

**Written:** 2026-10-01. **Verified against:** `master` at `ea580b4` (the defect-fix run, `PLAN_ORCAMON_DEFECTS.md`, is committed; every regex below was prototyped against the three real IRC/NEB outputs and the synthetic fixtures in § Facts, and every expected value was computed by hand from them).
**Branch:** `master` (the agents commit straight to `master`).
**Progress file:** `PLAN_ORCAMON_PATHS_PROGRESS.md` (template at the end).

orcamon reads an IRC or an NEB as if it were a single point. `orcamon show` on a finished 101-point IRC prints one energy, `energies` prints one row, and the TUI shows the transition state and nothing else. Yet ORCA writes every point's energy into the `.out` as it goes, and every point's geometry into trajectory files that the `.out` names. This plan reads both. It adds a signed point numbering, a one-line path summary in `show` and `ls`, an `irc_not_converged` flag, `energies` for paths, `geom`/`snapshot --point N`, and a TUI chart that scrubs through the path, so a person can watch the reaction happen. **The goal a deviation must still meet:** every point on an IRC or NEB path can be listed and drawn, each with the energy and geometry ORCA wrote for it; nothing about a job that is not an IRC or an NEB changes; and the JSON key sets do not change.

## Execution

The opencode orchestrator runs this plan. These are this project's rules; the loop itself belongs to the agents and is not repeated here.

### Verification commands

| What | Working directory | Command | Expected now | Time |
|---|---|---|---|---|
| renderer and parser tests | `/home/madsr2d2/masterThesis` | `.venv/bin/python test_monitor.py` | 161 `pass` lines, last line `0 failure(s)` | ~2 s |
| CLI, TUI and skill tests | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_paths/xdg .venv/bin/python test_orcamon.py` | 178 `pass` lines, last line `0 failure(s)` | ~10 s |
| duplicate-definition guard | `/home/madsr2d2/masterThesis` | `.venv/bin/python data/test_curve_metrics.py` | last line `0 failure(s)` | ~45 s |
| every gate | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_paths/xdg .venv/bin/python run_gates.py` | last block says `38 gates in … on 8 jobs, 0 failed` | 500–750 s; give the bash call a timeout of at least 900000 ms |

Count `pass` lines with `.venv/bin/python test_monitor.py | grep -c '^  pass'` (same for `test_orcamon.py`, with its environment variables). Each task states the counts it must end at.

Every task runs the first three commands. Only PB2 runs `run_gates.py`: every task touches only `tools/orcamon/`, `test_monitor.py`, `test_orcamon.py` and `.claude/skills/orcamon/SKILL.md`, which no other gate reads except the duplicate guard.

**Do not add a new root `test_*.py` file.** `test_root_documents.py` pins the gate count to `CLAUDE.md`, which is the user's; a new gate file would change the count and break that gate. New tests go in `test_monitor.py` or `test_orcamon.py`.

`test_monitor.py` and `test_orcamon.py` run their tests from an explicit list in `if __name__ == "__main__":` at the bottom of each file. A new test function that is not added to that list never runs.

### Stop conditions

Stop and report NEEDS_USER (implementer) or BLOCK (reviewer) when:
- a hand-derived expected value in this plan does not hold after the change is implemented as written (report the measured value; do not change the expected value, the regex or the design to make it pass);
- a regex written in this plan does not match the line of a real output quoted in § Facts (report the line and the regex);
- any existing test fails and the fix would mean changing that test's assertion, expected value or tolerance;
- a real-tree acceptance value differs from this plan's (report it; the tree under `computational/` is the user's and may have changed — say whether the job's files changed since `ea580b4`);
- in PB2, a gate that does not involve orcamon fails (report it; do not fix it);
- the working tree has changes not made by this plan's tasks.

A different name or path for the same thing is not a stop: adapt, and report it as a deviation.

### Never

- Never edit `CLAUDE.md`, the root `README.md` or `.claude/skills/run-orca/SKILL.md`; they are the user's. (`tools/orcamon/README.md`, `tools/orcamon/src/orcamon/SKILL.md` and the generated `.claude/skills/orcamon/SKILL.md` are orcamon's and may change where a task says so; the generated copy changes ONLY through `orcamon skill install --project .`.)
- Never `git add -A`, `git add .` or `git commit -a`; stage the files the task lists, by path.
- Never push.
- Never add a new root `test_*.py` gate file.
- Never run `orcamon` against the user's cache: every manual `orcamon` command sets `XDG_CACHE_HOME=/tmp/orcamon_paths/xdg`. Anything that starts the TUI sets `HERDR_ENV=0`.
- Never write into `computational/`, never start ORCA, and never kill or signal a process the task did not start. Tests use the synthetic fixtures in § Facts, written to temporary directories; the real tree appears only in read-only acceptance commands.
- Never change an existing test's assertion, expected value or tolerance.
- Never add a dependency. `core` and `cli` import only the standard library at module level.
- Never put thesis paths, job names or chemistry names into orcamon's package code or its shipped skill.
- Never add or remove a key in any `--json` document, and never change `REPORT_SCHEMA` (it stays `1`). New path data reaches JSON only through the `rows` of `energies --json` (D7).
- Never add a fenced code block containing an `orcamon …` line to `tools/orcamon/src/orcamon/SKILL.md`: `test_orcamon.test_the_skill_examples_run` executes every such line and would add checks.
- Never parse an ORCA marker this plan does not list (D8).

## Scope

**In scope:** parsing the IRC iteration rows and the NEB iteration rows; a multi-frame XYZ reader; `core/paths.py` (path points, the path view, the summary line); the `irc_not_converged` flag; the `show`/`ls`/`wait` text; `energies` for paths; `geom --point` and `snapshot --point`; the TUI chart for paths; the tests, README and shipped skill for all of it.
**Out of scope:** see § Out of scope.

## Decisions

| # | Decision | By, date |
|---|---|---|
| D1 | This plan is IRC and NEB only. Keeping real normal modes, the Act layer and the Steer layer are later plans. | user, 2026-10-01 |
| D2 | Live progress and every point's energy come from the iteration ROWS ORCA prints in the `.out` (IRC: one block per direction; NEB: one block per phase), read as they are written. Geometries come from the trajectory files the `.out` NAMES (`Storing … trajectory in …`, `Current trajectory will be written to …`), never from guessed file names. The `IRC PATH SUMMARY` block is not parsed: it is the same rows, reordered (verified on both real IRCs). | planner, 2026-10-01 |
| D3 | IRC points are numbered by a signed index: backward iteration `k` is `-(k+1)`, the TS is `0`, forward iteration `k` is `+(k+1)`. Labels keep ORCA's own iteration numbers: `IRC backward k`, `IRC TS`, `IRC forward k`. NEB points are the image numbers `0 … n-1`, labelled `NEB image i`, with ` (CI)` or ` (HEI)` on the image the latest row names. | planner, 2026-10-01 |
| D4 | IRC dE is ORCA's own `dE(kcal/mol)` column times `KCAL_TO_KJ = 4.184`: relative to the TS, needs no TS energy, and is what ORCA printed. The TS point's dE is `0.0` and its energy is the job's first printed combined energy. NEB dE is `(E_i - E_0) * EH_TO_KJ_PER_MOL`, from the energies in the MEP trajectory's comment lines. | planner, 2026-10-01 |
| D5 | IRC geometries: from the full trajectory when its frame count is `len(backward) + 1 + len(forward)` (path order); otherwise forward iteration `k` is frame `k` of the forward file, backward `k` frame `k` of the backward file, and the TS is the input's coordinate file when its atom count matches the trajectories'. A point whose frame is not written yet has no atoms and is said to be "not written yet", never borrowed from another point. | planner, 2026-10-01 |
| D6 | The point the chart rings when it is following (the "focus"): for an IRC, the TS once the job has ended (normal termination or an error line), otherwise the newest point (the backward end if any backward row exists, since ORCA runs forward first; otherwise the forward end). For an NEB, the image the latest row names. | planner, 2026-10-01 |
| D7 | No JSON key is added. `show --json`/`ls --json` are unchanged (path progress is text-only, through private `JobReport` fields as `_converged_reason` is). Path data reaches agents through `energies --json`, whose top-level keys stay `job, mode, truncated, rows`, with `mode` `irc` or `neb`. `irc_not_converged` is a new flag CODE (a value, not a key) and is added to `FLAG_CODES`, the README and the shipped skill. | planner, 2026-10-01 |
| D8 | Only markers with a real sample are parsed. There is none for an IRC that converged in a direction, for a converged NEB, for NEB-TS's closing TS optimization, or for `*_NEB-TS_converged.xyz` — so none of those is read. An IRC direction is reported as "hit MaxIter" only on ORCA's own MaxIter banner. | planner, 2026-10-01 |
| D9 | Commits go straight to `master`, one per task. Tests go in `test_monitor.py` (parser and core) and `test_orcamon.py` (CLI and TUI). | user, 2026-10-01 |

## Facts (verified 2026-10-01)

**Environment.** `/home/madsr2d2/masterThesis/.venv/bin/python` is Python 3.12.3 with numpy, PIL and textual; orcamon is installed editable, console script `/home/madsr2d2/masterThesis/.venv/bin/orcamon`. `git status --porcelain` was empty at `ea580b4`.

**The real jobs** (all under `computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/`; read-only):
- `geometry_b973c-xtb/irc_bridge/` — finished IRC (`! QM/XTB B97-3c ddCOSMO(Water) IRC TightSCF`, `%irc MaxIter 50 Direction both Monitor_Internals …`, `* xyzfile 0 1 job.xyz`, 137 atoms). `job.out`: line 2029–2031 `Storing full IRC trajectory in      .... job_IRC_Full_trj.xyz`, `Storing forward trajectory in       .... job_IRC_F_trj.xyz`, `Storing backward trajectory in      .... job_IRC_B_trj.xyz`; `Max. no of cycles        MaxIter    .... 50`; forward block from line 3018, backward from 3080, each ending with the banner `*  MAXIMUM NUMBER OF ITERATIONS REACHED - STOPPING IRC RUN  *`. Header: `Iteration    E(Eh)      dE(kcal/mol)  max(|G|)   RMS(G)  B(O 133,C 128) B(O 129,C 128) B(O 133,H 132)`, then `Convergence thresholds                0.002000  0.000500`, then 50 rows per direction such as `    0     -1011.280182    5.201526    0.074732  0.005858      1.66         1.29         0.92   ` (trailing spaces). Forward: 50 rows, last `49 -1011.290477 dE -1.258850`. Backward: 50 rows, first `0 -1011.298405 dE -6.233759`, last `49 -1011.386396 dE -61.448643`. The job's one combined energy before the IRC is `FINAL SINGLE POINT ENERGY (QM/QM2)    -1011.288471006717` (the TS; the summary block's `<= TS` row has the same energy). `job_IRC_Full_trj.xyz` has 101 frames, `job_IRC_F_trj.xyz` 50, `job_IRC_B_trj.xyz` 50, each frame 137 atoms with comment `Coordinates from ORCA-job job_IRC_Full E -1011.386395866175`. Full frame order = backward 49 … backward 0, TS, forward 0 … forward 49 (frame 1 energy = backward 49's, frame 51 = TS, frame 52 = forward 0's). F frame `k` = forward iteration `k`; B frame `k` = backward iteration `k`. The `.out` prints only a 23-atom coordinate block, so `orcamon show` currently says `(21 of 23 atoms)`.
- `geometry_r2scan3c-xtb/irc/` — finished IRC, `MaxIter .... 20`, both directions hit MaxIter; 20 rows each way; forward last dE `+1.413692`, backward last dE `-17.416051`; monitors `B(O 135,C 128) B(H 130,O 129) B(O 131,H 130)`; Full trajectory 41 frames of 140 atoms.
- `geometry_b973c-xtb/neb_anion/` — a STOPPED `NEB-TS` (no termination line). `job.out`: `Current trajectory will be written to    ....  job_MEP_trj.xyz`; a header `Optim.  Iteration  HEI  E(HEI)-E(0)  max(|Fp|)   RMS(Fp)    dS`, `Switch-on CI threshold               0.020000 `, 12 rows `   LBFGS     0      4    0.129304    0.067073   0.005816  24.3765       ` … iteration 11; a blank line; `Image  4 will be converted to a climbing image in the next iteration (max(|Fp|) < 0.0200) `; a blank line; header `Optim.  Iteration  CI   E(CI)-E(0)   max(|Fp|)   RMS(Fp)    dS     max(|FCI|)   RMS(FCI)`, `Convergence thresholds               0.020000   0.010000            0.002000    0.001000 `, rows 12–14 with two extra columns, the last `   LBFGS    14      4    0.064973    0.020692   0.001725  24.0491    0.017823    0.002398       `, which is the file's last line. `job_MEP_trj.xyz` has 10 frames (2 end points + 8 images) of 133 atoms, comments `Coordinates from ORCA-job job_MEP E -934.526630281048`; frame energies give dE (kJ/mol) `0.00, 13.10, 110.88, 151.54, 170.59, 159.10, 171.67, 173.02, 148.25, 125.51` (image 4 is the CI; it is NOT the highest on this path, which is the data, not a bug).

**Hand-computed real-tree values** (`EH_TO_KJ_PER_MOL = 2625.4996394799`, `KCAL_TO_KJ = 4.184`):
- `irc_bridge`: 101 points, indices -50 … 50; point -50 `IRC backward 49` dE `-61.448643 * 4.184 = -257.1011`; point 0 `IRC TS` energy `-1011.288471006717`, dE 0; point 50 `IRC forward 49` dE `-1.258850 * 4.184 = -5.2670`.
- `geometry_r2scan3c-xtb/irc`: 41 points; backward end `-72.8688`, forward end `+5.9149`.
- `neb_anion`: 10 images; image 4 dE `170.59` (from the frame energies) and the latest row's `E(CI)-E(0)` `0.064973 * 2625.4996 = 170.5866`.

**The regexes, prototyped** (each matched every line it must on the three real outputs and the synthetic fixtures below, and nothing else in them):
```python
_IRC_DIRECTION_RE = re.compile(r"^\s*\*\s+(FORWARD|BACKWARD) IRC\s+\*\s*$")
_IRC_HEADER_RE = re.compile(r"^Iteration\s+E\(Eh\)\s+dE\(kcal/mol\)\s+max\(\|G\|\)\s+RMS\(G\)(.*)$")
_IRC_ROW_RE = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)((?:\s+-?\d+\.\d+)*)\s*$")
_IRC_MONITOR_RE = re.compile(r"[A-Za-z]+\([^)]*\)")
_IRC_MAXITER_RE = re.compile(r"MAXIMUM NUMBER OF ITERATIONS REACHED - STOPPING IRC RUN")
_IRC_TRAJECTORY_RE = re.compile(r"Storing (full|forward|backward) (?:IRC )?trajectory in\s+\.+\s+(\S+)")
_NEB_HEADER_RE = re.compile(r"^Optim\.\s+Iteration\s+(HEI|CI)\s+E\((?:HEI|CI)\)-E\(0\)")
_NEB_ROW_RE = re.compile(r"^\s*([A-Za-z][A-Za-z-]*)\s+(\d+)\s+(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)(?:\s+(\d+\.\d+)\s+(\d+\.\d+))?\s*$")
_NEB_TRAJECTORY_RE = re.compile(r"Current trajectory will be written to\s+\.+\s+(\S+)")
```
`_IRC_HEADER_RE` and `_NEB_HEADER_RE` are matched against the STRIPPED line; the others against the line. The IRC row regex does not match `Convergence thresholds …` or `Max. no of cycles …`; the `IRC PATH SUMMARY` block's header starts `Step`, so `_IRC_HEADER_RE` never opens it.

**The synthetic IRC output** (`_IRC_OUT` in both test files; copy verbatim, leading spaces included). A finished IRC with MaxIter 3: forward 3 rows and the MaxIter banner, backward 2 rows and none:
```
Max. no of cycles        MaxIter    .... 3
Storing full IRC trajectory in      .... job_IRC_Full_trj.xyz
Storing forward trajectory in       .... job_IRC_F_trj.xyz
Storing backward trajectory in      .... job_IRC_B_trj.xyz
FINAL SINGLE POINT ENERGY      -100.000000000000

         *************************************************************
         *                          FORWARD IRC                      *
         *************************************************************

Iteration    E(Eh)      dE(kcal/mol)  max(|G|)   RMS(G)  B(O 0,H 1) B(O 0,H 2)
Convergence thresholds                0.002000  0.000500
    0     -100.001594   -1.000000    0.010000  0.001000      0.97         0.96
    1     -100.003187   -2.000000    0.008000  0.000800      0.98         0.95
    2     -100.004781   -3.000000    0.006000  0.000600      0.99         0.94

         *************************************************************
         *  MAXIMUM NUMBER OF ITERATIONS REACHED - STOPPING IRC RUN  *
         *************************************************************


         *************************************************************
         *                          BACKWARD IRC                     *
         *************************************************************

Iteration    E(Eh)      dE(kcal/mol)  max(|G|)   RMS(G)  B(O 0,H 1) B(O 0,H 2)
Convergence thresholds                0.002000  0.000500
    0     -100.000797   -0.500000    0.010000  0.001000      0.95         0.97
    1     -100.002390   -1.500000    0.008000  0.000800      0.94         0.98

                             ****ORCA TERMINATED NORMALLY****
```
Its input (`job.inp`): `"! B97-3c IRC\n%irc MaxIter 3 end\n* xyzfile 0 1 start.xyz\n"`.

**The synthetic IRC tree** (`_irc_tree(d: Path)` in both test files writes it into directory `d`): `job.inp` and `job.out` as above; and four XYZ files whose 3 atoms are `O (0, 0, z)`, `H (0.96, 0, z)`, `H (-0.24, 0.93, z)`, so a frame is identified by its O's z:
- `start.xyz`: one frame, z = 2.0;
- `job_IRC_Full_trj.xyz`: 6 frames, z = 0, 1, 2, 3, 4, 5;
- `job_IRC_F_trj.xyz`: 3 frames, z = 3, 4, 5 (forward iteration k is z = 3 + k);
- `job_IRC_B_trj.xyz`: 2 frames, z = 1, 0 (backward iteration k is z = 1 - k).
Every frame is `3`, a comment line `Coordinates from ORCA-job job E -100.000000`, then the three atom lines, each formatted `f"{el} {x:.6f} {y:.6f} {z:.6f}"`.

**Hand-derived values for the synthetic IRC:** indices `[-2, -1, 0, 1, 2, 3]`; labels `["IRC backward 1", "IRC backward 0", "IRC TS", "IRC forward 0", "IRC forward 1", "IRC forward 2"]`; dE (kJ/mol, 3 d.p.) `[-6.276, -2.092, 0.0, -4.184, -8.368, -12.552]`; O z per point `[0, 1, 2, 3, 4, 5]` (with or without the full trajectory); TS energy `-100.0`; focus index `2` (finished); last point's monitors `{"B(O 0,H 1)": 0.99, "B(O 0,H 2)": 0.94}`; `irc_maxiter == ["forward"]`; summary `IRC        backward 2 points, end -6.3 kJ/mol · forward 3 points (MaxIter), end -12.6 kJ/mol (from the TS)`; progress `irc B2/F3`; flag message `IRC forward hit MaxIter (3) before a minimum`.

**The synthetic NEB output** (`_NEB_OUT` in both test files; copy verbatim). A running NEB-CI, last block open:
```
Current trajectory will be written to    ....  job_MEP_trj.xyz

Starting iterations:

Optim.  Iteration  HEI  E(HEI)-E(0)  max(|Fp|)   RMS(Fp)    dS
Switch-on CI threshold               0.020000
   LBFGS     0      2    0.010000    0.050000   0.005000  2.0000
   LBFGS     1      2    0.009000    0.030000   0.004000  1.9000

Image  2 will be converted to a climbing image in the next iteration (max(|Fp|) < 0.0200)

Optim.  Iteration  CI   E(CI)-E(0)   max(|Fp|)   RMS(Fp)    dS     max(|FCI|)   RMS(FCI)
Convergence thresholds               0.020000   0.010000            0.002000    0.001000
   LBFGS     2      2    0.008000    0.015000   0.003000  1.8000    0.010000    0.002000
```
Its input: `"! B97-3c NEB-CI\n%neb NEB_End_XYZFile \"product.xyz\" end\n* xyzfile 0 1 start.xyz\n"`. **The synthetic NEB tree** (`_neb_tree(d)`): `job.inp`, `job.out`, `start.xyz` (z = 0) and `job_MEP_trj.xyz` with 4 frames, image `i` at z = i, comments `Coordinates from ORCA-job job_MEP E <energy>` with energies `-100.000000000000`, `-99.995000000000`, `-99.992000000000`, `-99.998000000000`.

**Hand-derived values for the synthetic NEB:** 3 rows; last row `NebRow(iteration=2, phase="CI", image=2, barrier_eh=0.008, max_fp=0.015, rms_fp=0.003)`; phases `["HEI", "HEI", "CI"]`; images `[0, 1, 2, 3]`; labels `["NEB image 0", "NEB image 1", "NEB image 2 (CI)", "NEB image 3"]`; dE (4 d.p.) `[0.0, 13.1275, 21.004, 5.251]` (0.005, 0.008 and 0.002 Eh times 2625.4996394799); focus 2; summary `NEB        iteration 2 · climbing image 2 · E(CI)-E(0) +21.0 kJ/mol`; progress `neb 2 CI2`.

**The code** (`tools/orcamon/src/orcamon/`, at `ea580b4`):
- `core/parser.py`: `_RARE_MARKERS_RE` is the literal prefilter for `_feed_marker`; `_LINE_OF_INTEREST_RE` includes its pattern, so a literal added there reaches both the line path and `feed_text`'s fast path. `JobState._in_block()` returns whether any block is open, and while one is, `feed_text` hands every line to `_process`. `_process` calls `_feed_marker` (behind the prefilter), the scan banner, the QM1 list, SCF rows, the convergence items, then `_feed_blocks`. `reset_for_restart` copies every field from a fresh instance. `cache.save` refuses while `_in_block()` is true. The cache key is the hash of `parser.py`, so editing it invalidates every entry.
- `core/geometry.py`: `read_xyz(path)` reads ONE frame (element letters kept, `.capitalize()`d); `FileGeometry` is the stand-in point for a file geometry (`source`, `scan_step=None`, `cycle=0`, `energy=None`, `energy_label=None`, `key=None`).
- `core/job.py`: `Job.__init__` sets `self.file_geometry` and `self.parsed`; `job.state.path` is the job directory; `job.input.coords_file` the file a `* xyzfile` line names.
- `core/report.py`: `JobReport`'s private fields end with `_liveness_note` and `_converged_reason` (`field(default=None, repr=False)`); `report_lines` emits keyed lines in order …, `"hessian"`, `"imaginary"`, `"energy"`, `"crash"`; `geometry_shown(job, point)` and `describe_point(point)` serve the panes and `geom`.
- `core/status.py`: `FLAG_CODES` lists `"failed", "stopped", "stalled", "quiet", "opt_not_converged", "scan_incomplete", "ts_hessian", "ts_imaginary", "minimum_imaginary", "qm2_errors"`; `attention()` appends the `scan_incomplete` flag right after the `opt_not_converged` ones; `maxiter = f" ({state.max_cycles})" if state.max_cycles else ""` is already computed there.
- `core/units.py`: `EH_TO_KJ_PER_MOL = 2625.4996394799`.
- `cli/commands.py`: `_progress(r)` (the `ls` PROGRESS cell and `wait`'s timeout lines), `cmd_energies`, `_find_point`, `cmd_geom`, `cmd_snapshot`, `UsageError`; `cli/__init__.py::build_parser` declares the flags (`geom`: `--cycle`, `--step`, `--region`; `snapshot`: `--mode`, `--phase`, …).
- `tui/app.py::ConvergencePlot`: `_series(state)` picks the plotted points (modes `scf`, `scan`, `all`, `opt`); `show_job(job)` replots only when `_signature_of(job)` changes; `_place_selection`, `selected_point`, `_select`, `action_select_prev/next/latest`; `MonitorApp._render_detail` passes `chart.selected_point()` to the geometry pane, and `_title_geometry` puts `describe_point(shown)` in its title.
- `SKILL.md` (package) has an "Attention flags" table and a "Looking closer" section; `.claude/skills/orcamon/SKILL.md` is generated from it by `.venv/bin/orcamon skill install --project /home/madsr2d2/masterThesis`, and `test_orcamon.test_this_repository_has_the_current_skill` fails until that is run after any change to the package `SKILL.md` or to a flag.
- `tools/orcamon/README.md` lists the flag codes in "Status and attention flags" and pastes every command's `--help` in "Command reference".

**Tests.** `test_monitor.py`: 161 `pass`; `test_orcamon.py`: 178 `pass`; both `0 failure(s)`. Helpers: `check(label, ok, detail="")` in both; `_feed(state, text)` and `_report_job(state, inp_text, liveness, now=None)` in `test_monitor.py`; `_orcamon(argv) -> (code, stdout, stderr)` in `test_orcamon.py`. `test_monitor.test_every_marker_reaches_its_parser` checks all of `validate.MARKER_CASES` as ONE check. `test_orcamon.test_the_tui_runs_headless` shows the headless pattern: `MonitorApp(root, liveness=…, graphics="text", notify_mode="off")`, `async with app.run_test(size=(180, 50)) as pilot`, wait until every job is parsed, move `app.query_one("#job_table")`'s cursor, `widget.focus()`, `await pilot.press(key)`, `await pilot.pause(0.2)`, read `app.query_one("#geometry").border_title`. The duplicate guard (`data/test_curve_metrics.py`) skips names that start with `_`, so the private fixtures may exist in both test files.

## Conventions

- Match the surrounding code: one-line comments where the reason is not obvious; docstrings say WHY.
- New tests go in the file the task names, use `check(...)` for every assertion, build their trees under `tempfile.TemporaryDirectory()`, and are added to that file's `__main__` list right after the test named in the task.
- The fixtures `_IRC_OUT`, `_NEB_OUT`, `_irc_tree`, `_neb_tree` are written once per test file (module-level, private), the first time a task needs them there, exactly as § Facts gives them.
- Expected values are the hand-derived ones in § Facts, never values copied from a run. Floats are compared after `round(…, 3)` (IRC dE) or `round(…, 4)` (NEB dE), as § Facts writes them.
- No network access.

---

# Phase 1 — IRC

### PA1 — Parse the IRC rows

**Why:** an IRC's `.out` prints every point's energy, ORCA's dE from the TS and the monitored internals, one block per direction, and names its trajectory files; orcamon reads none of it (`energies irc_bridge` prints one row).
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/core/parser.py`, `tools/orcamon/src/orcamon/validate.py`, `test_monitor.py`.

**Do:**
1. In `core/parser.py` add the six `_IRC_…` regexes from § Facts ("The regexes, prototyped") after `_CONVERGED_REASON_RE`, under a comment quoting one real header and row. Append `r"|FORWARD IRC|BACKWARD IRC|STOPPING IRC RUN|trajectory in|dE\(kcal/mol\)"` to `_RARE_MARKERS_RE`.
2. Add, before `JobState`:
   ```python
   @dataclass
   class IrcRow:
       """One IRC iteration as ORCA printed it: its energy, ORCA's dE from the
       TS in kcal/mol, and the values of the monitored internals."""
       iteration: int
       energy: float
       de_kcal: float
       monitors: list
   ```
3. Add `JobState` fields after `opt_converged_reason`: `irc_direction: str | None = None` ("forward"/"backward", from the banner); `irc_rows: dict = field(default_factory=dict)` (direction -> list of `IrcRow`, in printed order); `irc_monitors: list = field(default_factory=list)` (the header's monitor names, e.g. `"B(O 0,H 1)"`); `irc_maxiter: list = field(default_factory=list)` (directions whose MaxIter banner was printed, in order); `irc_files: dict = field(default_factory=dict)` (`"full"`/`"forward"`/`"backward"` -> file name). Private, beside `_pending_eig`: `_in_irc_rows: bool = False`, `_irc_seen_row: bool = False`.
4. In `_feed_marker`: on `_IRC_DIRECTION_RE` set `irc_direction` to the lower-cased group; on `_IRC_HEADER_RE.match(line.strip())` set `irc_monitors = _IRC_MONITOR_RE.findall(group 1)`, `_in_irc_rows = True`, `_irc_seen_row = False`; on `_IRC_MAXITER_RE` append `self.irc_direction or "?"` to `irc_maxiter`; on `_IRC_TRAJECTORY_RE` set `irc_files[group 1] = group 2`.
5. Add `_feed_irc_row(self, line)`: `_IRC_ROW_RE.match(line)` -> append `IrcRow(int(g1), float(g2), float(g3), [float(v) for v in g6.split()])` to `irc_rows.setdefault(self.irc_direction or "?", [])` and set `_irc_seen_row = True`; else, if the stripped line starts with `Convergence thresholds` or no row has been seen yet, do nothing; else close the block (`_in_irc_rows = False`). Call it from `_process` as `if self._in_irc_rows: self._feed_irc_row(line)` placed AFTER the `_feed_marker` call (so the header line itself opens the block and is then skipped as "no row yet"). Add `self._in_irc_rows` to `_in_block()`.
6. In `validate.py` append four `MARKER_CASES` rows:
   `("         *                          FORWARD IRC                      *", lambda s: s.irc_direction == "forward")`,
   `("Iteration    E(Eh)      dE(kcal/mol)  max(|G|)   RMS(G)  B(O 0,H 1)", lambda s: s._in_irc_rows and s.irc_monitors == ["B(O 0,H 1)"])`,
   `("         *  MAXIMUM NUMBER OF ITERATIONS REACHED - STOPPING IRC RUN  *", lambda s: len(s.irc_maxiter) == 1)`,
   `("Storing forward trajectory in       .... job_IRC_F_trj.xyz", lambda s: s.irc_files == {"forward": "job_IRC_F_trj.xyz"})`.
7. Add `_IRC_OUT` (§ Facts) to `test_monitor.py` and `test_irc_rows_are_parsed()`, registered right after `test_every_marker_reaches_its_parser`. Feed `_IRC_OUT` with `_feed` into `JobState(path=Path("/nonexistent"), stem="job")`. Exactly 7 `check` calls:
   - `irc_files == {"full": "job_IRC_Full_trj.xyz", "forward": "job_IRC_F_trj.xyz", "backward": "job_IRC_B_trj.xyz"}`;
   - `irc_monitors == ["B(O 0,H 1)", "B(O 0,H 2)"]`;
   - `len(irc_rows["forward"]) == 3 and len(irc_rows["backward"]) == 2 and irc_rows["forward"][2] == IrcRow(2, -100.004781, -3.0, [0.99, 0.94])`;
   - `irc_rows["backward"][0].monitors == [0.95, 0.97]`;
   - `irc_maxiter == ["forward"]`;
   - a second state fed by `state.feed_text(_IRC_OUT.strip("\n") + "\n")` has equal `irc_rows`, `irc_maxiter` and `irc_files`;
   - a third state fed (with `_feed`) only the lines of `_IRC_OUT` up to and including the backward row `    0     -100.000797 …` has `len(irc_rows["backward"]) == 1` and `_in_block()` true.
8. Run the first three verification commands, and the prototype check on the two real IRCs: `cd /home/madsr2d2/masterThesis && .venv/bin/python -c "from pathlib import Path; from orcamon.core.parser import JobState, read_appended; import sys
for d in sys.argv[1:]:
    s = JobState(path=Path(d), stem='job'); read_appended(s, Path(d) / 'job.out')
    print(d.split('/')[-1], {k: len(v) for k, v in s.irc_rows.items()}, s.irc_maxiter, s.irc_rows['backward'][-1].de_kcal)" computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_b973c-xtb/irc_bridge computational/C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_r2scan3c-xtb/irc`.

**Acceptance:**
- The real-tree command prints `irc_bridge {'forward': 50, 'backward': 50} ['forward', 'backward'] -61.448643` and `irc {'forward': 20, 'backward': 20} ['forward', 'backward'] -17.416051`.
- `test_irc_rows_are_parsed` passes its 7 checks; `validate.MARKER_CASES all read` still passes.
- `test_monitor.py`: 168 `pass`, `0 failure(s)`. `test_orcamon.py`: 178 `pass`. `data/test_curve_metrics.py`: `0 failure(s)`.

### PA2 — Build the IRC path from the rows and the trajectories

**Why:** the rows give each point's energy; the trajectory files give its geometry. Both readers — the CLI and the TUI — need the same ordered list of points, built once and cached, so they cannot disagree.
**Depends on:** PA1.
**Files:** `tools/orcamon/src/orcamon/core/geometry.py`, `tools/orcamon/src/orcamon/core/units.py`, `tools/orcamon/src/orcamon/core/paths.py` (new), `tools/orcamon/src/orcamon/core/job.py`, `tools/orcamon/src/orcamon/core/report.py`, `test_monitor.py`.

**Do:**
1. `core/units.py`: add `KCAL_TO_KJ = 4.184` with a comment (the thermochemical calorie; ORCA prints IRC dE in kcal/mol).
2. `core/geometry.py`: add `read_xyz_frames(path: Path) -> list[tuple[str, list]]` after `read_xyz`: every COMPLETE frame as `(comment, atoms)`, atoms as `read_xyz` makes them. Read the whole text (`errors="replace"`), walk it frame by frame (skip blank lines between frames; a first line that is not an integer ends the walk), and stop at a frame whose atom lines run past the end or fail to parse — a file ORCA is rewriting ends in a partial frame, which is dropped. A missing file gives `[]`.
3. Create `core/paths.py` (standard library only; module docstring: what a reaction path is here and that D2–D6 hold) with:
   ```python
   @dataclass(eq=False)
   class PathPoint:
       """One point of an IRC or NEB path. Stands in for a GeometryPoint in the
       panes and in `geom`, with no cycle: `label` is what it is called."""
       kind: str                 # "irc" | "neb"
       index: int                # IRC: signed (backward < 0, TS 0, forward > 0); NEB: the image
       label: str
       energy: float | None
       de_kj_mol: float          # IRC: from the TS; NEB: from image 0
       atoms: list
       monitors: dict = field(default_factory=dict)
       energy_label: str | None = None
       source: str = ""          # the file the atoms came from, "" when none
       scan_step: None = None
       cycle: int = 0
       key: None = None

   @dataclass
   class PathView:
       kind: str
       points: list              # PathPoint, in path order
       focus: int                # index into `points` the chart rings while following (D6)
   ```
   a private frame cache `_frames(path) -> list` (keyed on `(str(path), st_mtime_ns, st_size)`, at most 8 entries, oldest dropped first; a missing file gives `[]`), and `reaction_path(job) -> PathView | None`, cached on the job (step 5). For an IRC (`state.irc_rows` non-empty): build the points per D3–D5 — backward rows in REVERSE order, then the TS, then forward rows; `de_kj_mol = row.de_kcal * KCAL_TO_KJ`; `monitors = dict(zip(state.irc_monitors, row.monitors))`; `energy_label = state.final_energy_label` for every point; the TS energy is the first `state.points` entry with an energy (else `None`) and its dE `0.0`. Atoms per D5: the full file when `len(frames) == len(backward) + 1 + len(forward)`, frame `j` for point `j` of the list (`source` = that file's name); else the per-direction files and, for the TS, `read_xyz(job.state.path / job.input.coords_file)` when `job.input` names one and its atom count equals that of the first frame found in either direction file (`source` = that file's name); a missing frame gives `atoms=[]`. Focus per D6 (`state.normal_completion or state.crashed_marker` = ended). Return `None` when the job has no IRC rows (NEB is PB1's).
4. In `core/paths.py` also add `path_summary(state) -> tuple[str, str] | None` — from the parsed state alone, never a file (`ls` calls it for every job). For an IRC: per direction, in the order `backward`, `forward`, skipping a direction with no rows, the text `f"{direction} {n} point{'s' if n != 1 else ''}" + (" (MaxIter)" if direction in state.irc_maxiter else "") + f", end {last.de_kcal * KCAL_TO_KJ:+.1f} kJ/mol"`; joined with `" · "`, prefixed `"IRC        "` and suffixed `" (from the TS)"`. Progress: `f"irc B{len(backward)}/F{len(forward)}"`. Return `None` with no IRC rows.
5. `core/job.py`: add `self.path_view: tuple | None = None` in `Job.__init__` (comment: `paths.reaction_path`'s cache, `(signature, PathView)`). `reaction_path` returns the cached view when the signature `(len(forward), len(backward), state.normal_completion, state.crashed_marker, stamp(full), stamp(forward file), stamp(backward file), stamp(coords file))` is unchanged, where `stamp(p)` is `(st_mtime_ns, st_size)` or `None`; the SAME `PathView` object (and the same `PathPoint` objects) are returned until it changes.
6. `core/report.py`: in `geometry_shown`, right after the `if job is None` return, add `if isinstance(point, PathPoint): return (point.atoms, point) if point.atoms else ([], None)`; in `describe_point`, before the `FileGeometry` test, `if isinstance(point, PathPoint): return point.label`.
7. Add `_irc_tree` (§ Facts) to `test_monitor.py` and `test_the_irc_path_is_built()`, registered right after `test_irc_rows_are_parsed`. Build the tree in a temp dir `d`, make `job = Job(d, "job", d.parent, label="irc")`, `job.refresh(Liveness(False, "process"))`, `view = reaction_path(job)`. Exactly 11 `check` calls, against § Facts "Hand-derived values for the synthetic IRC":
   - `view.kind == "irc"` and the indices;
   - the labels;
   - the dE list (`round(p.de_kj_mol, 3)`);
   - the O z per point (`p.atoms[0][3]`), with the full trajectory present;
   - `view.focus == 2 and view.points[2].energy == -100.0`;
   - `view.points[-1].monitors == {"B(O 0,H 1)": 0.99, "B(O 0,H 2)": 0.94}`;
   - `reaction_path(job) is view` (nothing changed);
   - after deleting `job_IRC_Full_trj.xyz`, a fresh `reaction_path(job)` is a different object whose O z per point is still `[0, 1, 2, 3, 4, 5]`;
   - `read_xyz_frames` on a file holding 2 complete 3-atom frames and then `"3\ncomment\nO 0 0\n"` returns 2 frames;
   - `geometry_shown(job, p) == (p.atoms, p)` and `describe_point(p) == "IRC TS"` for the TS point, and `geometry_shown(job, PathPoint("irc", 9, "IRC forward 8", None, 0.0, [])) == ([], None)`;
   - a second tree in another temp dir whose `job.out` is `_IRC_OUT` cut after the forward block's last row (no backward block, no termination line) gives `focus == 3` (the newest forward point, of 4 points).
8. Run the first three verification commands.

**Acceptance:**
- `test_the_irc_path_is_built` passes its 11 checks.
- `test_monitor.py`: 179 `pass`, `0 failure(s)`. `test_orcamon.py`: 178 `pass`. `data/test_curve_metrics.py`: `0 failure(s)`.
- `core/paths.py` imports nothing outside the standard library and `orcamon.core` (`test_orcamon.test_the_core_needs_only_the_standard_library` still passes).

### PA3 — Say where the IRC is, and flag a direction that hit MaxIter

**Why:** `show irc_bridge` says `finished` with one energy and no flag, though both directions stopped at MaxIter 50 without reaching a minimum — exactly what the job's own input comment says went wrong with the previous IRC.
**Depends on:** PA2.
**Files:** `tools/orcamon/src/orcamon/core/status.py`, `tools/orcamon/src/orcamon/core/report.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/README.md`, `tools/orcamon/src/orcamon/SKILL.md`, `.claude/skills/orcamon/SKILL.md` (regenerated), `test_monitor.py`, `test_orcamon.py`.

**Do:**
1. `core/status.py`: insert `"irc_not_converged"` into `FLAG_CODES` right after `"scan_incomplete"`; in `attention()`, right after the `scan_incomplete` flag, add: when `state.irc_maxiter` is non-empty, `Flag("irc_not_converged", f"IRC {' and '.join(state.irc_maxiter)} hit MaxIter{maxiter} before a minimum")`.
2. `core/report.py`: add private fields `_path_summary: str | None = field(default=None, repr=False)` and `_path_progress: str | None = field(default=None, repr=False)` after `_converged_reason`; in `build_report`, `summary = path_summary(state)` and set both from it (or `None`); in `report_lines`, right before the `"energy"` line, add `("path", [(report._path_summary, None)])` when it is not None.
3. `cli/commands.py::_progress`: return `r._path_progress` first when it is not None.
4. `tools/orcamon/README.md`: add `irc_not_converged` after `scan_incomplete` in the flag-code list, and one sentence after it: "An IRC direction that stopped at MaxIter before reaching a minimum raises `irc_not_converged`." `tools/orcamon/src/orcamon/SKILL.md`: add the attention-table row `` | `irc_not_converged` | an IRC direction hit MaxIter before reaching a minimum | `orcamon energies JOB` | `` after the `scan_incomplete` row. Then run `cd /home/madsr2d2/masterThesis && .venv/bin/orcamon skill install --project /home/madsr2d2/masterThesis` (add `--force` only if it refuses because the installed copy is stale).
5. `test_monitor.py`: add `test_the_irc_is_summarised_and_flagged()`, registered right after `test_the_irc_path_is_built`. Exactly 5 `check` calls, on a state fed `_IRC_OUT`:
   - `path_summary(state) == ("IRC        backward 2 points, end -6.3 kJ/mol · forward 3 points (MaxIter), end -12.6 kJ/mol (from the TS)", "irc B2/F3")`;
   - `render_plain(build_report(_report_job(state, "! B97-3c IRC\n%irc MaxIter 3 end\n* xyzfile 0 1 start.xyz\n", Liveness(False, "process"))))` contains that summary line;
   - that report's `attention` contains `{"code": "irc_not_converged", "message": "IRC forward hit MaxIter (3) before a minimum"}`;
   - a state fed `_IRC_OUT` with the MaxIter banner's three lines removed raises no `irc_not_converged`;
   - `FLAG_CODES.index("irc_not_converged") == FLAG_CODES.index("scan_incomplete") + 1`.
6. `test_orcamon.py`: add `_IRC_OUT` and `_irc_tree` (§ Facts) and `test_ls_says_where_an_irc_is()`, registered right after `test_ls_prints_labels_an_agent_can_pass_back`. In a temp root with `_irc_tree(root / "irc")`, `_orcamon(["--root", str(root), "--liveness", "mtime", "--no-cache", "ls"])`; exactly 1 `check`: the `irc` row contains `irc B2/F3` and `irc_not_converged`.
7. Run the first three verification commands, and `cd /home/madsr2d2/masterThesis && XDG_CACHE_HOME=/tmp/orcamon_paths/xdg .venv/bin/orcamon --root computational show irc_bridge`.

**Acceptance:**
- The real `show irc_bridge` contains the line `IRC        backward 50 points (MaxIter), end -257.1 kJ/mol · forward 50 points (MaxIter), end -5.3 kJ/mol (from the TS)` and the flag `! IRC forward and backward hit MaxIter (50) before a minimum`.
- `XDG_CACHE_HOME=/tmp/orcamon_paths/xdg .venv/bin/orcamon --root computational ls --max-lines 200 | grep 'irc_bridge'` shows `irc B50/F50` in the PROGRESS column.
- `test_monitor.py`: 184 `pass`; `test_orcamon.py`: 179 `pass` (including `test_this_repository_has_the_current_skill` and `test_report_keys_are_stable`); both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

### PA4 — `energies`, `geom --point` and `snapshot --point` on a path

**Why:** an agent asked "what does the IRC connect" needs the path's energies and the end points' geometries; today `energies` prints one row and `geom` cannot name a path point.
**Depends on:** PA3.
**Files:** `tools/orcamon/src/orcamon/cli/__init__.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/README.md`, `tools/orcamon/src/orcamon/SKILL.md`, `.claude/skills/orcamon/SKILL.md` (regenerated), `test_orcamon.py`.

**Do:**
1. `cli/__init__.py::build_parser`: add to `geom` and to `snapshot` the option `--point N` (`type=int`, `metavar="N"`, help `"a point on an IRC or NEB path: the signed IRC point (0 is the TS, negative is backward) or the NEB image"`).
2. `cli/commands.py`: add `_path_point(job, n)` that calls `paths.reaction_path(job)` and raises `UsageError("--point needs an IRC or NEB job; this job has no reaction path")` when it is None, `UsageError(f"point {n} is not on the path (points {first.index} to {last.index})")` when no point has `index == n`, and `UsageError(f"{p.label}: geometry not written yet")` when that point has no atoms; else returns the point.
3. `cmd_geom`: when `args.point is not None`, raise `UsageError("--point cannot be combined with --cycle or --step")` if either is set, else `point = _path_point(job, args.point)` and skip `_find_point`; the rest of the command (region, energy, XYZ comment from `describe_point`, JSON with `getattr(point, "source", "output")`) is unchanged. `cmd_snapshot`: when `args.point is not None`, raise `UsageError("--point cannot be combined with --mode")` if `args.mode` is set, else draw `_path_point(job, args.point)` instead of the latest geometry.
4. `cmd_energies`: first, `view = paths.reaction_path(job)`; when it is not None, `rows = view.points[-args.last:] if args.last else view.points`. JSON: `_emit_json({"job": job.label, "mode": view.kind, "truncated": False, "rows": [{"point": p.index, "label": p.label, "energy_eh": p.energy, "dE_kj_mol": p.de_kj_mol, "monitors": p.monitors} for p in rows]})`. Text: a header `f"{'point':>6}  {'label':<20}  {'energy/Eh':>18}  {'dE (kJ/mol)':>12}"` plus `f"  {name:>16}"` per monitor name (`state.irc_monitors` for an IRC, none for an NEB); one line per point `f"{p.index:>6}  {p.label:<20}  {energy}  {p.de_kj_mol:>+12.2f}"` with `energy = f"{p.energy:>18.9f}"` or `f"{'-':>18}"`, plus `f"  {v:>16.2f}"` (or `f"  {'-':>16}"`) per monitor; then `dE from the TS` (IRC) or `dE from image 0` (NEB); then the existing `energies are the {label} total` line when the job has a label. Return through `_job_exit(job)` as the other branches do.
5. `tools/orcamon/src/orcamon/SKILL.md`, in "Looking closer", one prose paragraph (no fenced block): "On an IRC or NEB job, `orcamon energies JOB` is the path, one row per point: the dE is from the TS (IRC) or from image 0 (NEB), and an IRC row carries its monitored internals. `orcamon geom JOB --point N` gives one point's geometry: N is the signed IRC point (0 is the TS, negative is backward) or the NEB image. An NEB path is the CURRENT one: ORCA rewrites it every iteration." In `tools/orcamon/README.md`, replace the pasted `orcamon geom` and `orcamon snapshot` `--help` blocks with the new output of `.venv/bin/orcamon geom --help` and `.venv/bin/orcamon snapshot --help`, and add a short "Reaction paths" section after "Saving a picture" saying the same as the skill paragraph plus that the TUI chart scrubs the path. Then rerun `orcamon skill install --project /home/madsr2d2/masterThesis`.
6. `test_orcamon.py`: add `test_a_path_is_listed_and_drawn()`, registered right after `test_snapshot_writes_a_png`. Temp root with `_irc_tree(root / "irc")` and a plain job `root / "plain" / "job.inp"` = `"! B97-3c SP\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n"`; every call passes `["--root", str(root), "--liveness", "mtime", "--no-cache", …]`. Exactly 8 `check` calls:
   - `energies irc --json`: `mode == "irc"`, 6 rows, `rows[0] == {"point": -2, "label": "IRC backward 1", "energy_eh": -100.00239, "dE_kj_mol": <-6.276 after round 3>, "monitors": {"B(O 0,H 1)": 0.94, "B(O 0,H 2)": 0.98}}` (compare `dE_kj_mol` rounded);
   - `energies irc` (text) contains `IRC TS` and `dE from the TS`;
   - `geom irc --point 1`: exit 0, line 2 contains `IRC forward 0`, and line 3's z field is `3.00000000`;
   - `geom irc --point 9`: exit 2, stderr contains `not on the path (points -2 to 3)`;
   - `geom irc --point 0 --cycle 1`: exit 2;
   - `geom plain --point 0`: exit 2, stderr contains `no reaction path`;
   - `snapshot irc --point -2 -o <tmp>/p.png`: exit 0, the file exists, stdout contains `IRC backward 1`;
   - `snapshot irc --point 0 --mode 6 -o <tmp>/q.png`: exit 2.
7. Run the first three verification commands, and `cd /home/madsr2d2/masterThesis && XDG_CACHE_HOME=/tmp/orcamon_paths/xdg .venv/bin/orcamon --root computational energies irc_bridge --json | .venv/bin/python -c "import json,sys; d=json.load(sys.stdin); r=d['rows']; print(d['mode'], len(r), r[0]['point'], r[0]['label'], round(r[0]['dE_kj_mol'], 2), r[50]['label'], r[50]['energy_eh'], r[-1]['point'], round(r[-1]['dE_kj_mol'], 2))"` and `XDG_CACHE_HOME=/tmp/orcamon_paths/xdg .venv/bin/orcamon --root computational geom irc_bridge --point -50 | head -1`.

**Acceptance:**
- The real `energies` command prints `irc 101 -50 IRC backward 49 -257.1 IRC TS -1011.288471006717 50 -5.27`; the `geom --point -50` prints `137`.
- `test_a_path_is_listed_and_drawn` passes its 8 checks.
- `test_orcamon.py`: 187 `pass`; `test_monitor.py`: 184 `pass`; both `0 failure(s)` (the skill tests included). `data/test_curve_metrics.py`: `0 failure(s)`.

### PA5 — Scrub the path in the TUI

**Why:** the person following the agent should see the reaction: the chart's dots should be the path's points, and left/right should walk the geometry pane along it. Today an IRC's chart shows SCF iterations and the pane shows the 23-atom TS block.
**Depends on:** PA4.
**Files:** `tools/orcamon/src/orcamon/tui/app.py`, `test_orcamon.py`.

**Do:**
1. `ConvergencePlot.__init__`: add `self._mode: str | None = None` and `self._view = None`.
2. Change `_series(self, state)` to `_series(self, job)`: first `view = paths.reaction_path(job)`, store it in `self._view`, and when it is not None return `(view.kind, view.points)`; otherwise run the existing body on `job.state`. Update the one caller in `show_job`, and store the mode: `self._mode = mode`.
3. `_signature_of(job)`: append `id(paths.reaction_path(job))` to the tuple.
4. `show_job`: in the non-`scf` branch, for mode `irc` or `neb` use `xs = [p.index for p in points]`, `ys = [p.de_kj_mol for p in points]`; title `f"IRC: {len(points)} points, dE from the TS"` and xlabel `"point (backward < 0, TS 0, forward > 0)"` for `irc`; title `f"NEB: {len(points)} images, dE from image 0"` and xlabel `"image"` for `neb`; the existing ylabel, scatter and selection ring follow unchanged. The other modes keep their current `ref`/`ys`/titles.
5. Selection in path modes (`self._mode in ("irc", "neb")`):
   - `_place_selection`: when not following and `self._selected` is set, re-find it by `p.index == self._selected.index` (the points are rebuilt when the path grows, so identity cannot be used); otherwise `self.selected_index = self._view.focus`. Set `self._selected` to the chosen point.
   - `selected_point()`: return `self._points[self.selected_index]` (or None when there is no index) — also while following, so the pane draws the focus point and not the job's printed geometry.
   - `_select(index)`: set `self._following_latest = False` (a person who scrubbed stays where they scrubbed, even at the last point).
   - `action_select_latest` (`end`): set `self._following_latest = True`, then `self.app.update_detail()` (the focus is re-placed by `_place_selection`).
   The non-path modes keep their current behaviour exactly.
6. `test_orcamon.py`: add `test_the_tui_scrubs_a_path()`, registered right after `test_the_tui_runs_headless`. Temp root holding only `_irc_tree(root / "irc")`; `MonitorApp(root, liveness="mtime", graphics="text", notify_mode="off")`, run headless at `(180, 50)` as `test_the_tui_runs_headless` does, wait until the job is parsed and the chart has points. Let `title = lambda: app.query_one("#geometry").border_title`. Exactly 5 `check` calls:
   - `chart._mode == "irc" and len(chart._points) == 6`;
   - following, the title ends with `· IRC TS`;
   - after `chart.focus()` and `right`: the title ends with `· IRC forward 0`;
   - after `left`, `left`: the title ends with `· IRC backward 0`;
   - after `end`: the title ends with `· IRC TS`.
7. Run the first three verification commands.

**If unsure:** if the geometry title carries more after the point label than ` · <label>` (another suffix added since `ea580b4`), check `· IRC TS" in title` instead of `endswith`, and report it.

**Acceptance:**
- `test_the_tui_scrubs_a_path` passes its 5 checks; `test_the_tui_runs_headless` still passes unchanged (a non-path job's chart, ring and `end` behave as before).
- `test_orcamon.py`: 192 `pass`; `test_monitor.py`: 184 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

# Phase 2 — NEB

### PB1 — Parse the NEB rows and build its path

**Why:** an NEB prints one row per iteration (the highest or climbing image, its barrier and forces) and rewrites its current path, with every image's energy, to the MEP trajectory it names; orcamon reads neither (`show neb_anion` shows no progress at all).
**Depends on:** PA5.
**Files:** `tools/orcamon/src/orcamon/core/parser.py`, `tools/orcamon/src/orcamon/core/paths.py`, `tools/orcamon/src/orcamon/validate.py`, `test_monitor.py`.

**Do:**
1. `core/parser.py`: add `_NEB_HEADER_RE`, `_NEB_ROW_RE`, `_NEB_TRAJECTORY_RE` (§ Facts) after the IRC regexes; append `r"|E\(0\)|Current trajectory will be written"` to `_RARE_MARKERS_RE`. Add before `JobState`:
   ```python
   @dataclass
   class NebRow:
       """One NEB iteration: the highest (HEI) or climbing (CI) image, its energy
       above image 0 in Eh, and the perpendicular forces."""
       iteration: int
       phase: str          # "HEI" | "CI"
       image: int
       barrier_eh: float
       max_fp: float
       rms_fp: float
   ```
   and `JobState` fields `neb_rows: list = field(default_factory=list)`, `neb_file: str | None = None`, private `_in_neb_rows: bool = False`, `_neb_phase: str | None = None`, `_neb_seen_row: bool = False`.
2. In `_feed_marker`: on `_NEB_HEADER_RE.match(line.strip())` set `_neb_phase` to group 1, `_in_neb_rows = True`, `_neb_seen_row = False`; on `_NEB_TRAJECTORY_RE` set `neb_file`. Add `_feed_neb_row(line)`, called from `_process` after `_feed_marker` when `_in_neb_rows`: a `_NEB_ROW_RE` match appends `NebRow(int(g2), self._neb_phase, int(g3), float(g4), float(g5), float(g6))` and sets `_neb_seen_row`; else nothing before the first row, and the block closes after it. Add `_in_neb_rows` to `_in_block()`.
3. `core/paths.py`: extend `reaction_path(job)` for an NEB (`state.neb_file` set and the MEP file has frames): image `i` is frame `i`, its energy the float after `E` in the comment (`re.search(r"\bE\s+(-?\d+\.\d+)", comment)`); a frame with no energy ends the path there; frame 0 without one means `None` is returned. `de_kj_mol = (E_i - E_0) * EH_TO_KJ_PER_MOL`; label `f"NEB image {i}"` plus `f" ({last.phase})"` on the latest row's image; `source` = the MEP file name; `energy_label = state.final_energy_label`; focus = the latest row's image when it is within the path, else 0. Signature for the cache: `(len(state.neb_rows), stamp(MEP file))`. `path_summary`: for an NEB with rows, `f"NEB        iteration {r.iteration} · {'climbing' if r.phase == 'CI' else 'highest'} image {r.image} · E({r.phase})-E(0) {r.barrier_eh * EH_TO_KJ_PER_MOL:+.1f} kJ/mol"` and progress `f"neb {r.iteration} {r.phase}{r.image}"`, with `r` the latest row. IRC handling is unchanged.
4. `validate.py`: append `("Optim.  Iteration  CI   E(CI)-E(0)   max(|Fp|)   RMS(Fp)    dS", lambda s: s._in_neb_rows and s._neb_phase == "CI")` and `("Current trajectory will be written to    ....  job_MEP_trj.xyz", lambda s: s.neb_file == "job_MEP_trj.xyz")` to `MARKER_CASES`.
5. `test_monitor.py`: add `_NEB_OUT` and `_neb_tree` (§ Facts) and `test_the_neb_is_parsed_and_built()`, registered right after `test_the_irc_is_summarised_and_flagged`. Exactly 9 `check` calls, against § Facts "Hand-derived values for the synthetic NEB":
   - the state fed `_NEB_OUT` has 3 rows and the last equals `NebRow(2, "CI", 2, 0.008, 0.015, 0.003)`;
   - phases `["HEI", "HEI", "CI"]`;
   - `neb_file == "job_MEP_trj.xyz"`;
   - on the tree's job (`Job(d, "job", d.parent, label="neb")`, `refresh(Liveness(False, "process"))`), `reaction_path` has kind `neb`, images `[0, 1, 2, 3]` and the labels;
   - the dE list, `round(…, 4)`;
   - `focus == 2`;
   - `path_summary(state) == ("NEB        iteration 2 · climbing image 2 · E(CI)-E(0) +21.0 kJ/mol", "neb 2 CI2")`;
   - a MEP file whose 4th frame is cut after its first atom line gives 3 images;
   - with the MEP file deleted, `reaction_path` is None and `path_summary` is unchanged.
6. Run the first three verification commands.

**Acceptance:**
- `test_the_neb_is_parsed_and_built` passes its 9 checks; `validate.MARKER_CASES all read` still passes; `test_the_irc_path_is_built` is unchanged and passes.
- `test_monitor.py`: 193 `pass`; `test_orcamon.py`: 192 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

### PB2 — The NEB through the CLI and the TUI, then every gate

**Why:** PA3–PA5 wired the CLI and the chart to `reaction_path` and `path_summary`, which now answer for an NEB too; this task proves it end to end on the synthetic NEB and the real one, and runs every gate.
**Depends on:** PB1.
**Files:** `test_orcamon.py`; any of `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/src/orcamon/tui/app.py` only if a check below fails because PA3–PA5 assumed an IRC (say so in the report).

**Do:**
1. `test_orcamon.py`: add `_NEB_OUT` and `_neb_tree` (§ Facts) and `test_an_neb_is_listed_and_drawn()`, registered right after `test_a_path_is_listed_and_drawn`. Temp root holding `_neb_tree(root / "neb")`; every CLI call passes `["--root", str(root), "--liveness", "mtime", "--no-cache", …]`. Exactly 5 `check` calls:
   - `energies neb --json`: `mode == "neb"`, 4 rows, `[round(r["dE_kj_mol"], 4) for r in rows] == [0.0, 13.1275, 21.004, 5.251]` and `rows[2]["label"] == "NEB image 2 (CI)"`;
   - `geom neb --point 2`: exit 0 and line 3's z field is `2.00000000`;
   - `show neb` contains `NEB        iteration 2 · climbing image 2 · E(CI)-E(0) +21.0 kJ/mol`;
   - `ls`: the `neb` row contains `neb 2 CI2`;
   - headless TUI on that root (as `test_the_tui_scrubs_a_path`): the chart's `_mode == "neb"` and the geometry title ends with `· NEB image 2 (CI)`.
2. Run all four verification commands, and on the real tree (`cd /home/madsr2d2/masterThesis`, `XDG_CACHE_HOME=/tmp/orcamon_paths/xdg`):
   - `.venv/bin/orcamon --root computational energies neb_anion --json | .venv/bin/python -c "import json,sys; d=json.load(sys.stdin); r=d['rows']; print(d['mode'], len(r), r[4]['label'], round(r[4]['dE_kj_mol'], 2), [round(x['dE_kj_mol'], 2) for x in r])"`;
   - `.venv/bin/orcamon --root computational show neb_anion`;
   - `.venv/bin/orcamon --root computational ls --max-lines 200 | grep -E 'neb_anion|irc'`;
   - `.venv/bin/orcamon --root computational show geometry_r2scan3c-xtb/irc`.

**Acceptance:**
- The NEB `energies` prints `neb 10 NEB image 4 (CI) 170.59 [0.0, 13.1, 110.88, 151.54, 170.59, 159.1, 171.67, 173.02, 148.25, 125.51]`.
- `show neb_anion` contains `NEB        iteration 14 · climbing image 4 · E(CI)-E(0) +170.6 kJ/mol`; its `ls` row shows `neb 14 CI4`; the two IRC rows show `irc B50/F50` and `irc B20/F20`.
- `show geometry_r2scan3c-xtb/irc` contains `IRC        backward 20 points (MaxIter), end -72.9 kJ/mol · forward 20 points (MaxIter), end +5.9 kJ/mol (from the TS)` and `! IRC forward and backward hit MaxIter (20) before a minimum`.
- `test_an_neb_is_listed_and_drawn` passes its 5 checks.
- `test_orcamon.py`: 197 `pass`; `test_monitor.py`: 193 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`. `run_gates.py`: `38 gates`, `0 failed`.
- `git diff ea580b4 -- '*test*'` across the whole plan shows only added tests, fixtures, `__main__` lines and the six `MARKER_CASES` rows; no existing assertion or expected value changed.

---

## Out of scope

1. Keeping real normal modes as well as imaginary ones (the next plan, with picking a mode in the TUI).
2. Any marker without a real sample (D8): IRC convergence in a direction, NEB convergence, NEB-TS's closing TS optimization and its `*_NEB-TS_converged.xyz`, ZOOM-NEB specifics. A finished NEB-TS will show its NEB path; its TS stage is read as today.
3. The `IRC PATH SUMMARY` block, `*_MEP_ALL_trj.xyz` (the path history), `job.interp`, `job.NEB.log`, and the `*.QMRegion*` trajectory variants.
4. Path data in `show --json`/`ls --json` (would bump `REPORT_SCHEMA`; decide with the user later).
5. Fixing `n_atoms` for multilayer jobs whose `.out` prints only the QM block (`(21 of 23 atoms)` on `irc_bridge`).
6. Animating along the path, playing it as a movie, or overlaying two points.
7. GOAT, MD and scan-of-scan trajectories.
8. Any change to how a non-path job is read, listed, charted or drawn.

## Progress template

```markdown
# orcamon — reaction paths: IRC and NEB — Progress

Plan: `PLAN_ORCAMON_PATHS.md`. Branch `master` from `master` (`<sha of the plan commit>`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| PA1 Parse the IRC rows | TODO | | | |
| PA2 Build the IRC path | TODO | | | |
| PA3 Summary and irc_not_converged | TODO | | | |
| PA4 energies, geom/snapshot --point | TODO | | | |
| PA5 Scrub the path in the TUI | TODO | | | |
| PB1 Parse the NEB and build its path | TODO | | | |
| PB2 NEB end to end, every gate | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 161 | 178 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

## Backlog
```
