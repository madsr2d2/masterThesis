# orcamon — the agent changes a geometry, checks an input, and reads what a run did to the bonds (AG1–AG6)

**Written:** 2026-10-02. **Verified against:** `master` at `a825532` (the Steer-fixes run, `PLAN_ORCAMON_STEER_FIXES.md`, is committed and complete). Every rule below was prototyped on scratch copies against this commit and run over all 80 jobs under `computational/`; the hand-derived values in § Facts were measured on those prototypes.
**Branch:** `master` (the agents commit straight to `master`).
**Progress file:** `PLAN_ORCAMON_ACT_GEOMETRY_PROGRESS.md` (template at the end).

The person and the agent now share a view (`orcamon view`, `orcamon focus`). What the agent still cannot do safely is the thing the person asks for: "rotate this bond 55 degrees", "stretch this distance to 2.1", "is this input right before it runs", "did that run do something to the bonds it should not have". This plan adds three commands. `orcamon edit` rotates a bond, sets a dihedral or sets a distance, moving one side of the molecule rigidly, and writes a NEW XYZ file, never touching an existing one; it refuses an edit that is impossible (atoms not bonded, a bond in a ring) or that creates a new bond, and always says how many atoms moved, how many of them are QM and environment, and the closest new contact. `orcamon check` finds what would make an input fail or mislead before it is run: no geometry, an electron count the multiplicity cannot have, an atom number out of range anywhere in the input, more processes or memory than the machine has, an output that would be overwritten, and a QM atom bonded to an environment hydrogen. `orcamon connectivity` lists the bonds a run formed and broke, from the first geometry ORCA printed to the one asked for, or from an IRC's TS (an NEB's image 0) to a point on the path. **The goal a deviation must still meet:** an agent can turn the person's words into a checked geometry file and a checked input, and can say what a run did to the bonds, with every number reproducible from the job's files; nothing here starts, edits or stops an ORCA job, writes into a job's directory unless the agent names that path, or overwrites any file; and every existing command, key, output and JSON document behaves as before.

## Execution

The opencode orchestrator runs this plan. These are this project's rules; the loop itself belongs to the agents and is not repeated here.

### Verification commands

| What | Working directory | Command | Expected now | Time |
|---|---|---|---|---|
| renderer and parser tests | `/home/madsr2d2/masterThesis` | `.venv/bin/python test_monitor.py` | 246 `pass` lines, last line `0 failure(s)` | ~3 s |
| CLI, TUI and skill tests | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_act/xdg XDG_STATE_HOME=/tmp/orcamon_act/state .venv/bin/python test_orcamon.py` | 263 `pass` lines, last line `0 failure(s)` | ~30 s |
| duplicate-definition guard | `/home/madsr2d2/masterThesis` | `.venv/bin/python data/test_curve_metrics.py` | last line `0 failure(s)` | ~45 s |
| every gate | `/home/madsr2d2/masterThesis` | `HERDR_ENV=0 XDG_CACHE_HOME=/tmp/orcamon_act/xdg XDG_STATE_HOME=/tmp/orcamon_act/state .venv/bin/python run_gates.py` | last block says `38 gates in … on 8 jobs, 0 failed` | 500–750 s; give the bash call a timeout of at least 900000 ms |

Count `pass` lines with `.venv/bin/python test_monitor.py | grep -c '^  pass'` (the same for `test_orcamon.py`, with its environment variables). Each task states the counts it must end at.

Every task runs the first three commands. Only AG6 runs `run_gates.py`: every task touches only `tools/orcamon/`, `test_monitor.py`, `test_orcamon.py` and `.claude/skills/orcamon/SKILL.md`, which no other gate reads except the duplicate guard.

**Do not add a new root `test_*.py` file.** `test_root_documents.py` pins the gate count to `CLAUDE.md`, which is the user's; a new gate file would change the count and break that gate. New tests go in `test_monitor.py` or `test_orcamon.py`.

`test_monitor.py` and `test_orcamon.py` run their tests from an explicit list in `if __name__ == "__main__":` at the bottom of each file. A new test function that is not added to that list never runs.

### Stop conditions

Stop and report NEEDS_USER (implementer) or BLOCK (reviewer) when:
- a hand-derived expected value in this plan does not hold after the change is implemented as written (report the measured value; do not change the expected value or the design to make it pass);
- any existing test fails and the fix would mean changing that test's assertion, expected value or tolerance;
- a real-tree value differs from this plan's (report it, and say whether that job's files changed since `a825532`);
- in AG6, a gate that does not involve orcamon fails (report it; do not fix it);
- the working tree has changes not made by this plan's tasks.

A different name or path for the same thing is not a stop: adapt, and report it as a deviation.

### Never

- Never edit `CLAUDE.md`, the root `README.md` or `.claude/skills/run-orca/SKILL.md`; they are the user's. (`tools/orcamon/README.md` and `tools/orcamon/src/orcamon/SKILL.md` are orcamon's; `.claude/skills/orcamon/SKILL.md` changes ONLY through `orcamon skill install --project .`.)
- Never `git add -A`, `git add .` or `git commit -a`; stage the files the task lists, by path.
- Never push.
- Never add a new root `test_*.py` gate file.
- Never define the same top-level name in both `test_monitor.py` and `test_orcamon.py`, private or not: the duplicate guard reads private names too. A fixture both files need lives in `test_monitor.py` and `test_orcamon.py` imports it.
- Never let a test or a manual command read or write the user's own state or cache: every manual `orcamon` command sets `XDG_CACHE_HOME=/tmp/orcamon_act/xdg` and `XDG_STATE_HOME=/tmp/orcamon_act/state`. Anything that starts the TUI sets `HERDR_ENV=0`.
- Never write into `computational/`: every real-tree `orcamon edit` in this plan writes its `-o` file under `/tmp/orcamon_act/`. Never start ORCA, and never kill or signal a process.
- Never change an existing test's assertion, expected value or tolerance.
- Never change `geometry.bonds()` or `BOND_CUTOFF`: they are the renderer's, and the new bond rule (D3) is a separate function.
- Never add a dependency. `core` and `cli` import only the standard library at module level (`test_the_core_needs_only_the_standard_library` imports every `orcamon.core` and `orcamon.cli` module with the TUI's dependencies refused).
- Never put thesis paths, job names or chemistry names into orcamon's package code or its shipped skill.
- Never change the key set of an existing `--json` document, and never change `REPORT_SCHEMA` or `steer.STEER_SCHEMA`. `edit --json`, `check --json` and `connectivity --json` are new documents.
- Never add a fenced code block containing an `orcamon …` line to `tools/orcamon/src/orcamon/SKILL.md`: `test_orcamon.test_the_skill_examples_run` executes every such line. Mention commands in prose, in backticks.

## Scope

**In scope:** the edit bond rule and the element table; `core/edit.py` (rotate a bond, set a dihedral, set a distance, bond changes); `orcamon edit`; reading every atom number an input uses; `core/inputcheck.py` and `orcamon check`; `JobState.start_atoms`; `orcamon connectivity`; tests, README and the shipped skill.
**Out of scope:** see § Out of scope.

## Decisions

| # | Decision | By, date |
|---|---|---|
| D1 | An edit is never refused for its size. It is refused (exit 2, nothing written) when it is impossible — the two atoms are not bonded, or the bond is in a ring — or when it creates a new bond between a moved and an unmoved atom, unless `--allow-new-bonds` is given. It always reports the atoms moved (with the QM/environment split on a multilayer job), the side that moved and how to move the other, and the closest new contact. | user, 2026-10-02 |
| D2 | The smaller side moves (fewer atoms; a tie moves J's side). `--move ATOM` names an atom of the edited bond whose side moves instead. Fragments connected to neither side never move. | user, 2026-10-01 |
| D3 | The edit and connectivity bond rule is `edit.bonded_pairs`: two heavy atoms are bonded under 1.7 Å (`BOND_CUTOFF`, the renderer's); an atom and a hydrogen under 1.4 Å; two hydrogens under 1.0 Å. The renderer's `bonds()` (1.7 Å for everything except H–H) is unchanged. Why: under 1.7 Å for hydrogen, hydrogen bonds at 1.6–1.7 Å flicker in and out (`pc_bridge` read three "formed" O–H at 1.67–1.70 Å), and a rotation that merely makes a hydrogen bond would be refused. | planner, 2026-10-02 |
| D4 | `--rotate-bond I J DEG` turns the moving side about the axis from the fixed atom to the moving atom by DEG degrees (right-hand rule), so EVERY dihedral X-I-J-Y across the bond (X on I's side, Y on J's side) grows by DEG whichever side moves. `--set-dihedral A B C D DEG` rotates bond B-C by the difference wrapped into [-180, 180); A must lie on B's side and D on C's side. `--set-distance I J ANGSTROM` moves one side along the line from the fixed atom to the moving one; I and J are either bonded (and not in a ring) or in separate fragments — two atoms of one fragment that are not bonded are refused. | planner, 2026-10-02 |
| D5 | The command is `orcamon edit`, not new flags on `geom`: `geom` prints to stdout and changes nothing, `edit` writes a file. `-o FILE` is required and an existing FILE is refused: `edit` never overwrites. The source geometry is chosen as `geom` chooses it (latest, `--cycle`, `--step`, `--point`). | planner, 2026-10-02 |
| D6 | `check` tests the WHOLE system's electron count against the multiplicity (Charge_Total/Mult_Total on a multilayer job) and never the QM region's: the QM region's count depends on the link atoms ORCA creates, which orcamon cannot reproduce (two real jobs, `pc_relax` and `ts_earlier`, read 3 links by bond distance where ORCA's own log says 2). It exits 1 when it finds an error. | planner, 2026-10-02 |
| D7 | `connectivity` compares against the FIRST geometry the output printed (`JobState.start_atoms`, kept however long the run), never against the input's XYZ file: ORCA overwrites `<basename>.xyz` with the latest geometry, and this tree's inputs read `job.xyz` (`rc`'s `job.xyz` equals its final geometry to 8e-7 Å). On an IRC or NEB with `--point N`, it compares point 0 (the TS, image 0) with point N. | planner, 2026-10-02 |
| D8 | Commits go straight to `master`, one per task. | user, 2026-10-02 |

## Facts (verified 2026-10-02)

**Environment.** `/home/madsr2d2/masterThesis/.venv/bin/python` is Python 3.12.3; orcamon is installed editable and its script is on PATH as `orcamon` (`~/.local/bin/orcamon` links to the venv's). `git status --porcelain` was empty at `a825532`. Baseline: `test_monitor.py` 246 pass, `test_orcamon.py` 263 pass (both re-run for this plan), `run_gates.py` 38 gates, 0 failed. This machine has 32 CPUs (`os.cpu_count()`) and 128543 MB of memory (`os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") // 2**20`). README `--help` blocks are pasted at `COLUMNS=90`. `argparse` with `nargs=3` accepts a negative number as a value (`--rotate-bond 133 134 -60` parses to `['133', '134', '-60']`; probed).

**The code** (at `a825532`). `core/geometry.py`: `BOND_CUTOFF = 1.7`, `bonds(atoms)` (grid-binned, H–H skipped), `measure(atoms, indices) -> (kind, value)`, `_check_in_range`, `parse_atom_list`, `read_xyz(path)`, `MAX_PICKED`. `core/orca_input.py`: `parse_input(text) -> JobInput` with `charge`, `mult`, `multilayer`, `layers` (`{"total": (charge, mult)}`, lowercase names), `nprocs`, `maxcore_mb`, `coords_kind`, `coords_file`, `coords_atoms`; `_COMMENT_RE = re.compile(r"#.*$", re.M)` strips comments; it reads no atom numbers today. `core/parser.py`: `JobState.atoms` and `points` (a deque of `GeometryPoint`, `maxlen=HISTORY_LEN = 500`); the coordinate block closes on a blank line, where (lines ~770–774) `point = self.current_point()` and, `if len(self._pending_atoms) > len(point.atoms)`, `self.atoms = self._pending_atoms` and `point.atoms = self._pending_atoms` — the largest block of a cycle wins, because a multilayer job prints the whole system and then QM-only subsets. The state cache is keyed by the hash of `parser.py`'s source, so a new field invalidates every entry by itself. `cli/commands.py`: `_with_job`, `UsageError` (exit 2, `orcamon: <message>` on stderr), `_emit_json`, `_job_exit(job)`, `describe_point`, `_find_point`, `_path_point(job, n)`, `_requested_point(args, job)` (line 499: `--point`/`--cycle`/`--step`/latest/file geometry, refusing as `geom` does), `_point_doc(point)` (line 527), `cmd_measure` (line 538, the pattern for a command over one geometry). `cli/__init__.py`: `add(name, help_, example, func)`, `_job(p)`, `_output(p, lines=…)`, `EXIT_CODES` (`1`: `the command worked and the job (for ls: any listed job) is failed or stopped`). `core/paths.py`: `reaction_path(job)` returns a view with `points` (each `PathPoint` with `index`, `label`, `atoms`) or None.

**The bond rule on fixtures** (D3; prototype). `F = [("C", 0.0, 0.0, 0.0), ("C", 1.5, 0.0, 0.0), ("H", -0.5, 1.0, 0.0), ("O", 2.0, 1.0, 0.0), ("H", 2.0, 1.0, 1.0)]`: `bonded_pairs(F) == [(0, 1), (0, 2), (1, 3), (3, 4)]` (C1–H4 is 1.50 Å: not bonded to hydrogen). `HH = [("H", 0.0, 0.0, 0.0), ("H", 0.9, 0.0, 0.0), ("H", 5.0, 0.0, 0.0), ("H", 6.1, 0.0, 0.0), ("O", 5.0, 1.39, 0.0), ("O", 6.1, 1.41, 0.0)]`: `[(0, 1), (2, 4), (4, 5)]` (H–H 0.90 yes, 1.10 no; O–H 1.39 yes, 1.41 no; O–O 1.10 yes).

**Edits on fixtures** (prototype of `core/edit.py` exactly as AG2 specifies; coordinates rounded to 6, distances to 4):
- `rotate_bond(F, 0, 1, 90.0)`: sides {0, 2} and {1, 3, 4}, so `moved == [0, 2]`, `fixed == 1`, `mover == 0`; atom 2 becomes `("H", -0.5, 0.0, -1.0)`, every other atom unchanged; the dihedral 2-0-1-3 goes 0.0 → 90.0; `new_bonds == []`; `closest == (2.8723, 2, 3)`. With `move=1`: `moved == [1, 3, 4]`, atom 3 → `(2.0, 0.0, 1.0)`, atom 4 → `(2.0, -1.0, 1.0)`, and the dihedral 2-0-1-3 is again 90.0.
- `rotate_bond(F, 1, 3, 180.0)`: `moved == [3, 4]`, atom 4 → `(2.0, 1.0, -1.0)`, `closest == (2.4495, 0, 4)`.
- `set_dihedral(F, 2, 0, 1, 3, -60.0)`: `moved == [0, 2]`, `before == 0.0`, `after == -60.0`, atom 2 → `(-0.5, 0.5, 0.866025)`.
- `set_distance(F, 0, 1, 2.0)`: `moved == [0, 2]`, `before == 1.5`, `after == 2.0`, atom 0 → `(-0.5, 0.0, 0.0)`, atom 2 → `(-1.0, 1.0, 0.0)`.
- `G = [("C", 0.0, 0.0, 0.0), ("C", 1.5, 0.0, 0.0), ("O", 2.1, 1.2, 0.0), ("O", 2.1, -1.9, 0.0), ("H", -0.5, 0.9, 0.0), ("H", -0.5, -0.9, 0.0)]` (atom 3 is a fragment of its own): `bonded_pairs(G) == [(0, 1), (0, 4), (0, 5), (1, 2)]`. `rotate_bond(G, 0, 1, 180.0)`: `moved == [1, 2]`, atom 2 → `(2.1, -1.2, 0.0)`, `new_bonds == [(2, 3, 0.7)]`, `closest == (0.7, 2, 3)`. `rotate_bond(G, 0, 1, 90.0)`: atom 2 → `(2.1, 0.0, 1.2)`, `new_bonds == []`, `closest == (2.2472, 2, 3)`. `set_distance(G, 2, 3, 2.5)` (separate fragments): `moved == [3]`, `before == 3.1`, `after == 2.5`, atom 3 → `(2.1, -1.3, 0.0)`, `new_bonds == [(1, 3, 1.4318)]`, `closest is None` (no moved atom outside the edited pair).
- `RING = [("C", 0.0, 0.0, 0.0), ("C", 1.5, 0.0, 0.0), ("C", 0.75, 1.3, 0.0)]`.
- Refusals, with these exact messages: `rotate_bond(F, 0, 3, 10)` → `atoms 0 and 3 are not bonded: 2.24 Å apart`; `rotate_bond(F, 0, 1, 10, move=4)` → `--move takes 0 or 1, the atoms of the edited bond`; `set_dihedral(F, 3, 0, 1, 2, 10)` → `atoms 3 and 2 must lie on opposite sides of the bond 0-1`; `set_distance(F, 2, 3, 2.0)` → `atoms 2 and 3 are 2.50 Å apart in one connected fragment: moving either would bend it`; `set_distance(F, 0, 1, 0)` → `a distance must be positive`; `rotate_bond(RING, 0, 1, 10)` → `the bond 0-1 is in a ring: both sides stay connected, so neither can move alone`.

**Edits on the real tree** (`C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_b973c-xtb/ts/optts_freq_tight`, latest geometry = cycle 41, 137 atoms, QM atoms `{61, 64, 120, 121, 122, 124, 125, 126, 128, 129, 130, 131, 132, 133, 134, 135, 136}`; under D3 the structure is two fragments, of 130 and 7 atoms): rotating 133-134 by 60 moves `[134, 135]` (2 QM, 0 environment; 133's side has 5 atoms), makes no new bond, closest new contact 102 H–135 H 1.971 Å, and turns the dihedral 132-133-134-135 from -123.507 to -63.507; rotating 132-133 by 120 moves `[133, 134, 135]` and makes 5 new bonds (124 C–134 O 1.485, 124 C–135 H 1.036, 128 C–134 O 0.316, 128 C–135 H 1.041, 129 O–134 O 1.385); 124-128 is a ring bond; 129 and 55 are not bonded, 3.95 Å apart.

**Inputs in the real tree** (80 jobs). The atom numbers they use: `QMATOMS` with its list on the next line(s) in braces, with ranges (`{0:139}`, `{61 64 120 121 …}`); `ActiveAtoms {0:139} end`; `Hybrid_Hess {128 … 139} end`; `TS_Active_Atoms {128 … 136} end`; `TS_Mode {B 135 128} end`; `Constraints` with `{ C 0:60 C}`, `{ C 62 C}`; `Scan` lines `B 130 129 = 2.4, 1.5, 10` and `A 130 128 129  = 112, 75, 10`. A range `a:b` is inclusive. The prototype reader of AG4 finds NO out-of-range reference in any of the 80 inputs. The whole-system parity of D6 holds for all 74 inputs whose geometry can be read. `check` over every job (prototype): exit 1 on exactly six jobs, each with only a `geometry_file` error (their `job.xyz` is missing): `…/K+H2O2_water-relay_to_KP/geometry_b973c-xtb/neb_anion`, `…/K+H2O2_water-relay_to_KP/geometry_r2scan3c-xtb/rc/attempt1_no_recalc`, `…/rc/attempt1_qm12`, `…/geometry_r2scan3c-xtb/ts/attempt4_frozen_hybrid`, `…/K+peroxide+BnOH_bridged_UNIDENTIFIED/geometry_r2scan3c-xtb/ts/attempt1_stalled`, `…/ts/attempt2_bad_hessian`; one `qm_boundary_h` warning in the whole tree, on `…/geometry_b973c-xtb/pc_relax` (QM atom 135 H is 0.91 Å from environment atom 19 O); no `parity`, `index_range`, `nprocs` or `memory` finding anywhere.

**Connectivity on the real tree** (prototype of AG5 under D3; first printed geometry → latest unless a point is named). 9 of the 56 jobs that printed geometries report a change. `…/geometry_b973c-xtb/ts/scan` (failed): formed 58 C–129 O 2.39 → 1.53, broken 132 H–133 O 1.13 → 2.96. `…/geometry_r2scan3c-xtb/pc` (finished): formed 58 C–129 O 2.41 → 1.52, nothing broken. `…/geometry_b973c-xtb/irc_bridge --point 50`: IRC TS → IRC forward 49, broken 128 C–133 O 1.61 → 1.92 and 131 O–132 H 1.34 → 1.56; `--point -50`: broken 132 H–133 O 1.11 → 2.02. `…/geometry_b973c-xtb/ts/optts_freq_tight`: no change.

**Connectivity fixture** (parsed through the real `JobState`): with `block(atoms) = "CARTESIAN COORDINATES (ANGSTROEM)\n---------------------------------\n" + "".join(f"  {el}      {x:.6f}    {y:.6f}    {z:.6f}\n" for el, x, y, z in atoms)` and `G2` = `G` with atom 2 at `(2.1, 2.0, 0.0)` and atom 3 at `(2.9, 2.6, 0.0)`, the output `"\n".join([_CYCLE.format(n=1), block(G), "FINAL SINGLE POINT ENERGY     -80.000000000", _CYCLE.format(n=2), block(G2), "FINAL SINGLE POINT ENERGY     -80.100000000", _DONE]) + "\n"` gives points for cycles 1 and 2 with 6 atoms each; from cycle 1 to cycle 2, formed (2, 3) 3.100 → 1.000 and broken (1, 2) 1.342 → 2.088 (rounded to 3); 4 bonds before and after. With cycle 2's block replaced by `block(G[:2])` the latest geometry has 2 atoms. A cycle that prints `block(G)` and then `block(G[:2])` keeps the 6-atom block as its point.

**Test helpers.** `test_orcamon.py`: `check`, `_orcamon(argv)`, `_json(argv)`, `_CYCLE` (line 115), `_DONE` (line 121), `_coords`, `_water_job(root, name="water")` (line 187), `_OPT_FREQ`, `_irc_tree` (imported from `test_monitor` inside tests). `_irc_tree`'s path points are -2…3 with labels `IRC backward 1` … `IRC forward 2`, `0` = `IRC TS`; every frame is the water `O (0, 0, z)`, `H (0.96, 0, z)`, `H (-0.24, 0.93, z)`, so its bonds never change (2 bonds: both O–H). `test_monitor.py`: `check`, `_cycle(state, n)`, `_geometry(state, n_atoms, x=0.0)` (feeds `n_atoms` hydrogens at `x + i` and a closing blank line), `JobState`. Counts: `test_monitor.py` 246, `test_orcamon.py` 263.

## Conventions

- Match the surrounding code: one-line comments where the reason is not obvious; docstrings say WHY.
- New tests go in the file the task names, use `check(...)` for every assertion, build trees under `tempfile.TemporaryDirectory()`, and are added to that file's `__main__` list right after the test the task names. CLI calls in new tests pass `["--root", str(root), "--liveness", "mtime", "--no-cache", …]`.
- Expected values are the hand-derived ones in § Facts, never values copied from a run. Coordinates are compared after `round(…, 6)`, distances after `round(…, 4)` unless a task says otherwise.
- A new command's README entry is a `### \`orcamon <name>\`` section in `## Command reference` holding `COLUMNS=90 .venv/bin/orcamon <name> --help` in a fenced block like its neighbours; whenever a command is added, the reference's first block is replaced with `COLUMNS=90 .venv/bin/orcamon --help`. Every task that adds a command or edits `SKILL.md` ends with `cd /home/madsr2d2/masterThesis && .venv/bin/orcamon skill install --project /home/madsr2d2/masterThesis` (with `--force` only if it refuses), and mentions the new command in `SKILL.md` as `` `orcamon <name>` `` in prose (`test_the_skill_names_every_command` requires it).
- No network access.

---

# Phase 1 — the edit core

### AG1 — The edit bond rule and the element table

**Why:** edits and bond changes need a bond rule that does not count hydrogen bonds (D3), and `check` needs atomic numbers to count electrons (D6).
**Depends on:** none.
**Files:** `tools/orcamon/src/orcamon/core/geometry.py`, `tools/orcamon/src/orcamon/core/edit.py` (new), `test_monitor.py`.

**Do:**
1. `core/geometry.py`: add `ATOMIC_NUMBERS`, a dict from element symbol to atomic number for H (1) through Rn (86), written out from the periodic table in order (`"H": 1, "He": 2, "Li": 3, …, "Rn": 86`), with a comment that `check` counts electrons with it.
2. Create `core/edit.py` (standard library and `.geometry` only) with a module docstring: what it does (moves one side of a molecule rigidly, never in place: every function returns new atom lists), D2, D3 and D4 in prose. Constants `X_H_BOND_CUTOFF = 1.4` and `H_H_BOND_CUTOFF = 1.0`, each commented with why (D3).
3. `bonded_pairs(atoms) -> list[tuple[int, int]]`: the pairs of `geometry.bonds(atoms)` minus any pair with a hydrogen that is `>= X_H_BOND_CUTOFF` apart, plus every hydrogen–hydrogen pair closer than `H_H_BOND_CUTOFF` (found by binning the hydrogens into 1 Å cells, `math.floor` of each coordinate, and checking the 27 neighbouring cells), each pair `(a, b)` with `a < b`, the whole list sorted. Docstring: D3, and that the renderer's `bonds()` is unchanged on purpose.
4. `test_monitor.py`: `test_the_edit_bond_rule()`, registered right after `test_a_followed_path_point_is_addressed`. Exactly 3 checks: `bonded_pairs(F)` and `bonded_pairs(HH)` equal § Facts "The bond rule on fixtures"; and `ATOMIC_NUMBERS["H"] == 1`, `["C"] == 6`, `["O"] == 8`, `["K"] == 19`, `["Rn"] == 86`, `len(ATOMIC_NUMBERS) == 86`. Define `F` and `HH` inside the test.
5. Run the first three verification commands.

**Acceptance:**
- `test_the_edit_bond_rule` passes its 3 checks; every renderer test passes unchanged; `git diff` shows `bonds()` and `BOND_CUTOFF` untouched.
- `test_monitor.py`: 249 `pass`; `test_orcamon.py`: 263 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

### AG2 — Rotate a bond, set a dihedral, set a distance; list bond changes

**Why:** the person's "rotate bond 133-134 by 60" must become one rigid, checked motion of one side (D1, D2, D4), and "what did the run do to the bonds" must use the same rule (D3).
**Depends on:** AG1.
**Files:** `tools/orcamon/src/orcamon/core/edit.py`, `test_monitor.py`.

**Do:**
1. `class EditError(ValueError)` (docstring: the edit is impossible or ill-posed; the message says why). `@dataclass class EditResult` with fields `atoms: list` (the new geometry), `moved: list` (sorted atom numbers), `fixed: int`, `mover: int` (the edited bond's atom that stays and the one that moves), `fixed_side: int` (the number of atoms on the fixed atom's side), `new_bonds: list` (`[(a, b, distance)]`, `a < b`), `closest: tuple | None` (`(distance, a, b)`, `a < b`), `before: float | None = None`, `after: float | None = None`.
2. Private helpers: `_component(n, pairs, start, cut)` — the atoms reachable from `start` through `pairs` without using the pair `cut`; `_sides(atoms, pairs, i, j, need_bond)` — with `key = (min(i, j), max(i, j))` and `d` their distance: when `need_bond` and `key` is not in `pairs`, raise `EditError(f"atoms {i} and {j} are not bonded: {d:.2f} Å apart")`; `side_j = _component(n, pairs, j, key)`; when `i in side_j`, raise `EditError(f"the bond {i}-{j} is in a ring: both sides stay connected, so neither can move alone")` if `key` is a bond, else `EditError(f"atoms {i} and {j} are {d:.2f} Å apart in one connected fragment: moving either would bend it")`; return `(_component(n, pairs, i, key), side_j)`. `_choose(side_i, side_j, i, j, move)`: `move` None → `(i, j, side_j)` when `len(side_j) <= len(side_i)` else `(j, i, side_i)`; `move == j` → `(i, j, side_j)`; `move == i` → `(j, i, side_i)`; anything else → `EditError(f"--move takes {i} or {j}, the atoms of the edited bond")`. Returns `(fixed, mover, moved)`.
3. `_rotate(atoms, moved, fixed, mover, degrees)`: Rodrigues' rotation of each moved atom about the unit axis from `fixed` to `mover`, through `fixed`: `v' = v cos θ + (u × v) sin θ + u (u · v)(1 - cos θ)` with `v` the atom's position minus `fixed`'s. `_contacts(before_pairs, after, moved, edited)` with `edited = (min(i, j), max(i, j))`: `new_bonds` = the pairs of `bonded_pairs(after)` that are not in `before_pairs`, are not `edited`, and have exactly one atom in `moved`, each with its distance; `closest` = the smallest distance over `a` in `moved`, `b` not in `moved`, neither in `edited`, and `(min, max)` not in `before_pairs` — returned as `(distance, min(a, b), max(a, b))`, or None when there is no such pair.
4. Public functions, each computing `pairs = bonded_pairs(atoms)` once and returning an `EditResult`: `rotate_bond(atoms, i, j, degrees, move=None)` (`_sides(…, need_bond=True)`, `_choose`, `_rotate`, `_contacts`; `fixed_side` = the size of the fixed atom's side; `before`/`after` None); `set_dihedral(atoms, a, b, c, d, target, move=None)` (`_sides(atoms, pairs, b, c, True)`; `a` not in b's side or `d` not in c's side → `EditError(f"atoms {a} and {d} must lie on opposite sides of the bond {b}-{c}")`; `current = measure(atoms, [a, b, c, d])[1]`; `delta = ((target - current + 180.0) % 360.0) - 180.0`; the result of `rotate_bond(atoms, b, c, delta, move)` with `before = current`, `after = measure(result.atoms, [a, b, c, d])[1]`); `set_distance(atoms, i, j, target, move=None)` (`target <= 0` → `EditError("a distance must be positive")`; `_sides(…, need_bond=False)`; `_choose`; every moved atom shifted by `(target - current)` times the unit vector from `fixed` to `mover`; `_contacts`; `before = current`, `after` = the new i–j distance). And `bond_changes(before_atoms, after_atoms) -> tuple[list, list, int]`: `formed` = pairs in `bonded_pairs(after)` not in `bonded_pairs(before)`, `broken` = the reverse, each as `(a, b, distance_before, distance_after)`, both sorted, and the number of bonds after; it raises `EditError` when the two lists differ in length. Docstrings: D4's meaning of the angle, that nothing is changed in place, and that `bond_changes` counts the D3 rule.
5. `test_monitor.py`: `test_edits_move_one_side_and_report_contacts()`, registered right after `test_the_edit_bond_rule`, with `F`, `G` and `RING` from § Facts defined inside it. Exactly 11 checks, values from § Facts "Edits on fixtures":
   - `rotate_bond(F, 0, 1, 90.0)`: `moved == [0, 2]`, `fixed == 1`, `mover == 0`, `fixed_side == 3`, atom 2 rounds to `("H", -0.5, 0.0, -1.0)` and atoms 0, 1, 3, 4 equal `F`'s;
   - its dihedral 2-0-1-3 is 0.0 before and 90.0 after, and with `move=1` the result has `moved == [1, 3, 4]`, atoms 3 and 4 at `(2.0, 0.0, 1.0)` and `(2.0, -1.0, 1.0)`, and the same dihedral 90.0;
   - its `new_bonds == []` and `closest` rounds to `(2.8723, 2, 3)`;
   - `set_dihedral(F, 2, 0, 1, 3, -60.0)`: `moved == [0, 2]`, `before` 0.0, `after` -60.0, atom 2 at `(-0.5, 0.5, 0.866025)`;
   - `set_distance(F, 0, 1, 2.0)`: `moved == [0, 2]`, `before` 1.5, `after` 2.0, atoms 0 and 2 at `(-0.5, 0.0, 0.0)` and `(-1.0, 1.0, 0.0)`;
   - `rotate_bond(F, 1, 3, 180.0)`: `moved == [3, 4]`, atom 4 at `(2.0, 1.0, -1.0)`, `closest` `(2.4495, 0, 4)`;
   - `rotate_bond(G, 0, 1, 180.0)`: `moved == [1, 2]`, atom 2 at `(2.1, -1.2, 0.0)`, `new_bonds` rounds to `[(2, 3, 0.7)]`, `closest` to `(0.7, 2, 3)`;
   - `rotate_bond(G, 0, 1, 90.0)`: atom 2 at `(2.1, 0.0, 1.2)`, `new_bonds == []`, `closest` `(2.2472, 2, 3)`;
   - `set_distance(G, 2, 3, 2.5)`: `moved == [3]`, `before` 3.1, `after` 2.5, atom 3 at `(2.1, -1.3, 0.0)`, `new_bonds` `[(1, 3, 1.4318)]`, `closest is None`;
   - the six refusals of § Facts raise `EditError` with exactly those messages;
   - `bond_changes(G, rotate_bond(G, 0, 1, 180.0).atoms)` returns `formed` rounding to `[(2, 3, 3.1, 0.7)]`, `broken == []` and 5 bonds after (G's 4 plus the new one).
6. Run the first three verification commands, and on the real tree (read-only; `XDG_CACHE_HOME=/tmp/orcamon_act/xdg`): load `optts_freq_tight` as § Facts names it with `Job` + `refresh(Liveness(False, "process"))` and print `rotate_bond(atoms, 133, 134, 60.0)`'s `moved`, `fixed_side`, `new_bonds`, rounded `closest` and the dihedral 132-133-134-135 after, and `len(rotate_bond(atoms, 132, 133, 120.0).new_bonds)`.

**Acceptance:**
- The real tree prints `[134, 135]`, `5`, `[]`, `(1.971, 102, 135)`, `-63.507` and `5` (§ Facts "Edits on the real tree").
- `test_edits_move_one_side_and_report_contacts` passes its 11 checks.
- `test_monitor.py`: 260 `pass`; `test_orcamon.py`: 263 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

# Phase 2 — the commands

### AG3 — `orcamon edit` writes a checked geometry file

**Why:** the agent's half of "rotate this bond": one command that turns the person's words into a new XYZ file, refusing what would mislead (D1, D5).
**Depends on:** AG2.
**Files:** `tools/orcamon/src/orcamon/cli/__init__.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/src/orcamon/SKILL.md`, `.claude/skills/orcamon/SKILL.md` (regenerated), `tools/orcamon/README.md`, `test_orcamon.py`.

**Do:**
1. `cli/__init__.py`, after `measure`: `add("edit", "write a changed copy of a geometry: rotate a bond, set a dihedral or a distance", "orcamon edit opt/ts --rotate-bond 12 15 60 -o rotated.xyz", commands.cmd_edit)`; `_job(p)`; a required mutually exclusive group with `--rotate-bond` (`nargs=3`, `metavar=("I", "J", "DEG")`, help `"turn the bond I-J: every dihedral across it grows by DEG degrees"`), `--set-dihedral` (`nargs=5`, `metavar=("A", "B", "C", "D", "DEG")`, help `"turn the bond B-C until the dihedral A-B-C-D is DEG degrees"`), `--set-distance` (`nargs=3`, `metavar=("I", "J", "ANGSTROM")`, help `"move one side along I-J until they are ANGSTROM apart"`); `-o/--output` (required, `metavar="FILE"`, help `"the XYZ file to write; an existing file is refused"`); `--move` (`type=int`, `metavar="ATOM"`, help `"move the side of this atom of the edited bond (default: the smaller side)"`); `--allow-new-bonds` (`store_true`, help `"write the file even when the edit makes a new bond"`); `--cycle`, `--step`, `--point` exactly as `geom` declares them; `_output(p, lines=False)`.
2. `commands.py`, `@_with_job def cmd_edit(args, job)`: refuse an existing `args.output` first (`UsageError(f"{args.output} exists: orcamon edit never overwrites a file")`); `point = _requested_point(args, job)`; parse the chosen flag's values — atom numbers with `int()`, the angle or distance with `float()`, a bad value → `UsageError(f"not a number: {value!r}")`; check the atom numbers against `len(point.atoms)` with `geometry._check_in_range` (its `ValueError` → `UsageError`); call `edit.rotate_bond` / `set_dihedral` / `set_distance` with `move=args.move` (`EditError` → `UsageError`). `n = len(point.atoms)`; `qm = job.state.qm_atom_indices`, used for the split only when `qm and 0 < len(qm) < n`.
3. Text, one line each, in order: the action — `f"rotate bond {i}-{j} by {deg:+.1f}°"`, `f"set dihedral {a}-{b}-{c}-{d} to {after:.1f}° (was {before:.1f}°)"` or `f"set distance {i}-{j} to {after:.3f} Å (was {before:.3f} Å)"`; `f"moved {k} of {n} atoms, the {mover} side{split}; --move {fixed} moves the other {fixed_side}"` with `split = f" ({q} QM, {e} environment)"` on a multilayer job, else `""`; `f"closest new contact: {a} {el_a}-{b} {el_b} {d:.3f} Å"` or `closest new contact: none`; then one `f"new bond: {a} {el_a}-{b} {el_b} {d:.3f} Å"` per new bond. When there are new bonds and no `--allow-new-bonds`, print those lines and then `return _error(f"refused: the edit makes {k} new bond{'s' if k != 1 else ''}; nothing written (--allow-new-bonds writes it anyway)")`. Otherwise write the file — line 1 the atom count, line 2 `f"{job.label} · {describe_point(point) or 'latest'} · {action line}"`, then `f"{el:<2} {x:15.8f} {y:15.8f} {z:15.8f}"` per atom, as `geom` prints — and end with `f"wrote {args.output} ({n} atoms)"`. Return `EXIT_OK` (a written file is the result, whatever the job's status).
4. JSON instead of text: `{"job": job.label, "point": _point_doc(point), "edit": {"kind": "rotate-bond" | "set-dihedral" | "set-distance", "atoms": [...], "value": float}, "before": …, "after": …, "moved": [...], "moved_qm": q or None, "moved_environment": e or None, "fixed_side": int, "closest": {"atoms": [a, b], "elements": [..], "distance": round(d, 4)} or None, "new_bonds": [{"atoms": [a, b], "elements": [..], "distance": round(d, 4)}, …], "written": path or None}`; a refused edit prints the document with `"written": None`, the same `orcamon: refused: …` line on stderr, and exits 2.
5. `SKILL.md`, "Looking closer" is for reading; add a new section `## Changing a geometry` before "Working beside a person" with: "`orcamon edit JOB --rotate-bond I J DEG -o new.xyz` writes a copy of the job's geometry with one side of the bond I-J turned so that every dihedral across it grows by DEG degrees; `--set-dihedral A B C D DEG` and `--set-distance I J ANGSTROM` set one value instead. The smaller side moves (`--move ATOM` picks the other), the output names how many atoms moved — QM and environment on a multilayer job — and the closest new contact, and an existing file is never overwritten. It refuses atoms that are not bonded, a bond in a ring, and an edit that makes a new bond. Do not pass `--allow-new-bonds` unless the person asked for that contact. A bond here is under 1.7 Å between heavy atoms, 1.4 Å to a hydrogen and 1.0 Å between hydrogens, so a hydrogen bond, a long transition-state bond or an ion's contact is not one." README: a section `## Changing a geometry` after "Reaction paths" with the same content in the README's voice, and the `edit` reference section plus the refreshed first block (§ Conventions). Reinstall the skill.
6. `test_orcamon.py`: `test_edit_writes_a_checked_geometry()`, registered right after `test_the_tui_picks_atoms_by_click`. A temp root with jobs `f` and `g`, each a draft whose `job.inp` is `"! B97-3c Opt\n* xyz 0 1\n" + "".join(f"{el} {x} {y} {z}\n" for el, x, y, z in ATOMS) + "*\n"` with `F` and `G` of § Facts, and a job `ml` whose `job.inp` is `"! QM/XTB B97-3c Opt\n%qmmm\n  QMAtoms {0 1} end\n  Charge_Total 0\n  Mult_Total 1\nend\n* xyz 0 1\n*\n"` and whose `job.out` is `"QM1 Subsystem                    ...  0 1\n" + block(G) + "\n" + _DONE + "\n"` (`block` as in § Facts). Files are written to `tmp / "out"`. Exactly 8 checks:
   - `edit f --rotate-bond 0 1 90 -o <out>/r.xyz`: exit 0 and the stdout lines are exactly `rotate bond 0-1 by +90.0°`, `moved 2 of 5 atoms, the 0 side; --move 1 moves the other 3`, `closest new contact: 2 H-3 O 2.872 Å`, `wrote <out>/r.xyz (5 atoms)`;
   - `read_xyz(<out>/r.xyz)` has 5 atoms with atom 2 rounding to `("H", -0.5, 0.0, -1.0)`, and the file's second line is `f · input geometry · rotate bond 0-1 by +90.0°`;
   - `edit f --set-dihedral 2 0 1 3 -60 -o <out>/d.xyz --json`: `doc["after"] == -60.0` (rounded to 4), `doc["moved"] == [0, 2]`, `doc["written"] == str(<out>/d.xyz)`, `doc["new_bonds"] == []`;
   - `edit g --rotate-bond 0 1 180 -o <out>/x.xyz`: exit 2, stdout contains `new bond: 2 O-3 O 0.700 Å`, stderr contains `refused: the edit makes 1 new bond; nothing written (--allow-new-bonds writes it anyway)`, and `<out>/x.xyz` does not exist;
   - the same with `--allow-new-bonds`: exit 0 and the file exists;
   - with `<out>/r.xyz` existing, `edit f --rotate-bond 0 1 10 -o <out>/r.xyz`: exit 2, stderr contains `exists: orcamon edit never overwrites a file`, and the file's content is unchanged;
   - `edit f --rotate-bond 0 3 10 -o <out>/n.xyz` exits 2 with `atoms 0 and 3 are not bonded: 2.24 Å apart`, and `edit f --rotate-bond 0 1 abc -o <out>/n.xyz` exits 2 with `not a number: 'abc'`;
   - `edit ml --rotate-bond 0 1 90 -o <out>/m.xyz`: its second line is `moved 2 of 6 atoms, the 1 side (1 QM, 1 environment); --move 0 moves the other 3`.
7. Run the first three verification commands, and on the real tree (`cd /home/madsr2d2/masterThesis`, the two environment variables of § Never, `mkdir -p /tmp/orcamon_act/run`): `orcamon --root computational edit C8_perhydrate_trap/K+H2O2_water-relay_to_KP/geometry_b973c-xtb/ts/optts_freq_tight --rotate-bond 133 134 60 -o /tmp/orcamon_act/run/r60.xyz`, then the same with `--rotate-bond 132 133 120 -o /tmp/orcamon_act/run/r120.xyz`, `--rotate-bond 124 128 10 -o /tmp/orcamon_act/run/ring.xyz` and `--rotate-bond 129 55 10 -o /tmp/orcamon_act/run/nb.xyz`.

**Acceptance:**
- The first real-tree command exits 0 and prints `rotate bond 133-134 by +60.0°`, `moved 2 of 137 atoms, the 134 side (2 QM, 0 environment); --move 133 moves the other 5`, `closest new contact: 102 H-135 H 1.971 Å` and the `wrote` line; the second exits 2 after five `new bond:` lines (the first `new bond: 124 C-134 O 1.485 Å`) with `refused: the edit makes 5 new bonds`; the third and fourth exit 2 with `the bond 124-128 is in a ring` and `atoms 129 and 55 are not bonded: 3.95 Å apart`; only `r60.xyz` exists afterwards.
- `test_edit_writes_a_checked_geometry` passes its 8 checks; `test_the_skill_names_every_command`, `test_the_skill_reference_matches_the_parser` and `test_this_repository_has_the_current_skill` pass.
- `test_orcamon.py`: 271 `pass`; `test_monitor.py`: 260 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

### AG4 — Every atom number an input uses, and `orcamon check`

**Why:** an input is cheapest to fix before it runs: an atom number past the end of the geometry, an electron count the multiplicity cannot have, a missing geometry file or a machine too small should be found by the agent, not by ORCA hours into a queue (D6).
**Depends on:** AG1.
**Files:** `tools/orcamon/src/orcamon/core/orca_input.py`, `tools/orcamon/src/orcamon/core/inputcheck.py` (new), `tools/orcamon/src/orcamon/cli/__init__.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/src/orcamon/SKILL.md`, `.claude/skills/orcamon/SKILL.md` (regenerated), `tools/orcamon/README.md`, `test_monitor.py`, `test_orcamon.py`.

**Do:**
1. `core/orca_input.py`: `atom_references(text) -> list[tuple[str, list[int]]]` on the comment-stripped text (`_COMMENT_RE`): for every brace group `\{([^{}]*)\}`, split its content on whitespace and commas; when the first token is an integer or a range `a:b`, every token must be one (else the group is skipped) and the group is a list of atoms; when the first token is `B`, `A`, `D` or `C` (case-insensitive), the next 2, 3, 4 or 1 tokens are atoms (integers or ranges; else skipped) and the rest are values; any other group is skipped. A range `a:b` expands to `a…b` inclusive. The reference's name is the group as written, `f"{{{content.strip()}}}"`. Also every line matching `^\s*([BAD])\s+(\d+(?:\s+\d+)*)\s*=` is a scan coordinate, named `f"Scan {letter} {numbers}"` (e.g. `Scan B 0 5`). And `qm_atoms(text) -> set | None`: the atoms of the first `qmatoms\s*\{([^}]*)\}` (case-insensitive, across lines), expanded, or None.
2. `core/inputcheck.py` (standard library and `orcamon.core` only): `@dataclass class Finding(code: str, message: str)`, `@dataclass class InputCheck(atoms: int | None, electrons: int | None, charge: int | None, mult: int | None, errors: list, warnings: list)`, and `check_input(job) -> InputCheck`, reading `job.inp_path`'s text and `job.input`. The geometry is `input.coords_atoms` (an inline `* xyz`), else `read_xyz(job.state.path / input.coords_file)` for an `xyzfile`; no `*` geometry line at all → error `no_geometry` `the input has no * xyz or * xyzfile line`; an unreadable file → error `geometry_file` `f"{coords_file}: missing or not an XYZ file"`; any other form (`int`, `gzmt`, `pdbfile`, `%coords`) → warning `geometry_not_read` `f"a {kind} geometry is not read: the atom checks are skipped"`. With atoms: an element missing from `ATOMIC_NUMBERS` → error `unknown_element` `f"atom {i}: unknown element {el!r}"`; the whole system's `(charge, mult)` is `input.layers["total"]` on a multilayer job (when either half is None: warning `total_unset` `Charge_Total and Mult_Total are not both set: the whole system's electrons are not checked`) and `(input.charge, input.mult)` otherwise; `electrons = sum(Z) - charge`; `mult < 1` or `mult - 1 > electrons` → error `parity` `f"multiplicity {mult} is impossible with {electrons} electrons"`; `(electrons - (mult - 1)) % 2 != 0` → error `parity` `f"{electrons} electrons cannot have multiplicity {mult}: it needs an {'odd' if mult % 2 == 0 else 'even'} number"`; for each reference of `atom_references` with an atom `>= len(atoms)`, error `index_range` `f"atom {first such} in `{name}` is out of range: the geometry has {n} atoms (0-{n - 1})"` (one finding per reference).
3. Warnings, with atoms: on a multilayer job with `qm_atoms`, every pair of `edit.bonded_pairs(atoms)` with one atom in the QM set, one outside, and a hydrogen in it → `qm_boundary_h` `f"QM atom {q} {el} is {d:.2f} Å from environment atom {e} {el}: a bond to hydrogen crosses the QM boundary"`. Without atoms too: `input.nprocs > os.cpu_count()` → `nprocs` `f"nprocs {p} is more than this machine's {cpus} CPUs"`; when `maxcore_mb` is set, `nprocs * maxcore_mb` (nprocs 1 when unset) more than the machine's memory from `os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") // 2**20` (skipped when `os.sysconf` raises) → `memory` `f"{p} process{'es' if p != 1 else ''} x {mc} MB maxcore = {p * mc} MB is more than this machine's {mem} MB"`; `job.out_path.exists()` → `output_exists` `f"{stem}.out exists: running this job again overwrites it"`.
4. CLI: register `add("check", "what in a job's input would fail or mislead before it runs", "orcamon check opt/ts", commands.cmd_check)` after `edit`, `_job(p)`, `_output(p)`. `cmd_check`: text is one `f"error: {message}"` line per error, then one `f"warning: {message}"` per warning, then `f"{E} error(s), {W} warning(s) · {n} atoms, {electrons} electrons, charge {charge}, multiplicity {mult}"` (the part after `·` only when atoms and electrons are known; else `· atoms not read`); JSON `{"job", "atoms", "electrons", "charge", "mult", "errors": [{"code", "message"}], "warnings": [...]}`. Exit 1 when there is an error, else 0. In `EXIT_CODES`, exit 1's meaning becomes `the command worked and the job (for ls: any listed job) is failed or stopped, or (check) the input has an error`; update the docstring's exit sentence and the README table the same way.
5. `SKILL.md`, in "Changing a geometry": "`orcamon check JOB` reads a job's input before it runs: the geometry (inline or the `*xyzfile`), the whole system's electron count against the multiplicity, every atom number the input uses (QM atoms, active atoms, constraints, scans, TS modes) against the geometry, the machine's CPUs and memory, and an output that a run would overwrite. Exit 1 means an error: fix it before running. It does not check the QM region's electrons — they depend on the link atoms ORCA creates." README: the same in "Changing a geometry", the reference section, the refreshed first block. Reinstall the skill.
6. Tests. `test_monitor.py`: `test_the_input_names_every_atom_it_uses()`, registered right after `test_edits_move_one_side_and_report_contacts`, on `text = "! QM/XTB B97-3c OptTS\n# {99} in a comment\n%qmmm\n  QMATOMS\n    {0:2 5}\n  END\n  ActiveAtoms {0:9} end\nend\n%geom\n  TS_Mode {B 4 5} end\n  Constraints\n    { C 0:3 C }\n    { B 1 2 1.5 C }\n    { D 0 1 2 3 180 C }\n  end\n  Scan\n    B 1 2 = 1.0, 2.0, 5\n  end\nend\n* xyzfile 0 1 x.xyz\n"`. Exactly 3 checks: `atom_references(text) == [("{0:2 5}", [0, 1, 2, 5]), ("{0:9}", list(range(10))), ("{B 4 5}", [4, 5]), ("{C 0:3 C}", [0, 1, 2, 3]), ("{B 1 2 1.5 C}", [1, 2]), ("{D 0 1 2 3 180 C}", [0, 1, 2, 3]), ("Scan B 1 2", [1, 2])]`; `qm_atoms(text) == {0, 1, 2, 5}`; `atom_references("* xyz 0 1\n*\n") == []`. `test_orcamon.py`: `test_check_finds_what_would_fail()`, registered right after `test_edit_writes_a_checked_geometry`, on drafts whose geometry is the water `O 0.0 0.0 0.0`, `H 0.96 0.0 0.0`, `H -0.24 0.93 0.0` inline. Exactly 8 checks:
   - `! B97-3c Opt` with `* xyz 0 1`: exit 0, last line `0 error(s), 0 warning(s) · 3 atoms, 10 electrons, charge 0, multiplicity 1`;
   - with `* xyz 0 2`: exit 1, stdout contains `error: 10 electrons cannot have multiplicity 2: it needs an odd number`;
   - with `%geom Scan B 0 5 = 1.0, 2.0, 3 end end` (each on its own line): exit 1, stdout contains ``error: atom 5 in `Scan B 0 5` is out of range: the geometry has 3 atoms (0-2)``;
   - `* xyzfile 0 1 missing.xyz` and no atoms: exit 1, `error: missing.xyz: missing or not an XYZ file`;
   - `%pal nprocs 100000 end`: exit 0, a `warning:` line containing `nprocs 100000 is more than this machine's`;
   - the clean water draft with an empty `job.out` beside it: the `--json` document's `warnings` contains `{"code": "output_exists", "message": "job.out exists: running this job again overwrites it"}`;
   - `! QM/XTB B97-3c Opt\n%qmmm\n  QMAtoms {0 1} end\n  Charge_Total 0\n  Mult_Total 1\nend\n* xyz 0 1` + water: exit 0, stdout contains `warning: QM atom 0 O is 0.96 Å from environment atom 2 H: a bond to hydrogen crosses the QM boundary`;
   - the `* xyz 0 2` draft with `--json`: `doc["errors"][0]["code"] == "parity"` and `doc["electrons"] == 10`.
7. Run the first three verification commands, and on the real tree: `orcamon --root computational ls --json` gives the labels (`jobs[].label`); run `orcamon --root computational check LABEL --json` for each.

**Acceptance:**
- Over the real tree, exactly the six jobs of § Facts "Inputs in the real tree" exit 1, each with only `geometry_file` errors; the only `qm_boundary_h` warning is on `pc_relax`, naming atoms 135 H and 19 O at 0.91 Å; no job has a `parity`, `index_range`, `nprocs` or `memory` finding. Report the number of `output_exists` warnings.
- `test_the_input_names_every_atom_it_uses` passes its 3 checks and `test_check_finds_what_would_fail` its 8; `test_the_input_names_charge_and_multiplicity` passes unchanged.
- `test_monitor.py`: 263 `pass`; `test_orcamon.py`: 279 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

### AG5 — Keep the first geometry, and `orcamon connectivity`

**Why:** "did the run do something to the bonds it should not have" — the costliest past failure in this kind of work is a bond formed or broken that nobody noticed — and the baseline must be what ORCA started from, not the input file it overwrites (D7).
**Depends on:** AG2.
**Files:** `tools/orcamon/src/orcamon/core/parser.py`, `tools/orcamon/src/orcamon/cli/__init__.py`, `tools/orcamon/src/orcamon/cli/commands.py`, `tools/orcamon/src/orcamon/SKILL.md`, `.claude/skills/orcamon/SKILL.md` (regenerated), `tools/orcamon/README.md`, `test_monitor.py`, `test_orcamon.py`.

**Do:**
1. `core/parser.py`: `JobState.start_atoms: list = field(default_factory=list)` after `atoms`, commented: the first point's geometry — what ORCA started from — kept however long the run, because `points` drops its oldest after `HISTORY_LEN`, and ORCA overwrites the input's `<basename>.xyz`. A private `_start_point: object = None`. In the coordinate-block close (§ Facts "The code"), right after `point = self.current_point()`: `if self._start_point is None: self._start_point = point`; and inside the `if len(self._pending_atoms) > len(point.atoms):` branch, `if point is self._start_point: self.start_atoms = self._pending_atoms`.
2. `test_monitor.py`: `test_the_first_geometry_is_kept()`, registered right after `test_the_input_names_every_atom_it_uses`. Exactly 3 checks on fresh `JobState(path=Path("/nonexistent"))`:
   - for `c` in 1…502, `_cycle(state, c)` and `_geometry(state, 3, x=float(c))`: `state.start_atoms` is the three hydrogens at x = 1.0, 2.0, 3.0 (`[("H", 1.0, 0.0, 0.0), ("H", 2.0, 0.0, 0.0), ("H", 3.0, 0.0, 0.0)]`) and `state.points[0].cycle != 1`;
   - `_cycle(state, 1)`, `_geometry(state, 5)`, `_geometry(state, 2)`: `len(state.start_atoms) == 5` (the largest block of the first cycle);
   - `test_a_cached_state_equals_a_fresh_parse`'s assumptions hold: a `pickle.loads(pickle.dumps(state))` copy of the 502-cycle state has the same `start_atoms`.
3. CLI: `add("connectivity", "the bonds a run formed and broke", "orcamon connectivity opt/ts", commands.cmd_connectivity)` after `check`, `_job(p)`, `--cycle`, `--step`, `--point` as `geom` declares them, `_output(p)`. `cmd_connectivity`: with `--point N` (on a path): `before` = the path point with index 0 (`IRC TS` / `NEB image 0`) and `after = _path_point(job, N)`, named by their labels; otherwise `before = job.state.start_atoms` named `first printed geometry` (empty → `UsageError("this output has printed no geometry yet: nothing to compare")`) and `after = _requested_point(args, job)` named `describe_point(after) or "latest"`. Different atom counts → `UsageError(f"the first printed geometry has {len(before)} atoms and {name} has {len(after)}: ORCA printed different parts of the system, so they cannot be compared")` (on a path: the two labels). `formed, broken, n = edit.bond_changes(before, after)`.
4. Text: line 1 `f"{from} -> {to}: {F} bond{'s' if F != 1 else ''} formed, {B} broken"`, or `f"{from} -> {to}: no bond formed or broken ({n} bonds)"` when both are empty; then `f"formed  {a} {el}-{b} {el}  {d0:.2f} -> {d1:.2f} Å"` per formed and `f"broken  {a} {el}-{b} {el}  {d0:.2f} -> {d1:.2f} Å"` per broken; last line always `bonds: under 1.7 Å between heavy atoms, 1.4 Å to a hydrogen, 1.0 Å between hydrogens`. JSON `{"job", "from", "to", "bonds": n, "formed": [{"atoms": [a, b], "elements": [..], "before": round(d0, 3), "after": round(d1, 3)}], "broken": [...]}`. Exit `_job_exit(job)`.
5. `SKILL.md`, in "Changing a geometry": "`orcamon connectivity JOB` lists the bonds a run formed and broke, from the first geometry the output printed to the latest (or `--cycle`, `--step`); on an IRC or NEB, `--point N` compares the TS (image 0) with point N. Run it after an optimization or a scan and before trusting the result: a bond nobody asked for is the costliest mistake a run can hide. Distances are given before and after, so a bond that only crossed the threshold is visible as one." README: the same in "Changing a geometry", the reference section, the refreshed first block. Reinstall the skill.
6. `test_orcamon.py`: `test_connectivity_names_bonds_made_and_broken()`, registered right after `test_check_finds_what_would_fail`. Jobs: `c` (`job.inp` `"! B97-3c Opt\n* xyz 0 1\n*\n"`, `job.out` the § Facts connectivity fixture), `short` (the same with cycle 2's block `block(G[:2])`), `draft` (`job.inp` only, inline water), and `irc` (`_irc_tree(root / "irc")`). Exactly 6 checks:
   - `connectivity c`: exit 0; line 1 `first printed geometry -> cycle 2: 1 bond formed, 1 broken`; the lines include `formed  2 O-3 O  3.10 -> 1.00 Å` and `broken  1 C-2 O  1.34 -> 2.09 Å`; the last line is the `bonds:` line;
   - `connectivity c --cycle 1`: line 1 `first printed geometry -> cycle 1: no bond formed or broken (4 bonds)`;
   - `connectivity c --json`: `formed == [{"atoms": [2, 3], "elements": ["O", "O"], "before": 3.1, "after": 1.0}]` and `broken == [{"atoms": [1, 2], "elements": ["C", "O"], "before": 1.342, "after": 2.088}]`, `bonds == 4`;
   - `connectivity irc --point 1`: line 1 `IRC TS -> IRC forward 0: no bond formed or broken (2 bonds)`;
   - `connectivity draft`: exit 2, stderr contains `this output has printed no geometry yet`;
   - `connectivity short`: exit 2, stderr contains `the first printed geometry has 6 atoms and cycle 2 has 2`.
7. Run the first three verification commands, and on the real tree: `orcamon --root computational connectivity` for `…/geometry_b973c-xtb/ts/scan`, `…/geometry_r2scan3c-xtb/pc`, `…/geometry_b973c-xtb/irc_bridge --point 50` and `…/geometry_b973c-xtb/ts/optts_freq_tight` (labels prefixed `C8_perhydrate_trap/K+H2O2_water-relay_to_KP/`).

**Acceptance:**
- `ts/scan` exits 1 (a failed job) with `formed  58 C-129 O  2.39 -> 1.53 Å` and `broken  132 H-133 O  1.13 -> 2.96 Å`; `pc` exits 0 with `formed  58 C-129 O  2.41 -> 1.52 Å` and nothing broken; `irc_bridge --point 50` prints `IRC TS -> IRC forward 49: 0 bonds formed, 2 broken` with `broken  128 C-133 O  1.61 -> 1.92 Å` and `broken  131 O-132 H  1.34 -> 1.56 Å`; `optts_freq_tight` prints `no bond formed or broken` (report its bond count).
- `test_the_first_geometry_is_kept` passes its 3 checks and `test_connectivity_names_bonds_made_and_broken` its 6; `test_a_cached_state_equals_a_fresh_parse` passes unchanged.
- `test_monitor.py`: 266 `pass`; `test_orcamon.py`: 285 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`.

### AG6 — The workflow in the skill, then every gate

**Why:** the agent needs the order the commands go in when the person asks for a change, in one place, and the whole repository must still pass.
**Depends on:** AG3, AG4, AG5.
**Files:** `tools/orcamon/src/orcamon/SKILL.md`, `.claude/skills/orcamon/SKILL.md` (regenerated), `tools/orcamon/README.md`.

**Do:**
1. `SKILL.md`, at the end of "Changing a geometry", one paragraph: "When the person asks for a change ("rotate this bond 60 degrees", "set this distance to 2.1"): run `orcamon view` for the job, the geometry and the atoms they mean; `orcamon edit` it into a new file in a new job directory; write the input there (your own project's conventions decide how); `orcamon check` it; `orcamon focus` the draft so the person sees the new structure; and run it only when they say so. After the run, `orcamon connectivity` says what it did to the bonds." README: the same paragraph in "Changing a geometry". Reinstall the skill.
2. Run all four verification commands.

**Acceptance:**
- `test_monitor.py`: 266 `pass`; `test_orcamon.py`: 285 `pass`; both `0 failure(s)`. `data/test_curve_metrics.py`: `0 failure(s)`. `run_gates.py`: `38 gates`, `0 failed`.
- `git diff a825532 -- '*test*'` across the whole plan shows only added tests, fixtures, imports and `__main__` entries; no existing assertion or expected value changed.

---

## Out of scope

1. Running and stopping jobs — `orcamon run` (local, detached, after `check`) and `orcamon stop` — and the skill's rule for them: the Act-run plan.
2. The QM region's electron count (D6), link-atom prediction, and warnings about constrained or frozen atoms that an edit moves.
3. A covalent-radius rule for heavy atoms: the 1.7 Å rule misses longer heavy-atom bonds (S–S, C–S, C–Cl, C–Br), in the renderer and here alike.
4. Rotating about two atoms that are not bonded (a forming bond), angle edits, and internal-coordinate (`* int`, `gzmt`) geometries.
5. Showing an edit in the TUI other than through a draft input and `orcamon focus`.
6. Any judgement on the real-tree findings — `pc_bridge`'s broken C120–C128, the C58–O129 bond in `ts/scan` and `pc`, `pc_relax`'s boundary hydrogen: they are the user's to read.

## Progress template

```markdown
# orcamon — the agent changes a geometry, checks an input, reads bond changes — Progress

Plan: `PLAN_ORCAMON_ACT_GEOMETRY.md`. Branch `master` from `master` (`<sha of the plan commit>`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| AG1 the edit bond rule, element table | TODO | | | |
| AG2 rotate, set dihedral, set distance, bond changes | TODO | | | |
| AG3 orcamon edit | TODO | | | |
| AG4 atom references, orcamon check | TODO | | | |
| AG5 first geometry, orcamon connectivity | TODO | | | |
| AG6 the workflow in the skill, every gate | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 246 | 263 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

## Backlog
```
