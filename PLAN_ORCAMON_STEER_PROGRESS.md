# orcamon — the person and the agent share one view — Progress

Plan: `PLAN_ORCAMON_STEER.md`. Branch `master` from `master` (`9b4abeae965da17063bd4a866c93e8058ed29fd4`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| ST1 measure and the command | DONE | 5816182 | 1 | reviewer PASS r1: 216/218 pass, 0 failures; real-tree distance/angle/dihedral values matched |
| ST2 the state files | DONE | db8340d, 69a0d63 | 2 | reviewer PASS r2: 224/218 pass, 0 failures; r1 finding (PathPoint while following) fixed in 69a0d63 |
| ST3 renderers draw the selection | DONE | 04027ae | 1 | reviewer PASS r1: 231/218 pass, 0 failures |
| ST4 s picks atoms | DONE | 56fe103 | 1 | reviewer PASS r1: 231/224 pass, 0 failures |
| ST5 the TUI publishes; orcamon view | DONE | c908f4e | 1 | reviewer PASS r1: 231/231 pass, 0 failures |
| ST6 orcamon focus | TODO | | | |
| ST7 the TUI takes focus; b | TODO | | | |
| ST8 click picking, every gate | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 208 | 212 | 0 failure(s) | 38 gates, 0 failed |
| ST1 | 216 | 218 | 0 failure(s) | (not run; ST8 only) |
| ST2 | 224 | 218 | 0 failure(s) | (not run; ST8 only) |
| ST3 | 231 | 218 | 0 failure(s) | (not run; ST8 only) |
| ST4 | 231 | 224 | 0 failure(s) | (not run; ST8 only) |
| ST5 | 231 | 231 | 0 failure(s) | (not run; ST8 only) |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

- 2026-10-01T17:54Z ST1 DONE (5816182). Reviewer PASS round 1, no findings. Real tree: `distance 128-132: 2.299 Å (C-H)`, `angle 128-132-133: 46.3° (C-H-O)`, `dihedral 125-128-132-133: 115.0° (H-C-H-O)`, all exit 0.
- 2026-10-01T18:17Z ST2 DONE (db8340d, 69a0d63). Reviewer r1 REVISE on one finding: `geometry_address` returned `[]` for a `PathPoint` while following; fixed in 69a0d63 by testing `PathPoint` before the `following` guard. Reviewer PASS round 2, no findings. 224/218 pass, 0 failures.
- 2026-10-01T18:21Z ST3 DONE (04027ae). Reviewer PASS round 1, no findings. 231/218 pass, 0 failures.
- 2026-10-01T18:26Z ST4 DONE (56fe103). Reviewer PASS round 1, no findings. 231/224 pass, 0 failures.
- 2026-10-01T18:36Z ST5 DONE (c908f4e). Reviewer PASS round 1, no findings. 231/231 pass, 0 failures.

## Backlog

- An untracked `tools/orcamon/LANDSCAPE.md` appeared in the working tree during ST1's review (mtime 17:53:24), not from any plan task and not in ST1's commit. The user committed it as 9f0af83 before ST2 started.
- No test covers `geometry_address` for a `PathPoint` with `following=True`; the spec fixes the ST2 test at 8 checks, so add a check in ST5 or ST7 (reviewer, ST2 r2).
- The ST3 selection test does not check that a picked atom is labelled exactly once (that `_draw_labels` skips it); pixel-level double-drawing is hard to see (reviewer, ST3 r1).
- ST4's commit carries a Co-Authored-By line naming DeepSeek rather than the harness attribution; cosmetic (reviewer, ST4 r1). Same trailer on ST5's commit.
- ST5's commit message c908f4e lost its backticked identifiers to shell substitution (amending was denied); it reads awkwardly but contains `Task: ST5` (reviewer, ST5 r1).
