# orcamon — the person and the agent share one view — Progress

Plan: `PLAN_ORCAMON_STEER.md`. Branch `master` from `master` (`9b4abeae965da17063bd4a866c93e8058ed29fd4`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| ST1 measure and the command | DONE | 5816182 | 1 | reviewer PASS r1: 216/218 pass, 0 failures; real-tree distance/angle/dihedral values matched |
| ST2 the state files | TODO | | | |
| ST3 renderers draw the selection | TODO | | | |
| ST4 s picks atoms | TODO | | | |
| ST5 the TUI publishes; orcamon view | TODO | | | |
| ST6 orcamon focus | TODO | | | |
| ST7 the TUI takes focus; b | TODO | | | |
| ST8 click picking, every gate | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 208 | 212 | 0 failure(s) | 38 gates, 0 failed |
| ST1 | 216 | 218 | 0 failure(s) | (not run; ST8 only) |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

- 2026-10-01T17:54Z ST1 DONE (5816182). Reviewer PASS round 1, no findings. Real tree: `distance 128-132: 2.299 Å (C-H)`, `angle 128-132-133: 46.3° (C-H-O)`, `dihedral 125-128-132-133: 115.0° (H-C-H-O)`, all exit 0.

## Backlog

- An untracked `tools/orcamon/LANDSCAPE.md` appeared in the working tree during ST1's review (mtime 17:53:24), not from any plan task and not in ST1's commit. Reviewer flagged it; do not stage it.
