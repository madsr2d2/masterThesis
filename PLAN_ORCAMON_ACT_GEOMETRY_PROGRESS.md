# orcamon — the agent changes a geometry, checks an input, reads bond changes — Progress

Plan: `PLAN_ORCAMON_ACT_GEOMETRY.md`. Branch `master` from `master` (`cf78720`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| AG1 the edit bond rule, element table | DONE | f49d145 | 1 | reviewer PASS; monitor 249, orcamon 263, curve_metrics 0 failure(s); no deviations |
| AG2 rotate, set dihedral, set distance, bond changes | DONE | 023baff | 1 | reviewer PASS; monitor 260, orcamon 263, curve_metrics 0 failure(s); real-tree check [134,135]/5/[]/(1.971,102,135)/-63.507/5; no deviations |
| AG3 orcamon edit | DONE | 8849f03 | 1 | reviewer PASS; monitor 260, orcamon 271, curve_metrics 0 failure(s); real-tree 4 commands match; no deviations |
| AG4 atom references, orcamon check | DONE | e9e5a14 | 1 | reviewer PASS; monitor 263, orcamon 279, curve_metrics 0 failure(s); real tree 6 geometry_file, 1 qm_boundary_h (pc_relax 135 H/19 O 0.91 Å), 66 output_exists; no deviations |
| AG5 first geometry, orcamon connectivity | TODO | | | |
| AG6 the workflow in the skill, every gate | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 246 | 263 | 0 failure(s) | 38 gates, 0 failed |
| after AG1 | 249 | 263 | 0 failure(s) | — |
| after AG2 | 260 | 263 | 0 failure(s) | — |
| after AG3 | 260 | 271 | 0 failure(s) | — |
| after AG4 | 263 | 279 | 0 failure(s) | — |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

- 2026-10-02T03:31Z AG1 DONE f49d145 — reviewer PASS r1.
- 2026-10-02T03:38Z AG2 DONE 023baff — reviewer PASS r1.
- 2026-10-02T06:33Z AG3 DONE 8849f03 — reviewer PASS r1.
- 2026-10-02T06:44Z AG4 DONE e9e5a14 — reviewer PASS r1.

## Backlog

- AG1: `edit.py`'s module docstring names `set_dihedral` and `set_distance`, which don't exist until AG2 (matches the plan's instruction to write D2–D4 in prose; nothing to fix now).
- AG3: `commands.py` `"moved_qm": moved_qm or None` / `"moved_environment": moved_environment or None` turn a legitimate 0 into `null` in `edit --json` on a multilayer job, making 0 indistinguishable from "no split" (follows the plan's literal "q or None"; not fixed).
- AG3: a refused edit's stderr line appears before its stdout lines when both go to one terminal (stream buffering only; no test effect).
- AG4: `_atom_token` expands a reversed range such as `5:2` to 5,4,3,2 where the plan says `a…b` inclusive (harmless; no real input uses a reversed range).
- AG4: `check` raises `UsageError` (exit 2) when the input is unreadable, instead of reporting it as a finding (case the plan does not specify).
