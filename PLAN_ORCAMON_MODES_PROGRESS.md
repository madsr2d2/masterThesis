# orcamon — every normal mode, drawn and described — Progress

Plan: `PLAN_ORCAMON_MODES.md`. Branch `master` from `master` (`c0db5ef`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| NM1 Keep every mode, compactly | DONE | 1184c38 | 1 | reviewer PASS; test_monitor 200, test_orcamon 197, test_curve_metrics 0; deviation accepted: `_feed` strips the fixture's closing blank line, tests feed it explicitly |
| NM2 Find, list and rank modes | DONE | 6a5bf2a | 1 | reviewer PASS; test_monitor 208, test_orcamon 197, test_curve_metrics 0 |
| NM3 snapshot/freqs --mode for any mode | DONE | 9a110ab | 1 | reviewer PASS; test_orcamon 205, test_monitor 208, test_curve_metrics 0; deviations accepted: `cmd_snapshot` docstring wording, README `--help` block refreshed |
| NM4 I in the TUI, every gate | DONE | aba827f | 1 | reviewer PASS; test_orcamon 212, test_monitor 208, test_curve_metrics 0; run_gates 38 gates, 0 failed (reviewer saw one non-reproducible test_monitor failure under load, passed on rerun) |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 193 | 197 | 0 failure(s) | 38 gates, 0 failed |
| after NM1 | 200 | 197 | 0 failure(s) | — |
| after NM2 | 208 | 197 | 0 failure(s) | — |
| after NM3 | 208 | 205 | 0 failure(s) | — |
| after NM4 | 208 | 212 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

- NM1: plan step 7 said feed `_modes_block(...)` with `_feed`, but `_feed` does `text.strip("\n")`, so the fixture's closing blank line is stripped and `_finish_modes` never runs (at the change's base `state.normal_modes` stayed `None`) → tests feed the blank line explicitly with `feed_line("")`; every plan hand-derived value still holds, no assertion or expected value changed. Reviewer accepted.

## Log

- 2026-10-01T16:16Z NM1 DONE (1184c38), reviewer PASS round 1.
- 2026-10-01T16:20Z NM2 DONE (6a5bf2a), reviewer PASS round 1.
- 2026-10-01T16:31Z NM3 DONE (9a110ab), reviewer PASS round 1; real-tree mode 6/7 rankings match the plan.
- 2026-10-01T17:10Z NM4 DONE (aba827f), reviewer PASS round 1; run_gates 38 gates, 0 failed.

## Backlog

- Plan step 7 should note that `_feed` strips the trailing blank line; NM2+ tasks that feed `_modes_block` through `_feed` will hit the same problem. (reviewer, NM1)
- `cmd_freqs` raises the "no geometry holds mode" `UsageError` in two places (once when `mode_geometry` returns None, once when `participation` returns None); a shared message helper would be tidier. (reviewer, NM3)
- `test_monitor.py` failed once inside a loaded `run_gates.py` run (8 parallel jobs) during the NM4 review and passed on rerun; it may have a timing sensitivity under load. The failing check was not identified. (reviewer, NM4)
