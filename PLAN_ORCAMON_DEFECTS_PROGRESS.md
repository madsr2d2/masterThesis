# orcamon — fix the defects the end-to-end review found — Progress

Plan: `PLAN_ORCAMON_DEFECTS.md`. Branch `master` from `master` (`4a6f7bee83d077b25ab3228c6626743293e27e09`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| FX1 Liveness per job, orcamon is not ORCA | DONE | 3c422bc | 1 | reviewer: test_orcamon 170 pass / test_monitor 155 pass / curve_metrics 0 failure(s); accepted deviation: test uses the file's `_json(...)` helper (= `_orcamon([... "--json"])` + `json.loads`), same argv and assertions |
| FX2 Full labels in ls | TODO | | | |
| FX3 No matplotlib for the pixel panes | TODO | | | |
| FX4 Only ORCA inputs are jobs | TODO | | | |
| FX5 Convergence reason, every gate | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 155 | 164 | 0 failure(s) | 38 gates, 0 failed |
| after FX1 | 155 | 170 | 0 failure(s) | not run |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

- 2026-10-01T13:24Z FX1 DONE at 3c422bc (round 1, reviewer PASS). Liveness is now per-job via the `orca` driver's named input, and `orcamon`'s own console script no longer counts as ORCA. Suites: test_orcamon 170 pass, test_monitor 155 pass, data/test_curve_metrics 0 failure(s); run_gates.py not run (reserved for FX5).

## Backlog
