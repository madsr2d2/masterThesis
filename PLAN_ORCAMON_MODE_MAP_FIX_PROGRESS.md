# orcamon — draw a mode on a structure it fits — Progress

Plan: `PLAN_ORCAMON_MODE_MAP_FIX.md`. Branch `master` from `master` (`d3ac99e`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| MF1 Fit the QM map, name the cause | DONE | 9b3a302 | 1 | probe: IndexError before, 137 atoms / 0.2000 after; test_monitor 155, test_orcamon 164, run_gates 38/0 |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 153 | 163 | 0 failure(s) | 38 gates, 0 failed |
| MF1 | 155 | 164 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

## Deviations (plan said → evidence → what was done)

- `tui/app.py`: the plan's literal `self._job.state.frequencies` would raise `AttributeError` when no job is shown; the code guards `self._job is not None` first.
- `test_orcamon.py`: step 3 changes `opt_maxiter --mode 0` to the new "no frequency block" message, so the existing `"available: none"` check was retargeted to `opt_done` (a job with a frequency block but no imaginary mode). The assertion expression and expected value are unchanged.

## Log

- 2026-10-01T13:02Z MF1 DONE (9b3a302), reviewer PASS round 1; `optts_numfreq` maps via job.xyz (137 atoms, 0.2 Å peak); suites 155/164/0, run_gates 38/0.

## Backlog

- MF1 reviewer: the plan's step 5 / Never list / acceptance conflict over the `opt_maxiter` check; the plan should note the old "no frequencies" check moves to a job with a frequency block.
