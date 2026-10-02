# orcamon — what the first agent session in the Steer loop got wrong — Progress

Plan: `PLAN_ORCAMON_STEER_FIXES.md`. Branch `master` from `master` (`727f7f4851248dbc3e169b82ea61f9c6aeda2b92`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| FS1 show counts the whole system | DONE | 16fb732 | 1 | reviewer PASS: test_monitor 243, test_orcamon 248, test_curve_metrics 0 failure(s). Deviation accepted: plan's real-tree command nests `n_atoms` under `job` (plan said `d['n_atoms']`); values 137/140/133/137 unchanged. |
| FS2 view gives atoms, elements, measure | DONE | af0abb4 | 1 | reviewer PASS: test_orcamon 259, test_monitor 243, test_curve_metrics 0 failure(s). No deviations. |
| FS3 b stops the agent's mode | TODO | | | |
| FS4 events import, path address test, every gate | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 237 | 248 | 0 failure(s) | 38 gates, 0 failed |
| FS1 | 243 | 248 | 0 failure(s) | not run until FS4 |
| FS2 | 243 | 259 | 0 failure(s) | not run until FS4 |

## Gates

## Deviations (plan said → evidence → what was done)

- FS1: plan's real-tree command reads `d['n_atoms']` → `orcamon show --json` emits `{"schema":1,"job":{…}}`, so that literal command raises `KeyError: 'n_atoms'` → implementer read `d['job']['n_atoms']`; acceptance values (137/140/133/137; QM 21/20/17/17) unchanged. Reviewer accepted.

## Log

- 2026-10-02T02:39:03Z FS1 DONE, commit 16fb732, 1 round, reviewer PASS.
- 2026-10-02T02:44:08Z FS2 DONE, commit af0abb4, 1 round, reviewer PASS.

## Backlog

- FS1: the plan's real-tree command should read `d['job']['n_atoms']` (the `show --json` document nests the report under `job`). Not fixed here; plan is immutable.
