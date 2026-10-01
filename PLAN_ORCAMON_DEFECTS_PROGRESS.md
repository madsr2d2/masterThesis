# orcamon — fix the defects the end-to-end review found — Progress

Plan: `PLAN_ORCAMON_DEFECTS.md`. Branch `master` from `master` (`4a6f7bee83d077b25ab3228c6626743293e27e09`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| FX1 Liveness per job, orcamon is not ORCA | DONE | 3c422bc | 1 | reviewer: test_orcamon 170 pass / test_monitor 155 pass / curve_metrics 0 failure(s); accepted deviation: test uses the file's `_json(...)` helper (= `_orcamon([... "--json"])` + `json.loads`), same argv and assertions |
| FX2 Full labels in ls | DONE | 2b01c07 | 1 | reviewer: test_orcamon 173 pass / test_monitor 155 pass / curve_metrics 0 failure(s); no deviations |
| FX3 No matplotlib for the pixel panes | DONE | a6f9e1f | 1 | reviewer: test_orcamon 175 pass / test_monitor 155 pass / curve_metrics 0 failure(s); no deviations |
| FX4 Only ORCA inputs are jobs | DONE | 0b703af | 1 | reviewer: test_orcamon 178 pass / test_monitor 155 pass / curve_metrics 0 failure(s); real tree 85 -> 80 jobs, all 5 dropped are CREST `$constrain` files; no deviations |
| FX5 Convergence reason, every gate | DONE | 936fa84 | 1 | reviewer: test_monitor 161 pass / test_orcamon 178 pass / curve_metrics 0 failure(s) / run_gates 38 gates 0 failed; no deviations |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 155 | 164 | 0 failure(s) | 38 gates, 0 failed |
| after FX1 | 155 | 170 | 0 failure(s) | not run |
| after FX2 | 155 | 173 | 0 failure(s) | not run |
| after FX3 | 155 | 175 | 0 failure(s) | not run |
| after FX4 | 155 | 178 | 0 failure(s) | not run |
| after FX5 | 161 | 178 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

- FX1: the plan's new test calls `_orcamon([... "--json"])` and parses stdout. A different name/path for the same thing is not a stop, so the test calls the file's existing `_json(...)` helper, which is exactly `_orcamon` plus a guarded `json.loads` — same argv, same assertions. Reviewer accepted.

## Log

- 2026-10-01T13:24Z FX1 DONE at 3c422bc (round 1, reviewer PASS). Liveness is now per-job via the `orca` driver's named input, and `orcamon`'s own console script no longer counts as ORCA. Suites: test_orcamon 170 pass, test_monitor 155 pass, data/test_curve_metrics 0 failure(s); run_gates.py not run (reserved for FX5).
- 2026-10-01T13:29Z FX2 DONE at 2b01c07 (round 1, reviewer PASS). `ls` now prints every label whole; the real tree `ls` emits no `…` (74 before) and `show` of the full label exits 0. Suites: test_orcamon 173 pass, test_monitor 155 pass, data/test_curve_metrics 0 failure(s); run_gates.py not run (reserved for FX5). Reviewer backlog: test_monitor timing flake, FX2 round-trip guard.
- 2026-10-01T13:32Z FX3 DONE at a6f9e1f (round 1, reviewer PASS). `images_available` now checks numpy and PIL only, so a clean `[tui,images]` install keeps the pixel panes. Suites: test_orcamon 175 pass, test_monitor 155 pass, data/test_curve_metrics 0 failure(s); run_gates.py not run (reserved for FX5). Reviewer backlog: stale matplotlib in generated `PKG-INFO`.
- 2026-10-01T13:36Z FX4 DONE at 0b703af (round 1, reviewer PASS). A `.inp` is a job only when it reads as an ORCA input; the real tree drops 85 -> 80 and all five dropped labels are CREST `$constrain` files. Suites: test_orcamon 178 pass, test_monitor 155 pass, data/test_curve_metrics 0 failure(s); run_gates.py not run (reserved for FX5). No new backlog.
- 2026-10-01T14:01Z FX5 DONE at 936fa84 (round 1, reviewer PASS). The text report names ORCA's relaxed convergence rule; the JSON contract is unchanged. Suites: test_monitor 161 pass, test_orcamon 178 pass, data/test_curve_metrics 0 failure(s), run_gates.py 38 gates in 547s on 8 jobs, 0 failed. No new backlog.

## Backlog

- `test_monitor.py:1143` wall-clock check ("140 atoms at 940x900 in under 250 ms") flaked once under parallel load, passed 4/4 sequential and in the reviewer's sequential run. Possible timing flake; not fixed.
- The new FX2 test's third check (`show` with the full label) already passed at base; it is a round-trip guard, and the first two checks are the ones that fail at base. Noted, not changed.
- `tools/orcamon/src/orcamon.egg-info/PKG-INFO:48` still lists matplotlib for the `images` extra (generated metadata). Check whether it is tracked and regenerate it if so.
