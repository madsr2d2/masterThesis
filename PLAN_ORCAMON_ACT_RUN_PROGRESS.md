# orcamon — the agent starts and stops ORCA when the person says so — Progress

Plan: `PLAN_ORCAMON_ACT_RUN.md`. Branch `master` from `master` (`986aae3`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| AR1 the runner | DONE | 22018f3 | 1 | Reviewer PASS r1; 290/266/0; no deviations accepted |
| AR2 orcamon run | DONE | 5aab04b | 1 | Reviewer PASS r1; 299/266/0; accepted deviations: show check reads `doc["job"]["status"]`; --json rides warnings in the doc |
| AR3 orcamon stop, the skill's rule, every gate | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 266 | 285 | 0 failure(s) | 38 gates, 0 failed |
| after AR1 | 266 | 290 | 0 failure(s) | not run (AR3) |
| after AR2 | 266 | 299 | 0 failure(s) | not run (AR3) |

## Gates

## Deviations (plan said → evidence → what was done)

- AR2: plan wrote `show water --json` reads `doc["status"]`; the real document shape is `doc["job"]["status"]`. Reviewer accepted the test reading `doc["job"]["status"]`.
- AR2: in `--json` mode the `check_input` warnings are not printed as `warning:` lines (they ride in the `warnings` list) so stdout is one JSON document. Reviewer accepted as fitting step 5's "JSON instead of text".
- AR2: an extra refusal `job.input is None` → `no readable input` guards an AttributeError before `check_input`. Reviewer accepted (no existing behaviour changes).

## Log

- 2026-10-02T07:31Z — AR1 DONE (22018f3): reviewer PASS round 1; test_monitor 266 pass, test_orcamon 290 pass, test_curve_metrics 0 failure(s); run_gates deferred to AR3. No fake ORCA or `sleep 60` left behind.
- 2026-10-02T07:49Z — AR2 DONE (5aab04b): reviewer PASS round 1; test_monitor 266 pass, test_orcamon 299 pass, test_curve_metrics 0 failure(s); both real-tree dry runs matched (ts5 exit 0, optts_freq_tight exit 2); skill current.

## Backlog

- `runner.alive` treats a zombie as gone (`Z` counts as gone) while `steer.alive` counts a zombie as live. Intended by D7, but the two names could confuse a later reader.
- The README's "Changing a geometry" section still says "run it only when they say so" (the plan only asked for the SKILL.md wording change).
