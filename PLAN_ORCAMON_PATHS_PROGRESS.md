# orcamon — reaction paths: IRC and NEB — Progress

Plan: `PLAN_ORCAMON_PATHS.md`. Branch `master` from `master` (`ec5144624628357b384788f65278fa8d282a2c70`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| PA1 Parse the IRC rows | DONE | f99dc53 | 1 | suite: monitor 168/0, orcamon 178/0, curve_metrics 0 failures; no deviations |
| PA2 Build the IRC path | DONE | 6af91fa | 1 | suite: monitor 179/0, orcamon 178/0, curve_metrics 0 failures; reviewer accepted the row-energy and missing-file-cache-key interpretations |
| PA3 Summary and irc_not_converged | TODO | | | |
| PA4 energies, geom/snapshot --point | TODO | | | |
| PA5 Scrub the path in the TUI | TODO | | | |
| PB1 Parse the NEB and build its path | TODO | | | |
| PB2 NEB end to end, every gate | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 161 | 178 | 0 failure(s) | 38 gates, 0 failed |
| PA1 | 168 | 178 | 0 failure(s) | |
| PA2 | 179 | 178 | 0 failure(s) | |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

- 2026-10-01T14:21Z PA1 DONE (f99dc53); round 1 reviewer PASS. Real-tree prototype matched the plan exactly (irc_bridge forwards/backwards 50/50, last backward dE -61.448643; r2scan 20/20, -17.416051).
- 2026-10-01T14:27Z PA2 DONE (6af91fa); round 1 reviewer PASS. Real irc_bridge builds 101 points, 137 atoms each, focus 50; summary matches § Facts. Reviewer accepted `energy = row.energy` on the row points (consistent with D4/§ Facts) and the missing-file cache key.

## Backlog

- `test_monitor.py` has an extra blank line before `_TERMINATED` (three blank lines after the new test) — cosmetic, PA1 reviewer.
- The PA1 implementer's note that the `IRC PATH SUMMARY` block near `.out` line 3146 is closed by a blank line before the prefilter sees it could not be reproduced by the reviewer. It did not affect the exact real-tree row counts.
- `core/paths.py::_frames` caches a missing file as `[]` under a `None` stamp; those entries count against the 8 slots. Harmless (PA2 reviewer).

