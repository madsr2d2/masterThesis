# orcamon — depth cues for multilayer structures — Progress

Plan: `PLAN_ORCAMON_DEPTH_CUES.md`. Branch `master` from `master` (`3f799f4`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| DC1 Two layers, composited | DONE | c5717e8 | 1 | parity 12/12 identical; test_monitor 119, test_orcamon 154, curve_metrics 0 failure(s) |
| DC2 Per-layer depth fog | DONE | f041a88 | 1 | fades 0.036/0.564/0.275; test_monitor 122, test_orcamon 154, curve_metrics 0 failure(s) |
| DC3 Depth-gap halos | DONE | 34b280a | 1 | halo widths 1/2/5; test_monitor 128, test_orcamon 154, curve_metrics 0 failure(s) |
| DC4 Boundary host balls | DONE | c1cb836 | 1 | column heights 61/13/13; test_monitor 131, test_orcamon 154, curve_metrics 0 failure(s) |
| DC5 See-through toggle | TODO | | | |
| DC6 Docs, gates, pictures | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 119 | 154 | 0 failure(s) | 38 gates, 0 failed |
| DC1 | 119 | 154 | 0 failure(s) | 38 gates, 0 failed |
| DC2 | 122 | 154 | 0 failure(s) | 38 gates, 0 failed |
| DC3 | 128 | 154 | 0 failure(s) | 38 gates, 0 failed |
| DC4 | 131 | 154 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

- 2026-10-01T10:04Z DC1 DONE (c5717e8), reviewer PASS round 1; parity probe 12/12 identical, suites 119/154/0.
- 2026-10-01T10:05Z DC2 DONE (f041a88), reviewer PASS round 1; fades 0.036/0.564/0.275, suites 122/154/0.
- 2026-10-01T10:08Z DC3 DONE (34b280a), reviewer PASS round 1; halo widths 1/2/5 rows, suites 128/154/0.
- 2026-10-01T10:11Z DC4 DONE (c1cb836), reviewer PASS round 1; boundary ball 61/13/13 rows, suites 131/154/0.

## Backlog
