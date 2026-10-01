# orcamon — depth cues for multilayer structures — Progress

Plan: `PLAN_ORCAMON_DEPTH_CUES.md`. Branch `master` from `master` (`3f799f4`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| DC1 Two layers, composited | DONE | c5717e8 | 1 | parity 12/12 identical; test_monitor 119, test_orcamon 154, curve_metrics 0 failure(s) |
| DC2 Per-layer depth fog | DONE | f041a88 | 1 | fades 0.036/0.564/0.275; test_monitor 122, test_orcamon 154, curve_metrics 0 failure(s) |
| DC3 Depth-gap halos | DONE | 34b280a | 1 | halo widths 1/2/5; test_monitor 128, test_orcamon 154, curve_metrics 0 failure(s) |
| DC4 Boundary host balls | DONE | c1cb836 | 1 | column heights 61/13/13; test_monitor 131, test_orcamon 154, curve_metrics 0 failure(s) |
| DC5 See-through toggle | DONE | d9d7ec1 | 1 | blend (171,34,34); test_monitor 133, test_orcamon 157, curve_metrics 0 failure(s) |
| DC6 Docs, gates, pictures | DONE | 50826d7 | 1 | gates 38/0, suites 133/157, 5 PNGs in /tmp/orcamon_depth_cues/pictures |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 119 | 154 | 0 failure(s) | 38 gates, 0 failed |
| DC1 | 119 | 154 | 0 failure(s) | 38 gates, 0 failed |
| DC2 | 122 | 154 | 0 failure(s) | 38 gates, 0 failed |
| DC3 | 128 | 154 | 0 failure(s) | 38 gates, 0 failed |
| DC4 | 131 | 154 | 0 failure(s) | 38 gates, 0 failed |
| DC5 | 133 | 157 | 0 failure(s) | 38 gates, 0 failed |
| DC6 | 133 | 157 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

## Log

- 2026-10-01T10:04Z DC1 DONE (c5717e8), reviewer PASS round 1; parity probe 12/12 identical, suites 119/154/0.
- 2026-10-01T10:05Z DC2 DONE (f041a88), reviewer PASS round 1; fades 0.036/0.564/0.275, suites 122/154/0.
- 2026-10-01T10:08Z DC3 DONE (34b280a), reviewer PASS round 1; halo widths 1/2/5 rows, suites 128/154/0.
- 2026-10-01T10:11Z DC4 DONE (c1cb836), reviewer PASS round 1; boundary ball 61/13/13 rows, suites 131/154/0.
- 2026-10-01T10:17Z DC5 DONE (d9d7ec1), reviewer PASS round 1; see-through blend (171,34,34), README byte-identical, suites 133/157/0.
- 2026-10-01T10:38Z DC6 DONE (50826d7), reviewer PASS round 1; run_gates 38/0 in 530s (implementer 545s), suites 133/157/0, five 900x750 PNGs rendered.

## Backlog

- DC6 reviewer: the README's new pixel-pane sentence is a long comma chain of four clauses; it carries all four points but would read better as separate sentences.
- DC6 reviewer: the two render timings (78.5 ms at 900x750, 27.9 ms at 329x315) were not re-measured; they appear only in the commit message and the plan sets no threshold.
