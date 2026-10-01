# orcamon — reaction paths: IRC and NEB — Progress

Plan: `PLAN_ORCAMON_PATHS.md`. Branch `master` from `master` (`ec5144624628357b384788f65278fa8d282a2c70`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| PA1 Parse the IRC rows | DONE | f99dc53 | 1 | suite: monitor 168/0, orcamon 178/0, curve_metrics 0 failures; no deviations |
| PA2 Build the IRC path | DONE | 6af91fa | 1 | suite: monitor 179/0, orcamon 178/0, curve_metrics 0 failures; reviewer accepted the row-energy and missing-file-cache-key interpretations |
| PA3 Summary and irc_not_converged | DONE | 7d93dc9 | 1 | suite: monitor 184/0, orcamon 179/0, curve_metrics 0 failures; accepted deviation: test_orcamon imports `_irc_tree` from test_monitor (the duplicate guard reads private functions too, so the plan's "copy it" premise was wrong) |
| PA4 energies, geom/snapshot --point | DONE | cdd7399 | 1 | suite: monitor 184/0, orcamon 187/0, curve_metrics 0 failures; accepted deviations: private `_energies_path` extraction (output identical), README help blocks pasted at COLUMNS=90 to match the other blocks |
| PA5 Scrub the path in the TUI | BLOCKED | none (uncommitted changes discarded) | 1 | `end` does not re-place the focus: with step 5 as written `_following_latest` toggles but `_signature_of` is unchanged, so `show_job` early-returns and `_place_selection` is never called; the geometry pane stays on `IRC backward 0` where the plan's 5th check expects `IRC TS` (test_orcamon 191 pass / 1 failure). Baseline restored (187/0). |
| PB1 Parse the NEB and build its path | TODO | | | |
| PB2 NEB end to end, every gate | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 161 | 178 | 0 failure(s) | 38 gates, 0 failed |
| PA1 | 168 | 178 | 0 failure(s) | |
| PA2 | 179 | 178 | 0 failure(s) | |
| PA3 | 184 | 179 | 0 failure(s) | |
| PA4 | 184 | 187 | 0 failure(s) | |

## Gates

## Blocked

- **PA5** — implementing plan step 5 exactly as written, the plan's own hand-derived expected value does not hold: after `end` the geometry pane still shows `IRC backward 0`, not `IRC TS`. Cause: `action_select_latest` (path mode) sets `_following_latest = True` and calls `update_detail()`, but `_render_detail` → `show_job` early-returns because `_signature_of(job)` does not include `_following_latest`, so `_place_selection` is never reached and `selected_index` stays where `left,left` put it. Checks 1–4 pass; check 5 fails (test_orcamon 191/187+... see row). No commit was made; the uncommitted changes were discarded and the baseline re-verified (test_monitor 184/0, test_orcamon 187/0, curve_metrics 0 failures). Needs a user decision on how to amend PLAN_ORCAMON_PATHS.md PA5 (see the question in the orchestrator's report).

## Deviations (plan said → evidence → what was done)

## Log


- 2026-10-01T14:21Z PA1 DONE (f99dc53); round 1 reviewer PASS. Real-tree prototype matched the plan exactly (irc_bridge forwards/backwards 50/50, last backward dE -61.448643; r2scan 20/20, -17.416051).
- 2026-10-01T14:27Z PA2 DONE (6af91fa); round 1 reviewer PASS. Real irc_bridge builds 101 points, 137 atoms each, focus 50; summary matches § Facts. Reviewer accepted `energy = row.energy` on the row points (consistent with D4/§ Facts) and the missing-file cache key.
- 2026-10-01T14:33Z PA3 DONE (7d93dc9); round 1 reviewer PASS. Real `show irc_bridge` prints the 50/50 (MaxIter) summary and the `! IRC forward and backward hit MaxIter (50)` flag; `ls` row shows `irc B50/F50`.

- 2026-10-01T14:42Z PA4 DONE (cdd7399); round 1 reviewer PASS. Real `energies irc_bridge --json` prints `irc 101 -50 IRC backward 49 -257.1 IRC TS -1011.288471006717 50 -5.27`; `geom --point -50` prints `137`. No JSON key added (`schema` comes from `_emit_json` for every branch).
- 2026-10-01T14:52Z PA5 BLOCKED (no commit). Plan stop condition: a hand-derived expected value does not hold when step 5 is implemented as written (`end` leaves the pane on `IRC backward 0`, not `IRC TS`). The implementer's uncommitted changes were discarded and the baseline re-verified (184/187/0). Stopped for a user decision on amending PA5.

## Backlog

- `test_monitor.py` has an extra blank line before `_TERMINATED` (three blank lines after the new test) — cosmetic, PA1 reviewer.
- The PA1 implementer's note that the `IRC PATH SUMMARY` block near `.out` line 3146 is closed by a blank line before the prefilter sees it could not be reproduced by the reviewer. It did not affect the exact real-tree row counts.
- `core/paths.py::_frames` caches a missing file as `[]` under a `None` stamp; those entries count against the 8 slots. Harmless (PA2 reviewer).
- `data/test_curve_metrics.py::_defined_names` docstring says anything starting with `_` is skipped, but only private constants are; private functions/classes duplicated across modules are reported. Correcting the docstring is a separate change outside this plan (PA3 implementer/reviewer).
- PLAN_ORCAMON_PATHS.md PA3 acceptance names `test_report_keys_are_stable` under test_orcamon.py; it lives in test_monitor.py (plan slip, PA3 reviewer).
- Text `energies` of a 101-point path is cut at the default `--max-lines`, so `dE from the TS` and the label footer are not shown unless `--max-lines` is raised. Existing `Out` behaviour the plan asks for, but an agent may want the footer kept (PA4 reviewer).
- `_path_point` has no `-> int` return annotation; other helpers in the file are similarly loose (style, PA4 reviewer).

