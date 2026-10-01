# orcamon — animate imaginary normal modes — Progress

Plan: `PLAN_ORCAMON_VIBRATIONS.md`. Branch `master` from `master` (`6729d55`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| MO1 Parse NORMAL MODES | DONE | d933363 | 1 | real jobs parse (N=420/51); test_monitor 141, test_orcamon 157, curve_metrics 0 failure(s) |
| MO2 Map and displace a mode | DONE | f9c1e48 | 1 | plan revised (a22eb06, cf53402): draws the geometry the mode belongs to; parity 1e-6; test_monitor 153, test_orcamon 157, curve_metrics 0 failure(s) |
| MO3 `snapshot --mode` | DONE | 58b610c | 1 | --mode/--phase write a PNG and error cleanly; test_monitor 153, test_orcamon 160, curve_metrics 0 failure(s) |
| MO4 TUI `i` animation | DONE | bc6d436 | 1 | i cycles imaginary modes, title names it; test_monitor 153, test_orcamon 163, curve_metrics 0 failure(s) |
| MO5 Docs, gates, pictures | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 136 | 157 | 0 failure(s) | 38 gates, 0 failed |
| MO1 | 141 | 157 | 0 failure(s) | 38 gates, 0 failed |
| MO2 | 153 | 157 | 0 failure(s) | 38 gates, 0 failed |
| MO3 | 153 | 160 | 0 failure(s) | 38 gates, 0 failed |
| MO4 | 153 | 163 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

- Plan MO2 Facts said `verify_exact_hessian` maps directly because `len(job.xyz)=140`; measured, `state.atoms` is 22 (the `.out` prints only the QM1+link block; the 140-atom structure is in `job.xyz` only), so the direct map finds nothing to map onto. Open: which geometry a mode is drawn on.
- Plan MO2 probe's parity statistic `d/offs[i]` over `|offs|>1e-9` is ill-conditioned at six-decimal printing (spread 0.972). The convention D3 holds: `d/(2*printed)` is 1.000000 within printing precision for C, H and O alike. Statistic to be restated, not the convention.
- Resolved (user chose option a): plan revised in `a22eb06` (a mode is drawn on the geometry it belongs to, falling back to the file the input names) and `cf53402` (the convention is checked directly as `max |d - 2*printed| <= 2e-6`); MO2 then passed at `f9c1e48`.

## Log

- 2026-10-01T11:20Z MO1 DONE (d933363), reviewer PASS round 1; NORMAL MODES parsed from real jobs, suites 141/157/0.
- 2026-10-01T11:35Z MO2 NEEDS_USER, no commits; vibrations.py + 9 checks green and uncommitted. Two plan defects found by the ground-truth probe (see Deviations).
- 2026-10-01T11:39Z MO2 DONE (f9c1e48), reviewer PASS round 1; user chose the geometry fallback (plan a22eb06), parity restated (cf53402); suites 153/157/0.
- 2026-10-01T11:44Z MO3 DONE (58b610c), reviewer PASS round 1; snapshot --mode/--phase, README and skill regenerated; suites 153/160/0.
- 2026-10-01T11:50Z MO4 DONE (bc6d436), reviewer PASS round 1; `i` cycles imaginary modes and names them in the title; suites 153/163/0.

## Backlog

- MO2 reviewer: `phase_sine` has no docstring (trivial).
- MO3 reviewer: `--phase` is not range-checked despite the "0-360" help (harmless).
- MO3 reviewer: the two `UsageError` messages for "no geometry holds mode" are duplicated verbatim.
- MO4 reviewer: the commit message body has the typo `cm-l`; the code uses `cm-1`.
