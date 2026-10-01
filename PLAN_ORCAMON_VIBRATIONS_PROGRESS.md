# orcamon — animate imaginary normal modes — Progress

Plan: `PLAN_ORCAMON_VIBRATIONS.md`. Branch `master` from `master` (`6729d55`).

| Task | Status | Commits | Rounds | Note |
|---|---|---|---|---|
| MO1 Parse NORMAL MODES | DONE | d933363 | 1 | real jobs parse (N=420/51); test_monitor 141, test_orcamon 157, curve_metrics 0 failure(s) |
| MO2 Map and displace a mode | BLOCKED | none (uncommitted) | 1 | Probe cannot map: `state.atoms` is 22 for verify_exact_hessian (mode N=140), and the plan's parity statistic is ill-conditioned (spread 0.972) although D3 is confirmed by d/(2*printed)=1.000000. NEEDS_USER. |
| MO3 `snapshot --mode` | TODO | | | |
| MO4 TUI `i` animation | TODO | | | |
| MO5 Docs, gates, pictures | TODO | | | |

Statuses: TODO | DONE | BLOCKED

## Suite status

| After | test_monitor pass | test_orcamon pass | test_curve_metrics | run_gates |
|---|---|---|---|---|
| baseline | 136 | 157 | 0 failure(s) | 38 gates, 0 failed |
| MO1 | 141 | 157 | 0 failure(s) | 38 gates, 0 failed |

## Gates

## Deviations (plan said → evidence → what was done)

- Plan MO2 Facts said `verify_exact_hessian` maps directly because `len(job.xyz)=140`; measured, `state.atoms` is 22 (the `.out` prints only the QM1+link block; the 140-atom structure is in `job.xyz` only), so the direct map finds nothing to map onto. Open: which geometry a mode is drawn on.
- Plan MO2 probe's parity statistic `d/offs[i]` over `|offs|>1e-9` is ill-conditioned at six-decimal printing (spread 0.972). The convention D3 holds: `d/(2*printed)` is 1.000000 within printing precision for C, H and O alike. Statistic to be restated, not the convention.

- 2026-10-01T11:20Z MO1 DONE (d933363), reviewer PASS round 1; NORMAL MODES parsed from real jobs, suites 141/157/0.
- 2026-10-01T11:35Z MO2 NEEDS_USER, no commits; vibrations.py + 9 checks are green and uncommitted in the tree. Two plan defects found by the ground-truth probe (see Deviations). Awaiting the user's decision on which geometry a mode is drawn on.

## Backlog

- MO2: `verify_exact_hessian`'s displayed geometry is its 22-atom QM1+link block, so a 140-atom Hessian mode cannot map onto it; the full 140-atom structure is `job.xyz` only.
- MO2: for `optts_freq_tight`, `state.atoms` is one optimizer step behind `job.xyz` (worst gap 1.79e-4 A), while `job.xyz` matches `activeRegion.xyz` exactly.
