# Superseded by `../optts_freq/`

This was the C8 B97-3c TS done as **two jobs**: `optts/` (`! OptTS`, 46 cycles,
6:41) and `optts_numfreq/` (`! NumFreq` at the resulting geometry, 4:03). It
gave `transition state`, one imaginary mode at −379.28 cm⁻¹.

Superseded on 2026-09-25 by the single-job route (`! OptTS Freq`), which does
the optimisation and the frequency in one input and one `.out` — no geometry
file to carry between jobs, and ORCA refuses the Freq if the optimisation did
not converge.

**This two-job pair is also anomalous, and that is the interesting part.** Its
lowest real modes are 49.0 and 91.1 cm⁻¹, against 130–137 cm⁻¹ for the
single-job runs, and its G(298.15) differs from them by 30 kJ/mol — all of it in
the ZPE. It had converged to a different low-frequency well. The canonical C8 TS
is `../optts_freq_tight/` (`! OptTS Freq TightOpt`, 41 cycles, 9:15, imaginary
mode −369.0 cm⁻¹, G = −797.042444 Eh); `../optts_freq/` (default tolerance)
agrees with it to 0.09 kJ/mol on G.

This pair is kept only as the comparison. Do not quote a number from here.
