# Reference benzaldehyde UV spectrum (gas phase)

Source: [MPI-Mainz UV/VIS Spectral Atlas of Gaseous Molecules of Atmospheric
Interest](https://uv-vis-spectral-atlas-mainz.org/uvvis/cross_sections/Aromatic%20compounds/Aldehydes/C6H5CHO.spc)
(Keller-Rudek et al., *Earth Syst. Sci. Data* 2013), benzaldehyde (C6H5CHO)
gas-phase absorption cross sections, downloaded 2026-08-31.

- `iupac2021_benzaldehyde.txt` — IUPAC (2021) recommended cross sections,
  298 K, 220-395 nm, 5 nm steps. Two columns: wavelength (nm), cross section
  (cm^2/molecule).
- `xiang2009_benzaldehyde.txt` — Xiang, Zhu & Zhu (2009) individual cavity
  ring-down measurement, 295 K, 285-400 nm. Three columns: wavelength,
  cross section, uncertainty.
- `trost1997_285nm.txt` — single-point entry from Trost (1997), labelled
  "max" at 285 nm. **Treat as unreliable**: it converts to eps ~ 8900
  M^-1 cm^-1, the same order of magnitude as the *strong* 220-240 nm
  pi->pi* band, not the weak band this wavelength is supposed to represent
  — almost certainly a mislabelled/misdigitized database entry. Kept for
  the record but excluded from the comparison.

## Conversion to molar decadic extinction coefficient

`eps [M^-1 cm^-1] = sigma [cm^2/molecule] * NA / (1000 * ln10) = sigma * 2.6154e20`

| wavelength (nm) | IUPAC(2021) eps | Xiang(2009) eps |
|---|---|---|
| 275 | 594 | - |
| 280 | 450 | - |
| 285 | **479** | **555** |
| 290 | 31 | 24 |
| 295 | 14 | 12 |

Both independent sources agree: benzaldehyde's gas-phase eps at 285 nm is
**~480-555 M^-1 cm^-1**, sitting right at a steep edge — it drops by >10x
between 285 and 290 nm in both datasets (a real spectral feature, not a
compilation artifact, since it appears in both the recommended compilation
and the single-study Xiang data).

## Why this matters for the thesis calculation

This is **gas phase**, not the aqueous-buffer value the experiments actually
read (dataset reports 1230 M^-1 cm^-1 at 285 nm for BnOH's product, presumably
aqueous). It is the right *first* comparison point anyway: it isolates
method/functional error from solvent-model error. Our first-pass vertical
TD-DFT (wB97X-D3/def2-TZVP, CPCM(water), single-point stick spectrum
Gaussian-broadened at an assumed 3000 cm^-1 FWHM) gave eps(285 nm) of only
~1-4 M^-1 cm^-1 for benzaldehyde — roughly **100-500x too weak** even before
comparing to the aqueous number. Since this reference shows the *real*
gas-phase value is already ~500 M^-1 cm^-1, the deficit is not primarily a
solvent-model problem — it's that a bare vertical Franck-Condon transition
dipole moment radically underestimates this band's true intensity, which is
mostly vibronic (Herzberg-Teller) in origin. This is the direct evidence
motivating a vibronically-resolved (ORCA `%esd ABS`) spectrum simulation
instead of naive Gaussian-broadened stick spectra.
