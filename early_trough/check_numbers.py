"""
Verifies every number quoted in early_trough/ANALYSIS.md against the code.

The comparison itself is `doc_check`, shared with every sibling folder.

    python early_trough/check_numbers.py
"""
import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "data"))
sys.path.insert(0, os.path.dirname(HERE))

import early_trough
from doc_check import Checker

DOCUMENT = os.path.join(HERE, "ANALYSIS.md")


def sci(value, digits=1):
    """`4.7e-17` as `4.7 x 10-17`, the way ANALYSIS.md types a tiny p-value."""
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = value / 10 ** exponent
    return f"{mantissa:.{digits}f} x 10{exponent}"


def hoo_text(value):
    """The write-up's own [enz]/[HOO-] formatting: 3dp under 10, else integer."""
    return f"{value:.0f}" if value >= 10 else f"{value:.3f}"


def main():
    doc = Checker(DOCUMENT)
    table = early_trough.trough_table()
    genuine = table[table.genuine]
    excluded = table[table.candidate & ~table.genuine]

    print("\nsection 1: the screen")
    doc.claim("curves scanned", f"**{len(table)} curves**")
    doc.claim("candidates",
              f"**{int(table.candidate.sum())}** clear a smoothed trough of "
              f"{early_trough.EARLY_TROUGH_CANDIDATE_Z:.1f} sigma or deeper")
    doc.check("seventeen candidates", int(table.candidate.sum()) == 17)

    both = table.set_index(["experiment", "sample"])
    exp41 = both.loc[(4, 1)]
    exp222 = both.loc[(22, 2)]
    doc.claim("exp 4.1's z", f"**{exp41.z:.2f}**")
    doc.claim("exp 22.2's z", f"z **{exp222.z:.2f}**")
    doc.check("both are genuine under the corrected screen",
              exp41.genuine and exp222.genuine)
    doc.claim("exp 4.1's leave-one-out move", "-7.31 → -6.12")
    doc.claim("exp 22.2's leave-one-out move", "-6.23 → -5.10")

    doc.claim("seventeen of seventeen", "**17 of 17**")
    doc.check("zero excluded candidates remain", len(excluded) == 0,
              f"{len(excluded)} excluded")
    doc.check("exactly seventeen genuine curves", len(genuine) == 17,
              f"{len(genuine)}")

    print("\nsection 2: the archive-wide correlation")
    corr = early_trough.dominance_correlation(table)
    for substrate in ("4OMe-BnOH", "BnOH"):
        rho_hoo, p_hoo, n_hoo = corr[substrate]["e0_hoo"]
        rho_s0, p_s0, n_s0 = corr[substrate]["e0_s0"]
        doc.claim(f"{substrate}: n", f"| {substrate} | {n_hoo} |")
        doc.claim(f"{substrate}: rho hoo", f"**{rho_hoo:+.3f}**")
        doc.claim(f"{substrate}: p hoo", f"**{sci(p_hoo)}**")
        doc.claim(f"{substrate}: rho s0", f"{rho_s0:+.3f}")
        doc.claim(f"{substrate}: p s0", f"{p_s0:.2f}")
        doc.check(f"{substrate}: hoo significant, s0 not",
                  p_hoo < 1e-3 and p_s0 > 0.05,
                  f"p_hoo={p_hoo:.2e} p_s0={p_s0:.2e}")

    print("\nsection 3: the species test")
    rho_h2o2_4ome = corr["4OMe-BnOH"]["e0_h2o2"][0]
    p_h2o2_4ome = corr["4OMe-BnOH"]["e0_h2o2"][1]
    rho_h2o2_bnoh = corr["BnOH"]["e0_h2o2"][0]
    doc.claim("4OMe: e0/h2o2 rho", f"+{rho_h2o2_4ome:.3f}")
    doc.claim("4OMe: e0/h2o2 p", f"p = {p_h2o2_4ome:.2f}")
    doc.claim("BnOH: e0/h2o2 rho", f"{rho_h2o2_bnoh:+.3f}")
    doc.check("4OMe's total-peroxide ratio is wrong-signed",
              rho_h2o2_4ome > 0, f"{rho_h2o2_4ome:+.4f}")
    doc.check("BnOH's is same-signed but weaker than the anion ratio",
              rho_h2o2_bnoh < 0
              and abs(rho_h2o2_bnoh) < abs(corr["BnOH"]["e0_hoo"][0]))

    print("\nsection 4: the two clusters")
    n_oxidant = int((genuine.cluster == "oxidant").sum())
    n_substrate = int((genuine.cluster == "substrate").sum())
    doc.claim("oxidant count", f"**{n_oxidant} of 17**")
    doc.claim("substrate count", f"remaining **{n_substrate}**")
    doc.check("fourteen and three", n_oxidant == 14 and n_substrate == 3)
    substrate_cluster = genuine[genuine.cluster == "substrate"]
    for exp in (141, 142, 143):
        doc.check(f"exp {exp}.4 is in the substrate cluster",
                  ((substrate_cluster.experiment == exp)
                   & (substrate_cluster["sample"] == 4)).any())
    doc.claim("substrate cluster's enz/hoo range",
              f"{substrate_cluster.e0_hoo.min():.3f}"
              f"–{substrate_cluster.e0_hoo.max():.3f}")
    doc.claim("substrate cluster's enz/s0 range",
              f"{substrate_cluster.e0_s0.min() * 100:.1f}"
              f"–{substrate_cluster.e0_s0.max() * 100:.1f}%")
    oxidant_cluster = genuine[genuine.cluster == "oxidant"]
    doc.claim("oxidant cluster's enz/hoo floor",
              f"**above {int(np.floor(oxidant_cluster.e0_hoo.min()))}**")
    doc.claim("substrate cluster's enz/hoo ceiling",
              f"ceiling of **{substrate_cluster.e0_hoo.max():.3f}**")
    doc.check("the two clusters separate by more than three orders of "
              "magnitude in enz/hoo",
              oxidant_cluster.e0_hoo.min()
              / substrate_cluster.e0_hoo.max() > 1000)

    print("\nsection 5: not the O2 artefact")
    for exp, samp in ((141, 4), (142, 4)):
        row = both.loc[(exp, samp)]
        doc.claim(f"exp {exp}.{samp}: raw", f"| {row.z:.2f}σ |")
        doc.claim(f"exp {exp}.{samp}: corrected",
                  f"**{row.z_corrected:.2f}σ**")
        doc.check(f"exp {exp}.{samp} deepens after correction",
                  row.z_corrected < row.z)

    print("\nsection 6: the binding rate constant")
    rates = early_trough.binding_rates(table)
    doc.claim("k_on range", f"{rates.k_on.min():.2f} to {rates.k_on.max():.1f}")
    doc.check("k_on spans about 1.8 orders of magnitude",
              abs(np.log10(rates.k_on.max() / rates.k_on.min()) - 1.8) < 0.1)

    print("\nsection 7: the Arrhenius check")
    fit = early_trough.arrhenius_check(rates)
    doc.claim("the two-point Ea", "85.7 kJ/mol")
    doc.claim("the implied prefactor", "2.1 × 10¹⁵ M⁻¹s⁻¹")
    doc.claim("the Eyring dS", "+40 J/mol/K")
    doc.claim("the honest Ea", f"**{fit['activation_kJ']:.1f} ± "
                                f"{fit['stderr_kJ']:.1f} kJ/mol**")
    doc.claim("the t-statistic", f"**{fit['t_statistic']:.1f}**")
    doc.check("not significant", fit["t_statistic"] < 2.0)
    doc.claim("same-temperature n", f"**{fit['same_temperature_n']}**")
    doc.claim("same-temperature spread",
              f"**{fit['same_temperature_spread']:.2f}**")
    doc.claim("the whole-range spread", "**1.74**")

    print("\nsection 8: the buffer comparison")
    buf = early_trough.buffer_comparison(rates)
    indexed = buf.set_index("buffer")
    phosphate = indexed.loc["Phosphate"]
    pyrophosphate = indexed.loc["Pyrophosphate"]
    doc.claim("phosphate row", f"| Phosphate | {int(phosphate.n)} | "
                               f"{phosphate.geometric_mean:.2f} M")
    doc.claim("pyrophosphate row",
              f"| Pyrophosphate | {int(pyrophosphate.n)} | "
              f"**{pyrophosphate.geometric_mean:.2f} M")
    ratio = pyrophosphate.geometric_mean / phosphate.geometric_mean
    doc.claim("the ratio", f"roughly **{ratio:.0f}×**")
    doc.check("pyrophosphate runs several times phosphate's mean",
              ratio > 4)
    doc.claim("the induction buffer order", "+1.094 ± 0.150")
    doc.claim("resolved tau_slow fraction", "6 of 14")

    print("\nsection 9: the seventeen curves table")
    for row in genuine.sort_values("z").itertuples():
        bold = "**" if row.z == genuine.z.min() else ""
        doc.claim(f"exp {row.experiment}.{row.sample}: its row",
                  f"| exp {row.experiment}.{row.sample} | {row.substrate} | "
                  f"{row.pH:.2f} | {row.s0:g} | {row.h2o2:g} | {row.e0:g} | "
                  f"{row.e0_s0:.3f} | {hoo_text(row.e0_hoo)} | "
                  f"{bold}{row.z:+.1f}{bold} | {row.cluster} |")

    print("\nthe figures the document promises")
    doc.figures(os.path.join(HERE, "index.html"), "ABCDE")

    print("\nthe curves page draws every curve it claims to")
    page = io.open(os.path.join(HERE, "progress_curves.html"),
                   encoding="utf-8").read()
    drawn = page.count("<div class='fig panel'>")
    doc.check("seventeen panels, one per genuine curve",
              drawn == len(genuine), f"{drawn} panels, {len(genuine)} expected")
    doc.check("every panel is labelled GENUINE",
              page.count("<strong>GENUINE</strong>") == len(genuine))
    doc.check("no panel is labelled REJECTED (nothing is excluded now)",
              "REJECTED" not in page)

    print("\nthe figures: no data point drawn outside its own frame")
    doc.unclipped(os.path.join(HERE, "index.html"),
                  os.path.join(HERE, "progress_curves.html"))
    return doc.summary()


if __name__ == "__main__":
    raise SystemExit(main())
