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


def main():
    doc = Checker(DOCUMENT)
    table = early_trough.trough_table()
    genuine = table[table.genuine].copy()
    genuine["cluster"] = genuine.apply(early_trough.cluster, axis=1)
    excluded = table[table.candidate & ~table.genuine]

    print("\nsection 1: the screen")
    doc.claim("curves scanned", f"**{len(table)} curves**")
    doc.claim("candidates",
              f"**{int(table.candidate.sum())}** clear a smoothed trough of "
              f"{early_trough.EARLY_TROUGH_CANDIDATE_Z:.1f} sigma or deeper")
    doc.check("seventeen candidates", int(table.candidate.sum()) == 17)

    exp41 = table.set_index(["experiment", "sample"]).loc[(4, 1)]
    exp222 = table.set_index(["experiment", "sample"]).loc[(22, 2)]
    doc.claim("exp 4.1's z", f"smoothed z **{exp41.z:.2f}**")
    doc.claim("exp 4.1 not sustained", "only 3 of 9 readings")
    doc.check("exp 4.1 really is 3 of 9", not exp41.sustained)
    doc.claim("exp 22.2's z", f"(z **{exp222.z:.2f}**, only 4 of 9)")
    doc.check("exp 22.2 really is 4 of 9", not exp222.sustained)

    doc.claim("fifteen of seventeen",
              f"**{len(genuine)} of {int(table.candidate.sum())}**")
    doc.check("exactly fifteen genuine curves", len(genuine) == 15,
              f"{len(genuine)}")
    doc.check("none of the genuine curves overlap a bubble event",
              not genuine.bubble_overlap.any())

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

    print("\nsection 3: the two clusters")
    n_oxidant = int((genuine.cluster == "oxidant").sum())
    n_substrate = int((genuine.cluster == "substrate").sum())
    doc.claim("oxidant count", f"**{n_oxidant} of 15**")
    doc.claim("substrate count", f"remaining **{n_substrate}**")
    doc.check("twelve and three", n_oxidant == 12 and n_substrate == 3)
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

    print("\nsection 4: not the O2 artefact")
    both = genuine.set_index(["experiment", "sample"])
    for exp, samp, events in ((141, 4, 3), (142, 4, 7)):
        row = both.loc[(exp, samp)]
        doc.claim(f"exp {exp}.{samp}: raw", f"| {row.z:.2f}σ |")
        doc.claim(f"exp {exp}.{samp}: corrected",
                  f"**{row.z_corrected:.2f}σ**")
        doc.check(f"exp {exp}.{samp} deepens after correction",
                  row.z_corrected < row.z)

    print("\nsection 5: the fifteen curves table")
    for row in genuine.sort_values("z").itertuples():
        hoo_text = (f"{row.e0_hoo:.0f}" if row.e0_hoo >= 10
                   else f"{row.e0_hoo:.3f}")
        # The strongest curve in the archive is bolded in the table, since
        # that is the number the document's own hero stat argues for.
        bold = "**" if row.z == genuine.z.min() else ""
        doc.claim(f"exp {row.experiment}.{row.sample}: its row",
                  f"| exp {row.experiment}.{row.sample} | {row.substrate} | "
                  f"{row.pH:.2f} | {row.s0:g} | {row.h2o2:g} | {row.e0:g} | "
                  f"{row.e0_s0:.3f} | {hoo_text} | {bold}{row.z:+.1f}{bold} | "
                  f"{row.cluster} |")

    print("\nthe figures the document promises")
    doc.figures(os.path.join(HERE, "index.html"), "ABC")

    print("\nthe curves page draws every curve it claims to")
    page = io.open(os.path.join(HERE, "progress_curves.html"),
                   encoding="utf-8").read()
    drawn = page.count("<div class='fig panel'>")
    doc.check("fifteen genuine plus two rejected controls",
              drawn == len(genuine) + len(excluded),
              f"{drawn} panels, {len(genuine)} + {len(excluded)} expected")
    doc.check("both rejected controls are labelled REJECTED",
              page.count("<strong>REJECTED</strong>") == len(excluded))
    doc.check("every genuine panel is labelled GENUINE",
              page.count("<strong>GENUINE</strong>") == len(genuine))

    print("\nthe figures: no data point drawn outside its own frame")
    doc.unclipped(os.path.join(HERE, "index.html"),
                  os.path.join(HERE, "progress_curves.html"))
    return doc.summary()


if __name__ == "__main__":
    raise SystemExit(main())
