"""Change a geometry: turn or slide one side of a molecule about a bond,
never in place.

Every function here returns a NEW list of atoms and leaves the one it was
given alone -- the TUI may still be drawing that geometry, and a file only
ever gets the new one.

Which side moves is D2's: the smaller side, a tie moving the second atom's
side, or `move=ATOM` naming an atom of the edited bond whose side moves
instead. A fragment joined to neither side of the edited bond never moves.

Bonds in this module follow `bonded_pairs`, not the renderer's `bonds()` (D3):
two heavy atoms under the renderer's `BOND_CUTOFF`, an atom and a hydrogen
under `X_H_BOND_CUTOFF`, and two hydrogens under `H_H_BOND_CUTOFF`. Hydrogen
bonds sit at 1.6-1.7 Å, so the renderer's cutoff would count them and an edit
that merely brought a hydrogen near an acceptor would be refused as making a
new bond. `bonds()` is unchanged on purpose: it draws what a viewer expects to
see, and this rule is only for edits and connectivity.

An angle here is the change of EVERY dihedral X-I-J-Y across the edited bond
(D4) -- X on I's side, Y on J's side -- turning the moving side about the axis
from the fixed atom to the moving one by the right-hand rule, so the same
degrees arrive whichever side moves. `set_dihedral` turns the bond until one
named dihedral reaches its target, wrapped into [-180, 180), and
`set_distance` slides one side along the line between the two atoms.
"""
from __future__ import annotations

import math

from .geometry import bonds

# D3: a bond to a hydrogen is much shorter than a heavy-atom bond, and 1.4 Å
# keeps the 1.6-1.7 Å hydrogen bonds out -- under the renderer's 1.7 Å a
# rotation would refuse itself the moment it made one.
X_H_BOND_CUTOFF = 1.4
# D3: the renderer skips H-H outright, but two hydrogens can end up in one
# 1 Å cell in a rotation and must still be a real, close pair.
H_H_BOND_CUTOFF = 1.0


def bonded_pairs(atoms: list) -> list[tuple[int, int]]:
    """The edit and connectivity bond rule (D3): the renderer's pairs with
    the hydrogen bonds dropped, plus the hydrogen pairs it does not look for.

    Two heavy atoms bond under `geometry.BOND_CUTOFF`; an atom and a hydrogen
    under `X_H_BOND_CUTOFF`; two hydrogens under `H_H_BOND_CUTOFF`. Pairs are
    `(i, j)` with `i < j` and the list is sorted. The renderer's `bonds()` is
    unchanged on purpose -- an edit must not be refused because a rotation
    brought a hydrogen near an acceptor.
    """
    kept = []
    for i, j in bonds(atoms):
        if atoms[i][0] == "H" or atoms[j][0] == "H":
            if math.dist(atoms[i][1:], atoms[j][1:]) >= X_H_BOND_CUTOFF:
                continue
        kept.append((i, j))
    # H-H pairs are found by binning the hydrogens into 1 Å cells and looking
    # in the 27 neighbouring ones, the same grid idea as `bonds()`.
    grid: dict[tuple[int, int, int], list[int]] = {}
    for i, (element, x, y, z) in enumerate(atoms):
        if element == "H":
            grid.setdefault((math.floor(x), math.floor(y), math.floor(z)), []).append(i)
    cutoff2 = H_H_BOND_CUTOFF * H_H_BOND_CUTOFF
    for (cx, cy, cz), members in grid.items():
        near = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    near.extend(grid.get((cx + dx, cy + dy, cz + dz), ()))
        for i in members:
            for j in near:
                if j <= i:
                    continue
                if math.dist(atoms[i][1:], atoms[j][1:]) ** 2 < cutoff2:
                    kept.append((i, j))
    return sorted(kept)
