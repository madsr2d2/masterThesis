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
from dataclasses import dataclass

from .geometry import bonds, measure

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


class EditError(ValueError):
    """The edit is impossible or ill-posed: the atoms are not bonded, the
    bond is in a ring, the two atoms do not straddle it, `--move` names an
    atom that is not the edited bond's, or a distance was not positive. The
    message says why, and nothing is written."""


@dataclass
class EditResult:
    """What one edit did, with the geometry left alone: the new atom list,
    the atoms that moved, and the two sides of the edited bond.

    `fixed` is the bond atom that stayed and `mover` the one that moved;
    `fixed_side` is the number of atoms on the fixed atom's side, which is
    the count the report quotes to say how to move the other side instead
    (D1). `new_bonds` are the bonds the motion formed between the moving and
    the fixed side, and `closest` is the nearest approach between them
    whether or not the bond rule calls it a bond -- the pair a caller refuses
    an edit over, and the one a person wants to see. `before` and `after` are
    the value a setter aimed at (a dihedral in degrees, a distance in
    angstrom) and are None for a plain rotation, which aims at nothing.
    """

    atoms: list
    moved: list
    fixed: int
    mover: int
    fixed_side: int
    new_bonds: list
    closest: tuple | None
    before: float | None = None
    after: float | None = None


def _component(n: int, pairs: list, start: int, cut: tuple) -> list[int]:
    """The atoms reachable from `start` through `pairs` without using the
    pair `cut`, sorted.

    This is what makes one side of a bond: remove that one pair and the
    molecule falls into the pieces an edit can move independently."""
    adjacency: list[list[int]] = [[] for _ in range(n)]
    for a, b in pairs:
        if (a, b) == cut:
            continue
        adjacency[a].append(b)
        adjacency[b].append(a)
    seen = {start}
    stack = [start]
    while stack:
        for neighbour in adjacency[stack.pop()]:
            if neighbour not in seen:
                seen.add(neighbour)
                stack.append(neighbour)
    return sorted(seen)


def _sides(atoms: list, pairs: list, i: int, j: int, need_bond: bool) -> tuple[list, list]:
    """The `(i's side, j's side)` of the pair i-j, refusing what cannot move.

    With `need_bond` the pair must be a bond (a rotation or a dihedral about
    two atoms that are not bonded is undefined). Either way the pair's two
    atoms must fall apart when it is cut: if they stay connected it is a ring
    bond -- both sides include the other, so neither can move alone -- or,
    with `need_bond` false, two atoms of one connected fragment, where
    sliding one side would bend the fragment instead of moving it (D4)."""
    n = len(atoms)
    key = (min(i, j), max(i, j))
    d = math.dist(atoms[i][1:], atoms[j][1:])
    if need_bond and key not in pairs:
        raise EditError(f"atoms {i} and {j} are not bonded: {d:.2f} Å apart")
    side_j = _component(n, pairs, j, key)
    if i in side_j:
        if key in pairs:
            raise EditError(f"the bond {i}-{j} is in a ring: both sides stay connected, so neither can move alone")
        raise EditError(f"atoms {i} and {j} are {d:.2f} Å apart in one connected fragment: moving either would bend it")
    return _component(n, pairs, i, key), side_j


def _choose(side_i: list, side_j: list, i: int, j: int, move: int | None) -> tuple[int, int, list]:
    """`(fixed, mover, moved)`: D2's smaller side, a tie moving j's side, or
    `move=ATOM` naming an atom of the edited bond whose side moves instead.

    `move` must be i or j whatever the sizes say, so an atom picked outside
    the edited bond is told which two it may be rather than silently given
    the default."""
    if move is None:
        return (i, j, side_j) if len(side_j) <= len(side_i) else (j, i, side_i)
    if move == j:
        return i, j, side_j
    if move == i:
        return j, i, side_i
    raise EditError(f"--move takes {i} or {j}, the atoms of the edited bond")


def _rotate(atoms: list, moved: list, fixed: int, mover: int, degrees: float) -> list:
    """`atoms` with every atom in `moved` turned about the axis from `fixed`
    to `mover`, through `fixed`, by `degrees` under the right-hand rule (D4).

    Rodrigues' formula, `v' = v cos t + (u x v) sin t + u (u . v)(1 - cos t)`
    with `u` the unit axis and `v` the atom's position relative to `fixed`.
    Every dihedral X-I-J-Y across the bond then changes by `degrees`, and by
    the same angle whichever side is the moving one."""
    origin = atoms[fixed][1:]
    line = tuple(atoms[mover][k] - origin[k - 1] for k in (1, 2, 3))
    length = math.sqrt(sum(component * component for component in line))
    u = tuple(component / length for component in line)
    cos_t = math.cos(math.radians(degrees))
    sin_t = math.sin(math.radians(degrees))
    moving = set(moved)
    turned = []
    for index, atom in enumerate(atoms):
        if index not in moving:
            turned.append(atom)
            continue
        v = tuple(atom[k] - origin[k - 1] for k in (1, 2, 3))
        cross = (u[1] * v[2] - u[2] * v[1],
                 u[2] * v[0] - u[0] * v[2],
                 u[0] * v[1] - u[1] * v[0])
        dot = u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
        turned.append((atom[0],) + tuple(
            origin[k] + v[k] * cos_t + cross[k] * sin_t + u[k] * dot * (1.0 - cos_t)
            for k in range(3)))
    return turned


def _contacts(before_pairs: list, after: list, moved: list, edited: tuple) -> tuple[list, tuple | None]:
    """The bonds the motion formed between the two sides, and the closest new
    contact between them (D1).

    `new_bonds` are the pairs of `bonded_pairs(after)` that were not bonded
    before, are not the edited bond, and have exactly one atom in `moved`: a
    pair inside the moving side was already there. `closest` is the smallest
    distance over every moved/unmoved pair that is not already a bond,
    whether or not the bond rule calls it one, because a near miss is what
    says the edit is about to make a bond. The edited bond's own two atoms
    are excluded from `closest`: each must stay near the other, and that pair
    is not new."""
    known = set(before_pairs)
    moving = set(moved)
    new_bonds = []
    for a, b in bonded_pairs(after):
        if (a, b) in known or (a, b) == edited or (a in moving) == (b in moving):
            continue
        new_bonds.append((a, b, math.dist(after[a][1:], after[b][1:])))
    closest = None
    for a in moving:
        for b in range(len(after)):
            if b in moving or a in edited or b in edited:
                continue
            key = (min(a, b), max(a, b))
            if key in known:
                continue
            d = math.dist(after[a][1:], after[b][1:])
            if closest is None or d < closest[0]:
                closest = (d, key[0], key[1])
    return sorted(new_bonds), closest


def rotate_bond(atoms: list, i: int, j: int, degrees: float, move: int | None = None) -> EditResult:
    """Turn one side of the bond i-j by `degrees` (D4): the axis is fixed ->
    mover and the right-hand rule turns EVERY dihedral X-i-j-Y (X on i's
    side, Y on j's side) by the same `degrees`, whichever side moves.

    The input list is left alone; the new one is in the result. The fixed and
    moving sides are D2's (the smaller moves; `move=ATOM` picks the other),
    and `before`/`after` are None -- a rotation aims at no value."""
    pairs = bonded_pairs(atoms)
    side_i, side_j = _sides(atoms, pairs, i, j, True)
    fixed, mover, moved = _choose(side_i, side_j, i, j, move)
    after = _rotate(atoms, moved, fixed, mover, degrees)
    fixed_side = len(side_i) if fixed == i else len(side_j)
    new_bonds, closest = _contacts(pairs, after, moved, (min(i, j), max(i, j)))
    return EditResult(atoms=after, moved=sorted(moved), fixed=fixed, mover=mover,
                      fixed_side=fixed_side, new_bonds=new_bonds, closest=closest)


def set_dihedral(atoms: list, a: int, b: int, c: int, d: int, target: float,
                 move: int | None = None) -> EditResult:
    """Turn the bond b-c until the dihedral a-b-c-d is `target` degrees (D4),
    reporting the value it started from (`before`) and reached (`after`).

    The angle given to `rotate_bond` is the difference between target and
    current, wrapped into [-180, 180), so every dihedral across b-c turns
    with it and one choice reaches the target exactly. `a` must lie on b's
    side of the bond and `d` on c's; anything else names a dihedral that is
    not across b-c. The input list is left alone."""
    pairs = bonded_pairs(atoms)
    side_b, side_c = _sides(atoms, pairs, b, c, True)
    if a not in side_b or d not in side_c:
        raise EditError(f"atoms {a} and {d} must lie on opposite sides of the bond {b}-{c}")
    current = measure(atoms, [a, b, c, d])[1]
    delta = ((target - current + 180.0) % 360.0) - 180.0
    result = rotate_bond(atoms, b, c, delta, move)
    result.before = current
    result.after = measure(result.atoms, [a, b, c, d])[1]
    return result


def set_distance(atoms: list, i: int, j: int, target: float, move: int | None = None) -> EditResult:
    """Slide one side along the line from the fixed atom to the moving one
    until i and j are `target` angstrom apart (D4), reporting the distance it
    started from (`before`) and reached (`after`).

    i and j must be bonded (and not in a ring) or lie in separate fragments:
    two atoms of one fragment that are not bonded are refused, because
    sliding one side would bend the fragment rather than move it. The input
    list is left alone."""
    if target <= 0:
        raise EditError("a distance must be positive")
    pairs = bonded_pairs(atoms)
    side_i, side_j = _sides(atoms, pairs, i, j, False)
    fixed, mover, moved = _choose(side_i, side_j, i, j, move)
    origin = atoms[fixed][1:]
    line = tuple(atoms[mover][k] - origin[k - 1] for k in (1, 2, 3))
    current = math.dist(atoms[i][1:], atoms[j][1:])
    step = (target - current) / current
    moving = set(moved)
    pulled = []
    for index, atom in enumerate(atoms):
        if index not in moving:
            pulled.append(atom)
            continue
        pulled.append((atom[0],) + tuple(
            atom[k] + line[k - 1] * step for k in (1, 2, 3)))
    fixed_side = len(side_i) if fixed == i else len(side_j)
    new_bonds, closest = _contacts(pairs, pulled, moved, (min(i, j), max(i, j)))
    return EditResult(atoms=pulled, moved=sorted(moved), fixed=fixed, mover=mover,
                      fixed_side=fixed_side, new_bonds=new_bonds, closest=closest,
                      before=current, after=math.dist(pulled[i][1:], pulled[j][1:]))


def bond_changes(before_atoms: list, after_atoms: list) -> tuple[list, list, int]:
    """`(formed, broken, bonds after)` between two geometries of one size,
    under the same D3 bond rule an edit uses.

    Each entry is `(a, b, distance_before, distance_after)` with `a < b`, and
    both lists are sorted. The shared rule is the point: "what did the run do
    to the bonds" and "what will this edit make" are one answer, not two
    that can drift apart. ORCA may print different parts of a multilayer
    system at different times, so a size mismatch is refused rather than
    compared atom for atom."""
    if len(before_atoms) != len(after_atoms):
        raise EditError(f"the two geometries have {len(before_atoms)} and {len(after_atoms)} atoms: they cannot be compared")
    before = bonded_pairs(before_atoms)
    after = bonded_pairs(after_atoms)
    before_set, after_set = set(before), set(after)
    formed = [(a, b, math.dist(before_atoms[a][1:], before_atoms[b][1:]),
               math.dist(after_atoms[a][1:], after_atoms[b][1:]))
              for a, b in after if (a, b) not in before_set]
    broken = [(a, b, math.dist(before_atoms[a][1:], before_atoms[b][1:]),
               math.dist(after_atoms[a][1:], after_atoms[b][1:]))
              for a, b in before if (a, b) not in after_set]
    return sorted(formed), sorted(broken), len(after)
