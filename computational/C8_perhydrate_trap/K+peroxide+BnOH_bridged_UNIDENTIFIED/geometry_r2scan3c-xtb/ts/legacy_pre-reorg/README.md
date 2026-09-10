# Pre-2026-09-09 material, kept but not used

These files predate the `computational/` layout and were inherited as-is from
`orca_stuff/cat/catHO2(-)+BnOH/TS/`. They record an IRC (`job_IRC_F.xyz` E =
-278.497220987, `job_IRC_B.xyz` E = -278.500394747) and several optimisation
snapshots (`job.out.v000/v001/v006.xyz`, the `.scfgrad.inp` restarts), so a TS
was apparently found and characterised at some point.

Their energies (~-278.497 to -278.500) do not fall on the energy profile of
`../scan/` (~-278.526 to -278.558), so they belong to a different attempt --
a different QM region, constraint, or reaction coordinate -- and no `.out` log
survives to reconstruct which. Rather than guess, today's TS search
(`../job.inp`) starts fresh from `../scan/job.009.xyz`, the scan's own energy
maximum. If the two converge to the same TS, this material becomes a useful
independent cross-check; if not, it's most likely stale and can be retired
once that's confirmed.
