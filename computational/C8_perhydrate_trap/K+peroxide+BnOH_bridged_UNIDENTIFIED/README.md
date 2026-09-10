# Not K + HOO- -> KP. Renamed 2026-09-09; true identity unresolved.

Originally `orca_stuff/cat/catHO2(-)+BnOH/`, migrated here as C8's
`K + HOO- -> KP` (item 1: the perhydrate-formation step). It isn't that.

Direct connectivity analysis of both `geometry_r2scan3c-xtb/rc/job.xyz` and
`.../ts/job.xyz` (bond distances, not guessed) shows the peroxide unit
bridging TWO different carbons at once, not adding to the catalyst's ketone
alone:

    O131~C129 ---(1.32/1.56 A, RC/TS, elongating)--- O130 --(1.45 A)-- O132 ---(1.53/1.43 A)--- C138
    [catalyst ketone carbon]                    [peroxide]                    [BnOH benzylic carbon]

`C129` is the catalyst's own ketone carbon (bonded to two ring carbons plus
this system's O131/O130). `C138` is BnOH's benzylic carbon -- confirmed by
tracing an intact monosubstituted phenyl ring off it (`C139`) and its own
original hydroxyl (`O137-H134`) -- missing one of its original two benzylic
hydrogens, replaced by the new bond to `O132`.

So this is a three-body, concerted-looking TS: catalyst ketone ... peroxide
... BnOH benzylic carbon, all in one chain, not a simple two-body
`K + HOO- <=> KP` addition. It most likely belongs with the OTHER family
`orca_stuff/cat/cat+BnOH+H2O2/` was flagged for early on (COMPUTATIONAL.md
mapping, 2026-09-09) -- "a direct three-body TS, possibly skipping the
discrete KP intermediate" -- but that hasn't been confirmed, and neither has
which mechanism step (5, 7, or something not yet in `MECHANISM.md`) this
would correspond to.

**Do not spend more compute on this folder under the C8 framing.** It needs
its own mechanistic scoping (what bond is actually forming/breaking, is it
really concerted, does it connect to `cat+BnOH+H2O2/`) before any method
work continues. The real C8-item-1 pilot is
`../K+H2O2_water-relay_to_KP/`.
