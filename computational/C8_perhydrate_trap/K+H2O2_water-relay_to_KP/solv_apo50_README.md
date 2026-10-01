# solv_apo50_{a,b,c} -- three stochastic 50-water droplets around the bare host

Same input as `solv_apo6` (`! XTB ALPB(Water) SOLVATOR`, `CLUSTERMODE
STOCHASTIC`) with `NSOLVENTMOL 50`. The 130-atom host is read from
`../solv_apo6/job.xyz`, not copied, so the three runs share one starting
geometry. Three independent runs because the mode is random and the 3- and
6-water runs gave 1 and 2 cavity waters.

Next: count the cavity waters in each droplet (`crest_io.cavity_occupancy`),
trim to the waters within ~4 A of the host, then CREST with the host frozen
(`crest -constrain 1-130`, bias on the waters). The cycle this is for:
`cavity.(H2O)n + S.(H2O)m <=> cavity.S + (H2O)n.(H2O)m`, where n is the number
of waters the empty cavity holds at the free-energy minimum.
