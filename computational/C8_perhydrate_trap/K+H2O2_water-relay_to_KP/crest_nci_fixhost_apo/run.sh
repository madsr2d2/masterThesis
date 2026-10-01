#!/bin/bash
# Frozen-host NCI search: the 130 host atoms are restrained to the seed, the
# metadynamics bias acts on the 6 water atoms (131-136) only. See README.md.
cd "$(dirname "$0")"
/home/madsr2d2/crest-3.0.2/crest solvcluster.xyz -nci --gfnff --alpb water \
    -T 8 -cinp fixhost.inp > crest.out 2>&1
