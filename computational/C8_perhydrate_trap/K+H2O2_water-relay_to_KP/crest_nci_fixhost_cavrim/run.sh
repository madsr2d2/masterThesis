#!/bin/bash
# Host frozen, bias on the 54 water atoms; 5 ps MTDs (the default 30 ps cap took hours on 277 atoms).
cd "$(dirname "$0")"
/home/madsr2d2/crest-3.0.2/crest solvcluster.xyz -nci --gfnff --alpb water \
    -T 8 -mquick -mdlen 5 -cinp fixhost.inp > crest.out 2>&1
