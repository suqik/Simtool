#!/bin/bash
RSTAR=/public/home/suchen/applications/rockstar/rockstar
CFG=/public/home/suchen/Programs/Simtool/Pipeline/cfgs/cosmo0/rockstar/z0.30.cfg
ODIR=/public/share/ace66so15x/suchen/L1000_N1024_1000cosmo/cosmo0/a_0.7692/rstar

$RSTAR -c $CFG &
export ODIR
perl -e 'sleep 1 while (!(-e "$ENV{ODIR}/auto-rockstar.cfg"))'
$RSTAR -c ${ODIR}/auto-rockstar.cfg