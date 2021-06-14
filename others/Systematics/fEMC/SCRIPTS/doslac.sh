#!/bin/tcsh
cd ../
nice +10 run_extern slac_e139_cu
cat SCRIPTS/40degheader.txt OUT/slac_e139_cu.out >! OUT/slac_e139_cu6.out
