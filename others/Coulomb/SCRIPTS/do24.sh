#!/bin/tcsh
cd ../
nice +10 run_extern 24degc12
cat SCRIPTS/24degheader.txt OUT/24degc12.out >! OUT/5.011_24.0_12.out
nice +10 run_extern 24degd2
cat SCRIPTS/24degheader.txt OUT/24degd2.out >! OUT/5.011_24.0_2.out
