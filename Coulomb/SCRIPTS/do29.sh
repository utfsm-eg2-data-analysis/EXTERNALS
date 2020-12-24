#!/bin/tcsh
cd ../
nice +10 run_extern 29degc12
cat SCRIPTS/29degheader.txt OUT/29degc12.out >! OUT/5.011_29.0_12.out
nice +10 run_extern 29degd2
cat SCRIPTS/29degheader.txt OUT/29degd2.out >! OUT/5.011_29.0_2.out
