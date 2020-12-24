#!/bin/tcsh
cd ../
run_extern 36degc12
cat SCRIPTS/36degheader.txt OUT/36degc12.out >! OUT/5.011_36.0_12.out
run_extern 36degd2
cat SCRIPTS/36degheader.txt OUT/36degd2.out >! OUT/5.011_36.0_2.out
