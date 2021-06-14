#!/bin/tcsh
cd ../
nice +10 run_extern seamus_pvdis
cat SCRIPTS/18degheader.txt OUT/seamus_pvdis.out >! OUT/11.00_35.0_48.out
