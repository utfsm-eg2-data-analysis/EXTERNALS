#!/bin/tcsh
cd ../
run_extern 46degau
cat SCRIPTS/46degheader.txt OUT/46degau.out >! OUT/5.011_46.0_197.out
run_extern 46degcu
cat SCRIPTS/46degheader.txt OUT/46degcu.out >! OUT/5.011_46.0_63.out
run_extern 46degc12
cat SCRIPTS/46degheader.txt OUT/46degc12.out >! OUT/5.011_46.0_12.out
run_extern 46degbe
cat SCRIPTS/46degheader.txt OUT/46degbe.out >! OUT/5.011_46.0_9.out
run_extern 46degd2
cat SCRIPTS/46degheader.txt OUT/46degd2.out >! OUT/5.011_46.0_2.out
