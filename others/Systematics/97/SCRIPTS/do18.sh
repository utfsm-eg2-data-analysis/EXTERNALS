#!/bin/tcsh
cd ../
nice +10 run_extern 18degau
cat SCRIPTS/18degheader.txt OUT/18degau.out >! OUT/5.766_18.0_197.out
nice +10 run_extern 18degcu
cat SCRIPTS/18degheader.txt OUT/18degcu.out >! OUT/5.766_18.0_63.out
nice +10 run_extern 18degc12
cat SCRIPTS/18degheader.txt OUT/18degc12.out >! OUT/5.766_18.0_12.out
nice +10 run_extern 18degbe
cat SCRIPTS/18degheader.txt OUT/18degbe.out >! OUT/5.766_18.0_9.out
nice +10 run_extern 18deghe4
cat SCRIPTS/18degheader.txt OUT/18deghe4.out >! OUT/5.766_18.0_4.out
nice +10 run_extern 18deghe3
cat SCRIPTS/18degheader.txt OUT/18deghe3.out >! OUT/5.766_18.0_3.out
nice +10 run_extern 18degd2
cat SCRIPTS/18degheader.txt OUT/18degd2.out >! OUT/5.766_18.0_2.out
nice +10 run_extern 18degh2
cat SCRIPTS/18degheader.txt OUT/18degh2.out >! OUT/5.766_18.0_1.out
