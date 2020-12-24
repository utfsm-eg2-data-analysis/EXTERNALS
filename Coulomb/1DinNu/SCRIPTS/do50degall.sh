#!/bin/tcsh
cd ../
nice +10 run_extern 50degh2
cat SCRIPTS/50degheader.txt OUT/50degh2.out >! OUT/5.766_50.0_1.out
nice +10 run_extern 50degd2
cat SCRIPTS/50degheader.txt OUT/50degd2.out >! OUT/5.766_50.0_2.out
#nice +10 run_extern 50deghe3
#cat SCRIPTS/50degheader.txt OUT/50deghe3.out >! OUT/5.766_50.0_3.out
#nice +10 run_extern 50deghe4
#cat SCRIPTS/50degheader.txt OUT/50deghe4.out >! OUT/5.766_50.0_4.out
#nice +10 run_extern 50degbe
#cat SCRIPTS/50degheader.txt OUT/50degbe.out >! OUT/5.766_50.0_9.out
#nice +10 run_extern 50degc12
#cat SCRIPTS/50degheader.txt OUT/50degc12.out >! OUT/5.766_50.0_12.out
#nice +10 run_extern 50degcu
#cat SCRIPTS/50degheader.txt OUT/50degcu.out >! OUT/5.766_50.0_63.out
#nice +10 run_extern 50degau
#cat SCRIPTS/50degheader.txt OUT/50degau.out >! OUT/5.766_50.0_197.out
