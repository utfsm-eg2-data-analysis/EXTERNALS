#!/bin/tcsh
cp dummy.inp extern.inp
nice +10 externals_all
cat header.txt extern.out >! OUT/5.766_40.0_dummy.out
cp h2walls.inp extern.inp
nice +10 externals_all
cat header.txt extern.out >! OUT/5.766_40.0_h2walls.out
cp d2walls.inp extern.inp
nice +10 externals_all
cat header.txt extern.out >! OUT/5.766_40.0_d2walls.out


