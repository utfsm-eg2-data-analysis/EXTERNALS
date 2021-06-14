#!/bin/bash

###vamos a partir con el archivo hecho para C liq y H vertex cuts, de ahi voy a construir el resto

qsub job.sh
##sed -i 's/targtype="liq"/targtype="sol"/g' job.sh
sed -i 's/#PBS -N 1/#PBS -N 2/g' job.sh
sed -i 's/vc="H"/vc="KPA"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 2/#PBS -N 3/g' job.sh
sed -i 's/vc="KPA"/vc="RD"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 3/#PBS -N 4/g' job.sh
sed -i 's/vc="RD"/vc="OS"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 4/#PBS -N 5/g' job.sh
sed -i 's/vc="OS"/vc="TM"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 5/#PBS -N 6/g' job.sh
sed -i 's/vc="TM"/vc="TM_f"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 6/#PBS -N 7/g' job.sh
sed -i 's/vc="TM_f"/vc="OS_f"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 7/#PBS -N 8/g' job.sh
sed -i 's/vc="OS_f"/vc="RD_f"/g' job.sh
qsub job.sh

## done with carbon , liquid target

sed -i 's/#PBS -N 8/#PBS -N 9/g' job.sh
sed -i 's/targtype="liq"/targtype="sol"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 9/#PBS -N 10/g' job.sh
sed -i 's/vc="RD_f"/vc="OS_f"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 10/#PBS -N 11/g' job.sh
sed -i 's/vc="OS_f"/vc="TM_f"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 11/#PBS -N 12/g' job.sh
sed -i 's/vc="TM_f"/vc="H"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 12/#PBS -N 13/g' job.sh
sed -i 's/vc="H"/vc="KPA"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 13/#PBS -N 14/g' job.sh
sed -i 's/vc="KPA"/vc="OS"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 14/#PBS -N 15/g' job.sh
sed -i 's/vc="OS"/vc="RD"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 15/#PBS -N 16/g' job.sh
sed -i 's/vc="RD"/vc="TM"/g' job.sh
qsub job.sh


echo CARBON READY!
###solid




sed -i 's/#PBS -N 16/#PBS -N 17/g' job.sh
sed -i 's/metal="C"/metal="Fe"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 17/#PBS -N 18/g' job.sh
sed -i 's/vc="TM"/vc="RD"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 18/#PBS -N 19/g' job.sh
sed -i 's/vc="RD"/vc="OS"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 19/#PBS -N 20/g' job.sh
sed -i 's/vc="OS"/vc="KPA"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 20/#PBS -N 21/g' job.sh
sed -i 's/vc="KPA"/vc="H"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 21/#PBS -N 22/g' job.sh
sed -i 's/vc="H"/vc="TM_f"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 22/#PBS -N 23/g' job.sh
sed -i 's/vc="TM_f"/vc="RD_f"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 23/#PBS -N 24/g' job.sh
sed -i 's/vc="RD_f"/vc="OS_f"/g' job.sh
qsub job.sh
## now the liquid target

sed -i 's/#PBS -N 24/#PBS -N 25/g' job.sh
sed -i 's/targtype="sol"/targtype="liq"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 25/#PBS -N 26/g' job.sh
sed -i 's/vc="OS_f"/vc="RD"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 26/#PBS -N 27/g' job.sh
sed -i 's/vc="RD"/vc="KPA"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 27/#PBS -N 28/g' job.sh
sed -i 's/vc="KPA"/vc="H"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 28/#PBS -N 29/g' job.sh
sed -i 's/vc="H"/vc="TM"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 29/#PBS -N 30/g' job.sh
sed -i 's/vc="TM"/vc="TM_f"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 30/#PBS -N 31/g' job.sh
sed -i 's/vc="TM_f"/vc="RD_f"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 31/#PBS -N 32/g' job.sh
sed -i 's/vc="RD_f"/vc="OS"/g' job.sh
qsub job.sh


echo IRON READY!

###here the job is for Fe liq H
sed -i 's/#PBS -N 32/#PBS -N 33/g' job.sh
sed -i 's/metal="Fe"/metal="Pb"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 33/#PBS -N 34/g' job.sh
sed -i 's/vc="OS"/vc="KPA"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 34/#PBS -N 35/g' job.sh
sed -i 's/vc="KPA"/vc="RD"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 35/#PBS -N 36/g' job.sh
sed -i 's/vc="RD"/vc="H"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 36/#PBS -N 37/g' job.sh
sed -i 's/vc="H"/vc="OS_f"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 37/#PBS -N 38/g' job.sh
sed -i 's/vc="OS_f"/vc="RD_f"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 38/#PBS -N 39/g' job.sh
sed -i 's/vc="RD_f"/vc="TM"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 39/#PBS -N 40/g' job.sh
sed -i 's/vc="TM"/vc="TM_f"/g' job.sh
qsub job.sh
###solid target
sed -i 's/#PBS -N 40/#PBS -N 41/g' job.sh
sed -i 's/targtype="liq"/targtype="sol"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 41/#PBS -N 42/g' job.sh
sed -i 's/vc="TM_f"/vc="OS"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 42/#PBS -N 43/g' job.sh
sed -i 's/vc="OS"/vc="RD"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 43/#PBS -N 44/g' job.sh
sed -i 's/vc="RD"/vc="H"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 44/#PBS -N 45/g' job.sh
sed -i 's/vc="H"/vc="KPA"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 45/#PBS -N 46/g' job.sh
sed -i 's/vc="KPA"/vc="TM"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 46/#PBS -N 47/g' job.sh
sed -i 's/vc="TM"/vc="RD_f"/g' job.sh
qsub job.sh
sed -i 's/#PBS -N 47/#PBS -N 48/g' job.sh
sed -i 's/vc="RD_f"/vc="OS_f"/g' job.sh
qsub job.sh

echo LEAD READY!
