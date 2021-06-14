#!/bin/bash
#PBS -M sebastian.moran@alumnos.usm.cl
#PBS -N 48
#PBS -l walltime=50:00:00 
#PBS -o /dev/null # no stdout
#PBS -e /dev/null # no stderr
#PBS -q utfsm 
cd $TMPDIR    
echo working on `pwd`
metal="Pb"
targtype="sol"
vc="OS_f"

fout_root="/user/s/smoran/centroid/con_ceros/${metal}_${targtype}_${vc}_vc_centroid.root"
fout_txt1="/user/s/smoran/centroid/con_ceros/${metal}_${targtype}_${vc}_vc_centroid.txt"
fout_txt2="/user/s/smoran/centroid/con_ceros/${metal}_${targtype}_${vc}_vc_centroid_for_plot.txt"

cp /user/s/smoran/centroid/Centroid .
cp /data/user/s/smoran/elec_data/${metal}_data.root .
./Centroid ${metal} ${targtype} ${vc}

cp centroids_${metal}_${targtype}_${vc}_vc.txt $fout_txt1
cp centroids_${metal}_${targtype}_${vc}_vc_for_plot.txt  $fout_txt2
cp out_${metal}_${targtype}_${vc}_vc.root $fout_root
