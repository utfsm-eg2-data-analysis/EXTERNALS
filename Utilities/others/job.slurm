#!/bin/bash
#SBATCH --time=00:40:05
#SBATCH --mail-user=sebastian.moran@alumnos.usm.cl
#SBATCH --mail-type=FAIL
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

met_opt=${opt:0:2}    
tar_opt=${opt:3:3}     ## this takes values: "RD" or "HH"
vc_opt=${opt:7:2}    ## this takes values: "RD" or "HH"


if [ "$met_opt" = "Ca" ] ; then
	Metal="C"
elif [ "$met_opt" = "Fe" ] ; then
	Metal="Fe"
elif [ "$met_opt" = "Pb" ] ; then
	Metal="Pb"
fi

cd $TMPDIR
sourceDir="/work/smoran/externals/centroids"
exe="Centroid_AC_RC"
outdir="/work/smoran/externals/centroids/w_20Q2_10_11x12GridAC_RC_Rad_DIS"

cp $sourceDir/$exe exe
if [ ! -d $outdir ]; then mkdir -p $outdir; fi

./exe $Metal $tar_opt $vc_opt
cp *.root $outdir
cp *.txt $outdir



















