#!/bin/sh 
echo 'starting'
## starting state 1.5 and 1.8

#./Centroid C sol RD
#./Centroid C liq RD
#./Centroid Fe sol RD
#./Centroid Fe liq RD
#./Centroid Pb sol RD
#./Centroid Pb liq RD

#mkdir w_18Q2_15_35X35
#mv *.root w_18Q2_15_35X35
#mv *.txt w_18Q2_15_35X35

#sed -i 's|limQ2=1.5|limQ2=1.0|g' Centroid.cxx
#make
#./Centroid C sol RD
#./Centroid C liq RD
#./Centroid Fe sol RD
#./Centroid Fe liq RD
#./Centroid Pb sol RD
#./Centroid Pb liq RD

#mkdir w_18Q2_10_35X35
#mv *.root w_18Q2_10_35X35
#mv *.txt w_18Q2_10_35X35

#sed -i 's|limW2=3.24|limW2=4.0|g' Centroid.cxx
#make 
#./Centroid C sol RD
#./Centroid C liq RD
#./Centroid Fe sol RD
#./Centroid Fe liq RD
#./Centroid Pb sol RD
#./Centroid Pb liq RD

#mkdir w_20Q2_10_35X35
#mv *.root w_20Q2_10_35X35
#mv *.txt w_20Q2_10_35X35

#sed -i 's|limQ2=1.0|limQ2=1.5|g' Centroid.cxx
#make 
echo 'Carbon '
echo 'solid '
./Centroid C sol RD
echo 'liquid '
./Centroid C liq RD
echo 'Iron '
echo 'solid '
./Centroid Fe sol RD
echo 'liquid '
./Centroid Fe liq RD
echo 'Lead '
echo 'solid '
./Centroid Pb sol RD
echo 'liquid '
./Centroid Pb liq RD

#./script_changeNames.sh

#    mkdir w_20Q2_10_60x45GridAC
#mv *.root w_20Q2_10_60x45GridAC
#mv *.txt  w_20Q2_10_60x45GridAC


echo 'done! '




