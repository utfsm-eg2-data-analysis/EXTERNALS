#!/bin/sh 
echo 'starting'
mv centroids_C_sol_RD_vc_for_plot.txt centroids_C_st_RD_vc_for_plot.txt
mv centroids_C_liq_RD_vc_for_plot.txt centroids_C_lt_RD_vc_for_plot.txt
mv centroids_Fe_sol_RD_vc_for_plot.txt centroids_Fe_st_RD_vc_for_plot.txt
mv centroids_Fe_liq_RD_vc_for_plot.txt centroids_Fe_lt_RD_vc_for_plot.txt
mv centroids_Pb_sol_RD_vc_for_plot.txt centroids_Pb_st_RD_vc_for_plot.txt
mv centroids_Pb_liq_RD_vc_for_plot.txt centroids_Pb_lt_RD_vc_for_plot.txt

sed -i '1s/^/first line\n/' centroids_C_st_RD_vc_for_plot.txt
sed -i '1s/^/first line\n/' centroids_C_lt_RD_vc_for_plot.txt
sed -i '1s/^/first line\n/' centroids_Fe_st_RD_vc_for_plot.txt
sed -i '1s/^/first line\n/' centroids_Fe_lt_RD_vc_for_plot.txt
sed -i '1s/^/first line\n/' centroids_Pb_st_RD_vc_for_plot.txt
sed -i '1s/^/first line\n/' centroids_Pb_lt_RD_vc_for_plot.txt




