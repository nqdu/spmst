#!/bin/bash

if [ $# -ne 2 ];then 
    echo "Usage: ./plot_rel_model.sh model model_init"
    exit 1
fi 

# get file info
file=$1
file0=$2
line=`gmt gmtinfo $file -C`
xmin=`echo $line |awk '{print $1}'`
xmax=`echo $line |awk '{print $2}'`
zmin=`echo $line |awk '{print $3}'`
zmax=`echo $line |awk '{print $4}'`
bounds=-R$xmin/$xmax/$zmin/$zmax
proj=-JM12c
Ix=`echo "$xmin $xmax" |awk '{print ($2-$1)/127.}'`
Iy=`echo "$zmin $zmax" |awk '{print ($2-$1)/127.}'`
gmt surface $file $bounds -I$Ix/$Iy  -Gout.grd -Vq
gmt surface $file0 $bounds -I$Ix/$Iy  -Gout.init.grd -Vq
gmt grdmath out.grd out.init.grd SUB out.init.grd DIV 100 MUL = error.grd

# colorbar
#vmin=`gmt grdinfo -C out.grd | awk '{print $6}'`
#vmax=`gmt grdinfo -C out.grd | awk '{print $7}'`
vmin=-10
vmax=10
gmt makecpt -Cpolar -I -D -Z -T$vmin/$vmax/100+n > out.cpt 

gmt begin out jpg
gmt basemap $bounds $proj -Bxaf -Byaf -BWnSe 
gmt grdimage error.grd -Cout.cpt  $bounds $proj -E200
grep -v '^#' surfdata.txt |awk '{print $1,$2}' | gmt plot -St0.2c -Gblack
gmt colorbar $bounds $proj -I -Cout.cpt -Bxaf+l"@[\Delta c / c@[(\%)" 

gmt end 

\rm out.cpt out.grd out.init.grd error.grd