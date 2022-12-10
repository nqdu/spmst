#!/bin/bash
set -e 

file=models/mod_iter10.dat
line=`gmt gmtinfo $file -C`
xmin=`echo $line |awk '{print $1}'`
xmax=`echo $line |awk '{print $2}'`
zmin=`echo $line |awk '{print $3}'`
zmax=`echo $line |awk '{print $4}'`

bounds=-R$xmin/$xmax/$zmin/$zmax
proj=-JM12c
echo $bounds

# generate ray density
if [[ ! -f "white.dat" ]]; then
    python get_raydensity.py $xmin $xmax $zmin $zmax 1
fi

awk '{print $1,$2,$3}' $file | gmt surface $bounds -I128+n/128+n -Vq  -Gtmp.grd
gmt xyz2grd white.dat $bounds -I128+n/128+n -Vq  -Gwhite.grd
gmt grdmath white.grd tmp.grd MUL = model.grd 

# colorbar 
vmin=`gmt grdinfo -C model.grd | awk '{print $6}'`
vmax=`gmt grdinfo -C model.grd | awk '{print $7}'`
gmt makecpt -Cvik -I -D -Z -T$vmin/$vmax/100+n > out.cpt 

# plot
gmt begin out pdf 
gmt basemap $bounds $proj -Bxaf -Byaf -BWnSe 
gmt grdimage model.grd -Cout.cpt  $bounds $proj -E200
#grep -v '^#' surfdata.txt | awk '{print $1,$2}' |gmt plot -St0.15c -Gblue
grep '^#' surfdata.txt | awk '{print $2,$3}' |gmt plot -St0.15c -Gblue
gmt colorbar $bounds $proj -Cout.cpt -Bxaf+l"veloc,km/s"
gmt end 

rm  *.grd

