#!/bin/bash
#set -e 

file=models/mod_iter10.dat
line=`gmt gmtinfo $file -C`
xmin=`echo $line |awk '{print $1}'`
xmax=`echo $line |awk '{print $2}'`
zmin=`echo $line |awk '{print $3}'`
zmax=`echo $line |awk '{print $4}'`

bounds=-R$xmin/$xmax/$zmin/$zmax
proj=-JM12c
echo $bounds

awk '{print $1,$2,$3}' $file | gmt surface $bounds -I128+n/128+n -Vq  -Gmodel.grd
awk '{print $1}' veloctrue.in > tmp.dat 
awk '{print $1,$2}' $file > tmp1.dat 
paste tmp1.dat tmp.dat > tmp.out 
gmt surface tmp.out $bounds -I128+n/128+n -Vq  -Gtrue.grd

# colorbar 
vmin=`gmt grdinfo -C model.grd | awk '{print $6}'`
vmax=`gmt grdinfo -C model.grd | awk '{print $7}'`
gmt makecpt -Cseis -D -Z -T$vmin/$vmax/100+n > out.cpt 

# plot
gmt begin out pdf 
gmt basemap $bounds $proj -Bxaf -Byaf -BWnSe 
gmt grdimage model.grd -Cout.cpt  $bounds $proj -E200
grep -v '^#' surfdata.txt | awk '{print $1,$2}' |gmt plot -St0.3c -Gblue
grep '^#' surfdata.txt | awk '{print $1,$2}' |gmt plot -St0.3c -Gblue
gmt colorbar $bounds $proj -Cout.cpt -Bxaf+l"veloc,km/s"
#awk '{print $1,$2,$3,$5/5}' out_grad.dat| gmt plot $bounds $proj  -Sv0.05i -W0.5p -Gblack 
#awk '{print $1,$2,$4,$5/5}' out_grad.dat| gmt plot $bounds $proj  -Sv0.05i -W0.5p,-- -Gred 

# gmt basemap $bounds $proj -Bxaf -Byaf -BwnSe -X13c
# gmt grdimage true.grd -Cout.cpt  $bounds $proj -E200
# grep -v '^#' surfdata.txt | awk '{print $1,$2}' |gmt plot -St0.3c -Gblue
# gmt colorbar $bounds $proj -Cout.cpt -Bxaf+l"veloc,km/s"
gmt end 

rm  *.grd tmp*

