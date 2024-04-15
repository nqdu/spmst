#!/bin/bash 
set -e

# get travel time grd
file=time.out
line=`gmt gmtinfo $file -C`
xmin=`echo $line |awk '{print $1}'`
xmax=`echo $line |awk '{print $2}'`
zmin=`echo $line |awk '{print $3}'`
zmax=`echo $line |awk '{print $4}'`
bounds=-R$xmin/$xmax/$zmin/$zmax
proj=-JX12c
awk '{print $1,$2,$4}' $file | gmt surface $bounds -I128+n/128+n -Vq  -Gtime.grd
awk '{print $1,$2,$4}' time.true.out | gmt surface $bounds -I128+n/128+n -Vq   -Gtime.true.grd
gmt grdmath time.grd time.true.grd SUB time.true.grd DIV 100 MUL = error.grd
#gmt grdmath time.grd time.true.grd SUB = error.grd
gmt grdinfo -C *.grd

# read topo
nlon=`head -1 topo.txt|awk '{print $1}'`
nlat=`head -1 topo.txt |awk '{print $2}'`
echo -I$nlon+n/$nlat+n
sed -n '3,$p' topo.txt  > topo1.dat 
gmt xyz2grd topo1.dat $bounds -I$nlon+n/$nlat+n -ZBLa  -Gtopo.grd

# colorbar 
vmin=`gmt grdinfo -C topo.grd | awk '{print $6}'`
vmax=`gmt grdinfo -C topo.grd | awk '{print $7}'`
gmt makecpt -Ctopo  -D -Z -T$vmin/$vmax/100+n > out.cpt 
vmin=`gmt grdinfo -C error.grd | awk '{print $6}'`
vmax=`gmt grdinfo -C error.grd | awk '{print $7}'`
gmt makecpt -Cgray -I -D -Z -T0/$vmax/100+n > error.cpt 
 
# plot
gmt begin out jpg

gmt basemap $bounds $proj -Bxaf+l"x,km" -Byaf+l"y,km" -BwnSe
gmt grdimage error.grd -Cerror.cpt  $bounds $proj -E200
gmt colorbar $bounds $proj -I -Cerror.cpt -Bxaf+l"Error,%" 

gmt basemap $bounds $proj -Bxaf+l"x,km" -Byaf+l"y,km" -BWnSe  -X14c 
gmt grdimage topo.grd -Cout.cpt  $bounds $proj -E200
gmt grdcontour time.grd $bounds $proj  -W0.75p,black -l"Topo"
gmt grdcontour time.true.grd $bounds $proj  -W0.75p,black,. -l"No Topo"
grep -v '^#' surfdata.txt | awk '{print $1,$2}' |gmt plot -St0.2c -Gblack
grep '^#' surfdata.txt | awk '{print $2,$3}' |gmt plot -Sa0.2c -Gblack
gmt plot ray.dat -W0.4p,black
gmt plot ray.true.dat -W0.4p,black,.
gmt colorbar $bounds $proj -I -Cout.cpt -Bxaf+l"Topography,m" 

gmt end 

rm  *.grd *.cpt 

# error bar

