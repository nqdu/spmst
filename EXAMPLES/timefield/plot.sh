#!/bin/bash 

# get travel time grd
file=time.out
line=`gmt gmtinfo $file -C`
xmin=`echo $line |awk '{print $1}'`
xmax=`echo $line |awk '{print $2}'`
zmin=`echo $line |awk '{print $3}'`
zmax=`echo $line |awk '{print $4}'`
bounds=-R$xmin/$xmax/$zmin/$zmax
proj=-JM12c
echo $bounds
awk '{print $1,$2,$4}' $file | gmt surface $bounds -I128+n/128+n -Vq  -Gtime.grd

# read topo
nlon=`head -1 topo.dat|awk '{print $1}'`
nlat=`head -1 topo.dat |awk '{print $2}'`
echo -I$nlon+n/$nlat+n
sed -n '3,$p' topo.dat  > topo1.dat 
gmt xyz2grd topo1.dat $bounds -I$nlon+n/$nlat+n -ZBLa  -Gtopo.grd

# colorbar 
vmin=`gmt grdinfo -C topo.grd | awk '{print $6}'`
vmax=`gmt grdinfo -C topo.grd | awk '{print $7}'`
gmt makecpt -Crainbow -I -D -Z -T$vmin/$vmax/100+n > out.cpt 

# plot
gmt begin out pdf 
gmt basemap $bounds $proj -Bxaf -Byaf -BWnSe 
gmt grdimage topo.grd -Cout.cpt  $bounds $proj -E200
gmt grdcontour time.grd $bounds $proj  -W0.75p,black
grep -v '^#' surfdata.txt | awk '{print $1,$2}' |gmt plot -St0.2c -Gblack
grep '^#' surfdata.txt | awk '{print $2,$3}' |gmt plot -Sa0.2c -Gred
#gmt plot ray.dat -W0.2p,red
gmt colorbar $bounds $proj -Cout.cpt -Bxaf+l"Topography,m" 
gmt end 

rm  *.grd *.cpt 

