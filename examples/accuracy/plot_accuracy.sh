#!/bin/bash 
set -e 
gmt set FONT_TITLE 16p,0,Black
gmt set FONT_LABEL 16p,0,Black
gmt set FONT_ANNOT 16p,0

# get travel time grd
file=time.$1.out
file_50=time.true.out
line=`gmt gmtinfo $file -C`
xmin=`echo $line |awk '{print $1}'`
xmax=`echo $line |awk '{print $2}'`
zmin=`echo $line |awk '{print $3}'`
zmax=`echo $line |awk '{print $4}'`
bounds=-R$xmin/$xmax/$zmin/$zmax
proj=-JX12c
echo $bounds
# exit

awk '{print $1,$2,$4}' $file | gmt surface $bounds -I128+n/128+n -Vq  -Gtime_$1.grd
awk '{print $1,$2,$4}' $file_50 | gmt surface $bounds -I128+n/128+n -Vq  -Gtime_50.grd
#awk '{print $1,$2,$4}' time.out.notopo | gmt surface $bounds -I128+n/128+n -Vq  -Gnotop.grd

# read topo
nlon=`head -1 topo.txt |awk '{print $1}'`
nlat=`head -1 topo.txt |awk '{print $2}'`
echo -I$nlon+n/$nlat+n
sed -n '3,$p' topo.txt  > topo_tmp.dat 
gmt xyz2grd topo_tmp.dat $bounds -I$nlon+n/$nlat+n -ZBLa  -Gtopo.grd

# colorbar 
vmin=`gmt grdinfo -C topo.grd | awk '{print $6}'`
vmax=`gmt grdinfo -C topo.grd | awk '{print $7}'`
gmt makecpt -Cdem2 -I -D -Z -T$vmin/$vmax/100+n > out.cpt 

 
gmt grdmath time_$1.grd time_50.grd SUB time_50.grd DIV 100 MUL = time_error.grd

emin=`gmt grdinfo -C time_error.grd | awk '{print $6}'`
emax=`gmt grdinfo -C time_error.grd | awk '{print $7}'`

# echo "($emin*-1+$emax)/2" | bc
# error_cpt=`echo "($emin*-1+$emax)/2" | bc`

#gmt makecpt -Cpolar  -D -Z -T$emin/$emax/100+n > error.cpt 
gmt makecpt -Cgray -I -D -Z -T0/1/100+n > error.cpt 

# plot
gmt begin error_$1 png,pdf

gmt basemap $bounds $proj -Bxaf+l"x(km)"  -Byaf+l"y(km)"  -BWnSe  -X12c -Y13c
gmt grdimage time_error.grd -Cerror.cpt  $bounds -E200  
gmt colorbar -Cerror.cpt -Bxa1f1+l"Relative Error(%)"
echo "(a)" | gmt text -F+f16p,0,black+cTL -Dj-0.7i/-0.1i -N

gmt basemap $bounds $proj -Bxaf+l"x(km)"  -Byaf+l"y(km)"  -BWnSe  -X14c
gmt grdimage topo.grd -Cout.cpt  $bounds -E200 -t25
gmt grdcontour time_$1.grd $bounds -W0.5p,black -l"SPM Time" -An -t25
gmt grdcontour time_50.grd $bounds -W1p,darkred,. -l"True Time"
#gmt grdcontour notop.grd $bounds $proj  -W0.75p,black,. -l"No Topo"
grep -v '^#' surfdata.txt | awk '{print $1,$2}' |gmt plot -St0.3c -W0.5p,black -Gblack
grep '^#' surfdata.txt | awk '{print $2,$3}' |gmt plot -Sa0.4c -Gyellow -W0.5p,black
gmt plot ray.$1.dat -W0.5p,black -t35 -l"SPM Rays" -An
gmt plot ray.true.dat -W1p,darkblue,- -l"True Rays" -An
gmt colorbar -Cout.cpt -Bxaf+l"Topography(m)" 
echo "(b)" | gmt text -F+f16p,0,black+cTL -Dj-0.7i/-0.1i -N
gmt end 

rm  *.grd *.cpt  gmt.* topo_tmp.dat