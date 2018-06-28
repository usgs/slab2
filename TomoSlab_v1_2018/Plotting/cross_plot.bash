#!/bin/bash

### Created DEP.6.15.15 ###
### Edited DEP.5.26.16 ###
### Reworked DEP.7.28.16 to march across trench producing cross sections

###########################

out="sample_cross.eps"
tomo="../Tomo_files/SAM4_P_2017.tomo"
slabcent="../Output/sam_TO_SAM4.csv"
slabfile="../Slab_files/sam_slab1.0_augmented.grd"

rm $out

lat=-32
lon=-72
az=10

cptfile="tomo.cpt"

##########################

gmtset PS_MEDIA=b5
gmtset FONT_ANNOT_PRIMARY=8p,Helvetica,black
gmtset FONT_LABEL=10p,Helvetica,black

minx=-200
maxx=1400
minz=0
maxz=800
tot=$(($maxx - $minx))
L="-L$minx/$maxx"
G="-G$tot"
R="-R$minx/$maxx/$minz/$maxz"

A="-A$az"
az=`echo "$az" | awk '{print $1 + 90}'`
A2="-A$az"
C="-C$lon/$lat"

J="-JX4i/-2i"
F="-Fpz"
X="-X1i"
Y="-Y1i"
W="-W-5/5"
I="-I0.5/0.5"
w="-W-10/10"

project $tomo -: $C $A $L $W $F -N -Q | awk '{print $1,$2,$3}' | blockmean $R $I | surface $R $I -Gprofile.grd -T0

grdimage profile.grd $R $J -C$cptfile -Ba200f100g100:"X Distance (km)":/a100f1000g100:"Depth (km)":WSen -Q -nn -P --MAP_GRID_CROSS_SIZE_PRIMARY=0.05i $X $Y -U -K >> $out
echo "$minx $maxx" | awk '{for (i=$1 ; i<=$2+0.1 ; i+=0.1) print i,410} ; { print ">" } ; {for (i=$1 ; i<=$2+0.1 ; i+=0.1 ) print i,660 }' | gmt psxy -R -J -W1p,black,- -O -K >> $out
gmt grdcontour profile.grd -J -R -C+0.25 -W0.25,darkgray -O -K >> $out
gmt grdcontour profile.grd -J -R -C+0.5 -W0.5,darkgray -O -K >> $out
gmt grdcontour profile.grd -J -R -C+0.75 -W0.75,darkgray -O -K >> $out
gmt grdcontour profile.grd -J -R -C+1.0 -W1.0,darkgray -O -K >> $out
gmt project $C $A2 $L -G10 -N -Q | awk '{print $1+360,$2,$3}' | gmt grdtrack -G$slabfile | awk '{if($4 !="NaN") print $1-360,$2,-1*$4}' | gmt project $C $A $L -N $F -Q | gmt psxy -R -J -W1p,black -O -K >> $out
awk -F',' '{if(NR>1) print $2,$1,$3}' $slabcent | gmt project $C $A $L $F -N -Q $w | gmt psxy -R -J -Sd0.05i -W0.5,black -Gcyan -O -K >> $out

rm profile.grd
