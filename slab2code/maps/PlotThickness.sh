# GMT5
# KLH
# 10/09/2019
# Plot thickness of slab2 model output

gmt gmtset PS_PAGE_ORIENTATION landscape

# Input data:
Map="/Users/khaynie/Src/slab2-master2019cat/slab2code/library/forplotting/world_lr.grd"
ThickGrd="/Users/khaynie/Src/slab2/slab2code/Output/exp_slab2_10.09.19/exp_slab2_thk_10.09.19.grd"
clipmask="/Users/khaynie/Src/slab2/slab2code/Output/exp_slab2_10.09.19/exp_slab2_clp_10.09.19.csv"

# output file:
psfile="Thickness.ps"

rm -rf $psfile
rm -rf thk.cpt

# boundsion parameters:
bounds="-R122/134/-6/7"
proj="-JM6i"
ticks="-Ba4g4/a2g2WESN"

# miscellaneous parameters:
miscStart="-K -V"
miscMid="-O -K -V"
miscEnd="-O -V"
#AT="-C2 -A10" # contouring
AT="-C5 -A5" # contouring

# cpt files
Mapcpt="/Users/khaynie/Src/slab2-master2019cat/slab2code/library/forplotting/ghayes2.cpt"
gmt makecpt -Cno_green -T105/125/5 -Z >> thk.cpt

# Plot
gmt psbasemap $bounds $proj $ticks $miscStart > $psfile
gmt grdimage $Map -C$Mapcpt $bounds $proj $miscMid -BSnWe+t"Thickness Grid" >> $psfile
gmt psclip $clipmask $proj $bounds $miscMid >> $psfile
gmt grdimage $ThickGrd -Cthk.cpt $bounds $proj $miscMid >> $psfile
gmt grdcontour $ThickGrd $AT $bounds $proj $miscMid >> $psfile
gmt psxy $clipmask $bounds $proj -W2p,black $miscMid >> $psfile
#gmt pscoast -J -Df -W0.75p,100 $bounds $miscMid >> $psfile
gmt psscale -Cthk.cpt -D0.5/3.0/5/0.5 -B10 $miscEnd >> $psfile

open $psfile
