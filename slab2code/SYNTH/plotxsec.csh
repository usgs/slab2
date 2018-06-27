awk -F, '{if(NR>1) print $2,$1,$3}' exp_04-18_input.csv | gmt project -C126.4/1.3 -A105 -L-50/400 -W-350/350 -Q | awk '{print $4,$3}' | gmt psxy -JX15/-10 -R-50/400/0/300 -Sc0.2 -W0.25p -G200 -K -BSnWe -Bx40g20 -By20g10 > temp2.ps
awk '{print $1,$2}' ../library/avprofiles/global_as_av.txt | gmt psxy -J -R -Sc0.2 -G100 -W0.5p -O -K >> temp2.ps
gmt project -C126.4/1.3 -A105 -L-400/400 -G10 -Q | gmt grdtrack -Gexp_guide2.grd | gmt project -C126.4/1.3 -A105 -Fpz -Q | awk '{print $2,-$3}' | gmt psxy -J -R -W2p,blue,- -K -O >> temp2.ps
gmt project -C126.4/1.3 -A105 -L-400/400 -G10 -Q | gmt grdtrack -G../Output/exp_slab2_06.14.18/exp_slab2_dep_06.14.18.grd  | gmt project -C126.4/1.3 -A105 -Fpz -Q | awk '{print $2,-$3}' | gmt psxy -J -R -W2p,red -K -O >> temp2.ps
