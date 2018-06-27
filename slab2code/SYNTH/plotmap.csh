gmt grdimage /sw/share/gebco/gebco08.grd -C/Users/ghayes/Desktop/CPTs/ghayes2.cpt -R123/132/-4/5 -JM15 -K > temp.ps                      
gmt pscoast -J -R -O -K -Df -W0.6p,50 -BSnWe -Bx4g1 -By4g1 >> temp.ps                                                 
gmt grdimage exp_guide2.grd -C/Users/ghayes/Desktop/CPTs/slab.cpt -O -K -J -R >> temp.ps                                                 
gmt grdcontour -J -R -O -K exp_guide2.grd -C40 -W1p,200 >> temp.ps                                                                       
gmt psxy exp_poly.txt -J -R -O -K -W1.5p,50 >> temp.ps                                                                                   
gmt grdimage ../Output/exp_slab2_06.14.18/exp_slab2_dep_06.14.18.grd -C/Users/ghayes/Desktop/CPTs/slab.cpt -O -K -J -R >> temp.ps
gmt grdcontour ../Output/exp_slab2_06.14.18/exp_slab2_dep_06.14.18.grd -J -R -O -K -C40 -W1p,20 >> temp.ps
awk -F, '{if(NR>1) print $2,$1,$3}' exp_04-18_input.csv | gmt psxy -J -R -O -K -Sc0.1 -W0.2p,100 -C/Users/ghayes/Desktop/CPTs/dep.cpt >> temp.ps
