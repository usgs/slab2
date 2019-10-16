gmt grdimage ../../misc/bathymetry/gebco5x5.grd -Cghayes2.cpt -R123/132/-4/5 -JM15 -K > temp2019Map.ps                      
gmt pscoast -J -R -O -K -Df -W0.6p,50 -BSnWe -Bx4g1 -By4g1 >> temp2019Map.ps                                                  
gmt grdimage exp_guide2.grd -Cdep.cpt -O -K -J -R >> temp2019Map.ps                                                  
gmt grdcontour -J -R -O -K exp_guide2.grd -C40 -W1p,200 >> temp2019Map.ps                                                                        
gmt psxy exp_poly.txt -J -R -O -K -W1.5p,50 >> temp2019Map.ps                                                                                    
gmt grdimage ../Output/exp_slab2_10.03.19/exp_slab2_dep_10.03.19.grd -Cdep.cpt -O -K -J -R >> temp2019Map.ps 
gmt grdcontour ../Output/exp_slab2_10.03.19/exp_slab2_dep_10.03.19.grd -J -R -O -K -C40 -W1p,20 >> temp2019Map.ps 
awk -F, '{if(NR>1) print $2,$1,$3}' exp_10-19_input.csv | gmt psxy -J -R -O -K -Sc0.1 -W0.2p,100 -Cdep.cpt >> temp2019Map.ps 
