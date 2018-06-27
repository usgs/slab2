#!/usr/bin/perl -w
####
if (@ARGV != 4){print "Must input 4 string arguments: \n";
                print "    1) slab (e.g. alu, sam, phi) \n";
                print "    2) model date (MM.DD.YY) \n";
                print "    3) model directory name and location \n";
                print "    4) input file name and location \n";
exit
}
($ID,$folder,$folderloc,$unsorted)=@ARGV;
####
`rm -R $folderloc/$ID\_slab2_figs_$folder`;
`mkdir $folderloc/$ID\_slab2_figs_$folder`;
#####
# SLAB2 GRID
###
$raw="$folderloc/$ID\_slab2_dep_$folder.grd";
$raws="$folderloc/$ID\_slab2_str_$folder.grd";
$rawd="$folderloc/$ID\_slab2_dip_$folder.grd";
$rawu="$folderloc/$ID\_slab2_unc_$folder.grd";
$rawt="$folderloc/$ID\_slab2_thk_$folder.grd";
$gmt="/Applications/GMT-5.4.3.app/Contents/Resources/bin/gmt";
$clipmask="$folderloc/$ID\_slab2_clp_$folder.csv";
$trenches="library/forplotting/trenches_usgs_2017_depths.csv";
#####
$minlon=999;
$maxlon=-999;
$uncmax=0;
$minlat=$minlon;$maxlat=$maxlon;

#####
# READS UNFILTERED INPUT FILE
#####
#lat,lon,depth,unc,etype,ID,S1,D1,R1,S2,D2,R2,src,distance
$unfiltered="$folderloc/$ID\_slab2_dat2_$folder.csv";
`sort -t, -nk3 $unsorted > $unfiltered`;
open(DATs,"$unfiltered");
   @dta=<DATs>;close(DATs);
$unfilt="$folderloc/$ID\_slab2_figs_$folder/tempUNFL.out";
open(UNFILT,">$unfilt");
$n=0;
$nn=0;
foreach(@dta){
    if($nn>0){
	($dlat[$n],$dlon[$n],$ddep[$n],$dunc[$n],$dtype[$n],$ID[$n],$S1[$n],$D1[$n],$R1[$n],$S2[$n],$D2[$n],$R2[$n],$dsrc[$n],$dist[$n])=split(',',$_);
    	chomp($dist[$n]);
	$ddep2[$n]=-$ddep[$n];
	if($dlon[$n]<0){$dlon[$n]=$dlon[$n]+360;}
	if($dlon[$n]<$minlon){$minlon=$dlon[$n];}
	if($dlon[$n]>$maxlon){$maxlon=$dlon[$n];}
        if($dlat[$n]<$minlat){$minlat=$dlat[$n];}
        if($dlat[$n]>$maxlat){$maxlat=$dlat[$n];}
	if($ddep[$n]){
	    print UNFILT "$dlon[$n] $dlat[$n] $ddep2[$n] $dtype[$n] $dunc[$n]\n";
	}
	$n++;
    }
    $nn++;
}
close(UNFILT);

####
# READS FILTERED INPUT FILE
$EQID = "EQ";
$ERID = "ER";
$ASID = "AS";
$TOID = "TO";
$AAID = "AA";
$BAID = "BA";
$RFID = "RF";
$CPID = "CP";
$data="$folderloc/$ID\_slab2_dat_$folder.csv";
open(DAT,"$data");
@din=<DAT>;close(DAT);
#lat,lon,depth,unc,etype,ID,S1,D1,R1,S2,D2,R2,src,distance
$oginput="$folderloc/$ID\_slab2_figs_$folder/tempOGIN.out";
open(OGINPUT,">$oginput");
$indataEQu="$folderloc/$ID\_slab2_figs_$folder/indataEQu.dat";
open(OUTEQu,">$indataEQu");
$indataERu="$folderloc/$ID\_slab2_figs_$folder/indataERu.dat";
open(OUTERu,">$indataERu");
$indataASu="$folderloc/$ID\_slab2_figs_$folder/indataASu.dat";
open(OUTASu,">$indataASu");
$indataTOu="$folderloc/$ID\_slab2_figs_$folder/indataTOu.dat";
open(OUTTOu,">$indataTOu");
$indataAAu="$folderloc/$ID\_slab2_figs_$folder/indataAAu.dat";
open(OUTAAu,">$indataAAu");
$indataBAu="$folderloc/$ID\_slab2_figs_$folder/indataBAu.dat";
open(OUTBAu,">$indataBAu");
$indataRFu="$folderloc/$ID\_slab2_figs_$folder/indataRFu.dat";
open(OUTRFu,">$indataRFu");
$indataCPu="$folderloc/$ID\_slab2_figs_$folder/indataCPu.dat";
open(OUTCPu,">$indataCPu");
$n=0;
$nn=0;
foreach(@din){
    if($nn>0){
	($dlat[$n],$dlon[$n],$ddep[$n],$dunc[$n],$dtype[$n],$dID[$n],$S1[$n],$D1[$n],$R1[$n],$S2[$n],$D2[$n],$R2[$n],$dsrc[$n],$dist[$n])=split(',',$_);
    	chomp($dist[$n]);
	$ddep2[$n]=-$ddep[$n];
	if($dlon[$n]<0){$dlon[$n]=$dlon[$n]+360;}
	if($ddep[$n]){
	    print OGINPUT "$dlon[$n] $dlat[$n] $ddep2[$n] $dtype[$n] $dunc[$n]\n";
        if($dtype[$n]eq$EQID){print OUTEQu "$dlon[$n] $dlat[$n] $ddep[$n] $dID[$n]\n";}
        if($dtype[$n]eq$ERID){print OUTERu "$dlon[$n] $dlat[$n] $ddep[$n] $dID[$n]\n";}
        if($dtype[$n]eq$ASID){print OUTASu "$dlon[$n] $dlat[$n] $ddep[$n] $dID[$n]\n";}
        if($dtype[$n]eq$TOID){print OUTTOu "$dlon[$n] $dlat[$n] $ddep[$n] $dID[$n]\n";}
        if($dtype[$n]eq$AAID){print OUTAAu "$dlon[$n] $dlat[$n] $ddep[$n] $dID[$n]\n";}
        if($dtype[$n]eq$BAID){print OUTBAu "$dlon[$n] $dlat[$n] $ddep[$n] $dID[$n]\n";}
        if($dtype[$n]eq$RFID){print OUTRFu "$dlon[$n] $dlat[$n] $ddep[$n] $dID[$n]\n";}
        if($dtype[$n]eq$CPID){print OUTCPu "$dlon[$n] $dlat[$n] $ddep[$n] $dID[$n]\n";}
	}
	$n++;
    }
    $nn++;
}
close(OGINPUT);

####
# READS NODE FILE (LON,LAT,DEP)
$nodes="$folderloc/$ID\_slab2_nod_$folder.csv";
open(DAT,"$data");
@din=<DAT>;close(DAT);
#lon,lat,dep_raw,err,shiftmag,shifterr,avstr,avdip,avrk,psdep,psstr,psdip,ID,prelon,prelat
$posshift="$folderloc/$ID\_slab2_figs_$folder/posShift.dat";
$preshift="$folderloc/$ID\_slab2_figs_$folder/preShift.dat";
open(OUTPOS,">$posshift");
open(OUTPRE,">$preshift");
open(DAT,"$nodes");
   @dta=<DAT>;close(DAT);
   $n=0;
   foreach(@dta) {
        ($nlon[$n],$nlat[$n],$ndep[$n],$nstdv[$n],$nsmag[$n],$nshiftstd[$n],$navstr[$n],$navdip[$n],$navrke[$n],$npsdep[$n],$nsstr[$n],$nsdip[$n],$nnID[$n],$npslon[$n],$npslat[$n],$bzlon[$n],$bzlat[$n],$ncentsurf[$n],$nthickness[$n],$nalen[$n],$nblen[$n],$nclen[$n],$nogstr[$n],$nogdip[$n])=split(',',$_);
        #($nlon[$n],$nlat[$n],$ndep[$n],$nerr[$n],$nshiftmag[$n],$nshifterr[$n],$navstr[$n],$navdip[$n],$navrk[$n],$npsdep[$n],$npsstr[$n],$npsdip[$n],$nID[$n],$nprelon[$n],$nprelat[$n],$bzlon[$n],$bzlat[$n])=split(',',$_);
        #chomp($nprelat[$n]);
        chomp($nogdip[$n]);
        if($nlon[$n]<0){$nlon[$n]=$nlon[$n]+360;}
        if($ndep[$n]!="nan"){
        print OUTPOS "$nlon[$n] $nlat[$n] -$ndep[$n]\n";
        #print OUTPRE "$nprelon[$n] $nprelat[$n] -$npsdep[$n]\n";
        print OUTPRE "$bzlon[$n] $bzlat[$n] -$npsdep[$n]\n";
        }
      $n++;
   }
close(OUTPOS);
close(OUTPRE);



$thistrench="$folderloc/$ID\_slab2_figs_$folder/outtrench1";
$plottrench="$folderloc/$ID\_slab2_figs_$folder/outtrench2";
open(PLOT,">$plottrench");
open(THIS,">$thistrench");

############################
# Get trench data into an array
open(DAT,"$trenches");
   @dta=<DAT>;close(DAT);
   $n=0;
   $nn=0;
   foreach(@dta) {
        ($lon1[$n],$lat1[$n],$str1[$n],$bound[$n],$slab[$n])=split(',',$_);
        chomp($slab[$n]);
        if($lon1[$n]<0){$lon1[$n]=$lon1[$n]+360;}
        if($slab[$n]eq$ID){print THIS "$lon1[$n] $lat1[$n] $str1[$n]\n";
                            $lon[$nn]=$lon1[$n];
                            $lat[$nn]=$lat1[$n];
                            $str[$nn]=$str1[$n];
                            $nn++;}
        print PLOT "$lon1[$n] $lat1[$n] $str1[$n]\n";
      $n++;
   }
close(THIS);
close(PLOT);

####
# GET BOUNDS OF GRID
###
@string=`$gmt grdinfo $raw -C`;
foreach(@string){
    ($a,$we,$ea,$so,$no,$z0,$z1,$dx,$dy,$nx,$ny)=split;
    $minlon2=$we;
    $maxlon2=$ea;
    $minlat2=$so;
    $maxlat2=$no;
    $minlon=$we;
    $maxlon=$ea;
    $minlat=$so;
    $maxlat=$no;
}
#####
# CREATES BOUNDS OF SLAB
#####
$minlon=&round_to_nths($minlon-1,1);
$maxlon=&round_to_nths($maxlon+1,1);
$minlat=&round_to_nths($minlat-1,1);
$maxlat=&round_to_nths($maxlat+1,1);
$bounds="-R$minlon/$maxlon/$minlat/$maxlat";
print "BOUNDS = $bounds\n";

@string=`$gmt grdinfo $raws -C`;
foreach(@string){
    ($a,$we,$ea,$so,$no,$s0,$s1,$dx,$dy,$nx,$ny)=split;
}
@string=`$gmt grdinfo $rawd -C`;
foreach(@string){
    ($a,$we,$ea,$so,$no,$d0,$d1,$dx,$dy,$nx,$ny)=split;
}
@string=`$gmt grdinfo $rawu -C`;
foreach(@string){
    ($a,$we,$ea,$so,$no,$u0,$u1,$dx,$dy,$nx,$ny)=split;
}
@string=`$gmt grdinfo $rawt -C`;
foreach(@string){
    ($a,$we,$ea,$so,$no,$t0,$t1,$dx,$dy,$nx,$ny)=split;
}

if($ID eq "ker" or $ID eq "izu" or $ID eq "sol" or $ID eq "man"){
    $tilted = "$folderloc/$ID\_slab2_sup_$folder.csv";
    $conts = "$folderloc/$ID\_slab2_cdep20_$folder.txt";
    $tiltsurfp="$folderloc/$ID\_slab2_figs_$folder/tiltedp.dat";
    $tiltsurfs="$folderloc/$ID\_slab2_figs_$folder/tilteds.dat";
    $tiltsurfd="$folderloc/$ID\_slab2_figs_$folder/tiltedd.dat";
    $tiltsurft="$folderloc/$ID\_slab2_figs_$folder/tiltedt.dat";
    $tiltsurfu="$folderloc/$ID\_slab2_figs_$folder/tiltedu.dat";
    open(TILTSURFP,">$tiltsurfp");
    open(TILTSURFS,">$tiltsurfs");
    open(TILTSURFD,">$tiltsurfd");
    open(TILTSURFU,">$tiltsurfu");
    open(TILTSURFT,">$tiltsurft");

    open(DAT,"$tilted");
       @dta=<DAT>;close(DAT);
       $n=0;
       #lat,lon,depth,unc,etype,ID,S1,D1,R1,S2,D2,R2,src,distance
       foreach(@dta) {
            ($ilon[$n],$ilat[$n],$idep[$n],$istr[$n],$idip[$n],$idz1[$n],$idz2[$n],$idz3[$n],$ithk[$n])=split(',',$_);
            chomp($ithk[$n]);
            if($ilon[$n]<0){$ilon[$n]=$ilon[$n]+360;}
            print TILTSURFP "$ilon[$n] $ilat[$n] $idep[$n]\n";
            print TILTSURFS "$ilon[$n] $ilat[$n] $istr[$n]\n";
            print TILTSURFD "$ilon[$n] $ilat[$n] $idip[$n]\n";
            print TILTSURFU "$ilon[$n] $ilat[$n] $idz1[$n]\n";
            print TILTSURFT "$ilon[$n] $ilat[$n] $ithk[$n]\n";
            if ($idep[$n]>$z0*-1){$z0=-1*$idep[$n];}
            if ($istr[$n]>$s1){$s1=$istr[$n];}
            if ($ithk[$n]>$t1){$t1=$ithk[$n];}
            if ($idz1[$n]>$u1){$u1=$idz1[$n];}
            if ($idip[$n]>$d1){$d1=$idip[$n];}
            $n++;
       }
    close(TILTSURFP);
    close(TILTSURFS);
    close(TILTSURFD);
    close(TILTSURFU);
    close(TILTSURFT);
    
}
$z02 = -1*$z0;
$z12 = -1*$z1;
print "$z0,$z1,$s0,$s1,$d0,$d1,$t0,$t1,$u0,$u1";
$depcpt="$folderloc/$ID\_slab2_figs_$folder/dep.cpt";
`$gmt makecpt -Cseis -I -D -Fr -T$z0/$z1/10 -Z > $depcpt`;
$depcpt2="$folderloc/$ID\_slab2_figs_$folder/dep2.cpt";
`$gmt makecpt -Cseis -D -Fr -T$z12/$z02/10 -Z > $depcpt2`;
$strcpt="$folderloc/$ID\_slab2_figs_$folder/str.cpt";
`$gmt makecpt -Cno_green -I -D -Fr -T$s0/$s1/10 -Z > $strcpt`;
$dipcpt="$folderloc/$ID\_slab2_figs_$folder/dip.cpt";
`$gmt makecpt -Cno_green -I -D -Fr -T$d0/$d1/10 -Z > $dipcpt`;
$unccpt="$folderloc/$ID\_slab2_figs_$folder/unc.cpt";
`$gmt makecpt -Cgray -I -D -Fr -T$u0/$u1/10 -Z > $unccpt`;
$thkcpt="$folderloc/$ID\_slab2_figs_$folder/thk.cpt";
`$gmt makecpt -Cno_green -I -D -Fr -T$t0/$t1/10 -Z > $thkcpt`;
$bath="library/forplotting/world_lr.grd";
$clipmask="$folderloc/$ID\_slab2_clp_$folder.csv";
$ghayscpt="library/forplotting/ghayes2.cpt";

$rlo1=$minlon;
$rlo2=$maxlon;
$rla1=$minlat;
$rla2=$maxlat+1;
$reg="-R$rlo1/$rlo2/$rla1/$rla2";

$londiff = $maxlon-$minlon;
$latdiff = $maxlat-$minlat;
$depdiff = $z0-$z1;

$B="-Bx10g5 -By10g5";
$CP="-B40";
$AP="-C10 -A20";

$CS="-B40";
$AS="-C5 -A10";

$CD="-B10";
$AD="-C5 -A10";

$CT="-B5";
$AT="-C2 -A10";

$CU="-B10";
$AU="-C5 -A10";
if ($londiff < 10){$B="-Bx10g2 -By10g2";}
if ($latdiff < 10){$B="-Bx10g2 -By10g2";}
if ($depdiff > 200){
    $CP="-B50";
    $AP="-C20 -A40";
    }
if ($depdiff < -200){
    $CP="-B50";
    $AP="-C20 -A40";
    }
if ($depdiff > 500){$C="-B100";}
if ($depdiff < -500){$C="-B100";}

$JM="-JM12.5i";
if($latdiff>(1.5*$londiff)){$JM="-JM8.2i";}
if($ID eq "sco"){$JM="-JM8.2i";}



####
# PLOT Final surface
####
`$gmt grdimage $bath -C$ghayscpt $reg $JM -P -K -BSnWe+t"Depth Grid $ID\_slab2_$folder" >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_depth.ps`;
`$gmt psclip $clipmask $JM $reg -P -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_depth.ps`;
`$gmt grdimage $raw -C$depcpt -R -J -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_depth.ps`;
`$gmt grdcontour $raw $AP -J -R -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_depth.ps`;
`$gmt psclip $clipmask -J -R -P -C -O -K -V >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_depth.ps`;
if($ID eq "ker" or $ID eq "izu" or $ID eq "sol" or $ID eq "man"){
`$gmt psxy $tiltsurfp -J -R -O -K -Sc0.08 -C$depcpt2 >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_depth.ps`;
`$gmt psxy -J -R -O -K -Wblack $conts -P >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_depth.ps`;
}
`$gmt pscoast -J -R -O -K -Df -W0.75p,100 $B >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_depth.ps`;
`$gmt psxy -J -R -O -K -W2p $clipmask -P >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_depth.ps`;
`$gmt psscale -D0.5/3.0/5/0.5 -C$depcpt $CP -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_depth.ps`;

`$gmt grdimage $bath -C$ghayscpt $reg $JM -P -K -BSnWe+t"Strike Grid $ID\_slab2_$folder" >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_strike.ps`;
`$gmt psclip $clipmask $JM $reg -P -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_strike.ps`;
`$gmt grdimage $raws -C$strcpt -R -J -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_strike.ps`;
`$gmt grdcontour $raws $AS -J -R -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_strike.ps`;
`$gmt psclip $clipmask -J -R -P -C -O -K -V >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_strike.ps`;
if($ID eq "ker" or $ID eq "izu" or $ID eq "sol" or $ID eq "man"){
`$gmt psxy $tiltsurfs -J -R -O -K -Sc0.08 -C$strcpt >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_strike.ps`;
`$gmt psxy -J -R -O -K -Sc0.01,black $folderloc/$ID\_slab2_cstr10_$folder.txt -P >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_strike.ps`;
}
`$gmt pscoast -J -R -O -K -Df -W0.75p,100 $B >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_strike.ps`;
`$gmt psxy -J -R -O -K -W2p $clipmask -P >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_strike.ps`;
`$gmt psscale -D0.5/3.0/5/0.5 -C$strcpt $CS -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_strike.ps`;

`$gmt grdimage $bath -C$ghayscpt $reg $JM -P -K -BSnWe+t"Dip Grid $ID\_slab2_$folder" >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_dip.ps`;
`$gmt psclip $clipmask $JM $reg -P -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_dip.ps`;
`$gmt grdimage $rawd -C$dipcpt -R -J -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_dip.ps`;
`$gmt grdcontour $rawd $AD -J -R -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_dip.ps`;
`$gmt psclip $clipmask -J -R -P -C -O -K -V >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_dip.ps`;
if($ID eq "ker" or $ID eq "izu" or $ID eq "sol" or $ID eq "man"){
`$gmt psxy $tiltsurfd -J -R -O -K -Sc0.08 -C$dipcpt >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_dip.ps`;
`$gmt psxy -J -R -O -K -Sc0.01,black $folderloc/$ID\_slab2_cdip5_$folder.txt -P >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_dip.ps`;
}
`$gmt pscoast -J -R -O -K -Df -W0.75p,100 $B >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_dip.ps`;
`$gmt psxy -J -R -O -K -W2p $clipmask -P >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_dip.ps`;
`$gmt psscale -D0.5/3.0/5/0.5 -C$dipcpt $CD -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_dip.ps`;


`$gmt grdimage $bath -C$ghayscpt $reg $JM -P -K -BSnWe+t"Unc Grid$ID\_slab2_$folder" >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_unc.ps`;
`$gmt psclip $clipmask $JM $reg -P -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_unc.ps`;
`$gmt grdimage $rawu -C$unccpt -R -J -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_unc.ps`;
`$gmt grdcontour $rawu $AU -J -R -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_unc.ps`;
`$gmt psclip $clipmask -J -R -P -C -O -K -V >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_unc.ps`;
if($ID eq "ker" or $ID eq "izu" or $ID eq "sol" or $ID eq "man"){
`$gmt psxy $tiltsurfu -J -R -O -K -Sc0.08 -C$unccpt >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_unc.ps`;
}
`$gmt pscoast -J -R -O -K -Df -W0.75p,100 $B >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_unc.ps`;
`$gmt psxy -J -R -O -K -W2p $clipmask -P >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_unc.ps`;
`$gmt psscale -D0.5/3.0/5/0.5 -C$unccpt $CU -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_unc.ps`;

`$gmt grdimage $bath -C$ghayscpt $reg $JM -P -K -BSnWe+t"Thickness Grid$ID\_slab2_$folder" >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_thk.ps`;
`$gmt psclip $clipmask $JM $reg -P -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_thk.ps`;
`$gmt grdimage $rawt -C$thkcpt -R -J -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_thk.ps`;
`$gmt grdcontour $rawt $AT -J -R -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_thk.ps`;
`$gmt psclip $clipmask -J -R -P -C -O -K -V >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_thk.ps`;
if($ID eq "ker" or $ID eq "izu" or $ID eq "sol" or $ID eq "man"){
`$gmt psxy $tiltsurft -J -R -O -K -Sc0.08 -C$thkcpt >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_thk.ps`;
}
`$gmt pscoast -J -R -O -K -Df -W0.75p,100 $B >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_thk.ps`;
`$gmt psxy -J -R -O -K -W2p $clipmask -P >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_thk.ps`;
`$gmt psscale -D0.5/3.0/5/0.5 -C$thkcpt $CT -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_thk.ps`;

`$gmt grdimage $bath -C$ghayscpt $reg $JM -P -K -BSnWe+t"Shifted Nodes $ID\_slab2_$folder" >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_nodes.ps`;
`$gmt psxy $posshift -J -R -O -K -Sc0.08 -W0.05p -C$depcpt >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_nodes.ps`;
`$gmt psxy -J -R -O -K -W2p $clipmask -P >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_nodes.ps`;
`$gmt pscoast -J -R -O -K -Df -W0.75p,100 $B >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_nodes.ps`;
`$gmt psscale -D0.5/3.0/5/0.5 -C$depcpt $CP -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_nodes.ps`;

`$gmt grdimage $bath -C$ghayscpt $reg $JM -P -K -BSnWe+t"Filtered Dataset $ID\_slab2_$folder" >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_filtered.ps`;
`$gmt psxy $oginput -J -R -O -K -Sc0.08 -W0.05p -C$depcpt >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_filtered.ps`;
`$gmt psxy -J -R -O -K -W2p $clipmask -P >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_filtered.ps`;
`$gmt pscoast -J -R -O -K -Df -W0.75p,100 $B >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_filtered.ps`;
`$gmt psscale -D0.5/3.0/5/0.5 -C$depcpt $CP -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_filtered.ps`;

`$gmt grdimage $bath -C$ghayscpt $reg $JM -P -K -BSnWe+t"Unfiltered Dataset $ID\_slab2_$folder" >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_unfiltered.ps`;
`$gmt psxy $unfilt -J -R -O -K -Sc0.08 -W0.05p -C$depcpt >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_unfiltered.ps`;
`$gmt psxy -J -R -O -K -W2p $clipmask -P >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_unfiltered.ps`;
`$gmt pscoast -J -R -O -K -Df -W0.75p,100 $B >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_unfiltered.ps`;
`$gmt psscale -D0.5/3.0/5/0.5 -C$depcpt $CP -O -K >> $folderloc/$ID\_slab2_figs_$folder/$ID\_slab2_$folder\_unfiltered.ps`;


#unlink "$folderloc/$ID\_slab2_figs_$folder/nodes.out";
unlink "$folderloc/$ID\_slab2_figs_$folder/$folderloc/$ID\_slab2_figs_$folder/dep.cpt";
unlink "$folderloc/$ID\_slab2_figs_$folder/outtrench1";
unlink "$folderloc/$ID\_slab2_figs_$folder/outtrench2";
unlink "$folderloc/$ID\_slab2_figs_$folder/posShift.dat";
unlink "$folderloc/$ID\_slab2_figs_$folder/preShift.dat";
unlink "$folderloc/$ID\_slab2_figs_$folder/tempOGIN.out";
unlink "$folderloc/$ID\_slab2_figs_$folder/tempUNFL.out";
unlink "$folderloc/$ID\_slab2_dat2_$folder.csv";
unlink "$folderloc/$ID\_slab2_dat2_$folder.csv";
unlink "$folderloc/$ID\_slab2_figs_$folder/dep.cpt";
unlink "$folderloc/$ID\_slab2_figs_$folder/dep2.cpt";
unlink "$folderloc/$ID\_slab2_figs_$folder/dip.cpt";
unlink "$folderloc/$ID\_slab2_figs_$folder/str.cpt";
unlink "$folderloc/$ID\_slab2_figs_$folder/unc.cpt";
unlink "$folderloc/$ID\_slab2_figs_$folder/thk.cpt";
unlink "$folderloc/$ID\_slab2_figs_$folder/indataAAu.dat";
unlink "$folderloc/$ID\_slab2_figs_$folder/indataASu.dat";
unlink "$folderloc/$ID\_slab2_figs_$folder/indataBAu.dat";
unlink "$folderloc/$ID\_slab2_figs_$folder/indataCPu.dat";
unlink "$folderloc/$ID\_slab2_figs_$folder/indataEQu.dat";
unlink "$folderloc/$ID\_slab2_figs_$folder/indataERu.dat";
unlink "$folderloc/$ID\_slab2_figs_$folder/indataRFu.dat";
unlink "$folderloc/$ID\_slab2_figs_$folder/indataTOu.dat";
unlink "$folderloc/$ID\_slab2_figs_$folder/titledd.dat";
unlink "$folderloc/$ID\_slab2_figs_$folder/titledp.dat";
unlink "$folderloc/$ID\_slab2_figs_$folder/titleds.dat";
unlink "$folderloc/$ID\_slab2_figs_$folder/titledt.dat";
unlink "$folderloc/$ID\_slab2_figs_$folder/titledu.dat";

#####################
#####################
sub round_to_nths {
    my ($num, $n) = @_;
    (int $num*$n)/$n
}
######################







