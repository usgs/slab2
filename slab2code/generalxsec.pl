#!/usr/bin/perl -w
############################

use Math::Trig;
use File::Copy qw(move);
use File::Copy qw(copy);

############################

####
if (@ARGV < 6){print "Must input 6 or 14 string arguments: \n";
                print "    1) slab (e.g. alu, sam, phi) \n";
                print "    2) model date (MM.DD.YY) \n";
                print "    3) model directory name and location \n";
                print "    4) input file name and location \n";
                print "    5) Green grid to compare \n";
                print "    6) spacing (km) between consecutive sections \n";
                print "    OPTIONALLY include all of these arguments: \n";
                print "      7) trench file to use other than default. Enter na for default \n";
                print "               default is trenches_usgs_2017_depths.csv \n";
                print "      8) constant azimuth direction. Enter na for trench azimuth + 90º\n";
                print "      9) cyan grid to compare. Enter na for no extra grid \n";
                print "      10) length of cross section. Enter na for default \n";
                print "               default is maximum data point depth (km) \n";
                print "      11) back length of cross section. Enter na for default \n";
                print "               default is -50 (km) \n";
                print "      12) depth extent of cross section. Enter na for default \n";
                print "               default is maximum data point depth (km) \n";
                print "      13) additional supplementary dataset to compare (green) \n";
                print "                Enter na for no extra model \n";
                print "      14) additional supplementary dataset to compare (cyan) \n";
                print "                Enter na for no extra model \n";
exit
}
if (@ARGV > 14){print "Must input 6 or 14 string arguments: \n";
                print "    1) slab (e.g. alu, sam, phi) \n";
                print "    2) model date (MM.DD.YY) \n";
                print "    3) model directory name and location \n";
                print "    4) input file name and location \n";
                print "    5) Green grid to compare \n";
                print "    6) spacing (km) between consecutive sections \n";
                print "    OPTIONALLY include all of these arguments: \n";
                print "      7) trench file to use other than default. Enter na for default \n";
                print "               default is trenches_usgs_2017_depths.csv \n";
                print "      8) constant azimuth direction. Enter na for trench azimuth + 90º\n";
                print "      9) cyan grid to compare. Enter na for no extra grid \n";
                print "      10) length of cross section. Enter na for default \n";
                print "               default is maximum data point depth (km) \n";
                print "      11) back length of cross section. Enter na for default \n";
                print "               default is -50 (km) \n";
                print "      12) depth extent of cross section. Enter na for default \n";
                print "               default is maximum data point depth (km) \n";
                print "      13) additional supplementary dataset to compare (green) \n";
                print "                Enter na for no extra model \n";
                print "      14) additional supplementary dataset to compare (cyan) \n";
                print "                Enter na for no extra model \n";
exit
}
if (@ARGV < 14 & @ARGV > 6){print "Must input 6 or 14 string arguments: \n";
                print "    1) slab (e.g. alu, sam, phi) \n";
                print "    2) model date (MM.DD.YY) \n";
                print "    3) model directory name and location \n";
                print "    4) input file name and location \n";
                print "    5) Green grid to compare \n";
                print "    6) spacing (km) between consecutive sections \n";
                print "    OPTIONALLY include all of these arguments: \n";
                print "      7) trench file to use other than default. Enter na for default \n";
                print "               default is trenches_usgs_2017_depths.csv \n";
                print "      8) constant azimuth direction. Enter na for trench azimuth + 90º\n";
                print "      9) cyan grid to compare. Enter na for no extra grid \n";
                print "      10) length of cross section. Enter na for default \n";
                print "               default is maximum data point depth (km) \n";
                print "      11) back length of cross section. Enter na for default \n";
                print "               default is -50 (km) \n";
                print "      12) depth extent of cross section. Enter na for default \n";
                print "               default is maximum data point depth (km) \n";
                print "      13) additional supplementary dataset to compare (green) \n";
                print "                Enter na for no extra model \n";
                print "      14) additional supplementary dataset to compare (cyan) \n";
                print "                Enter na for no extra model \n";
exit
}


$reg_trench="library/forplotting/trenches_usgs_2017_depths.csv";
$saveflag=2;
$xclen="na";
$caz="na";
$guide="na";
$xclen2="na";
$maxdepthplot="na";
$tilted2="na";
$tilted3="na";
if (@ARGV == 6){($ID,$folder,$folderloc,$unfiltered,$tgs,$sp)=@ARGV;}
if (@ARGV == 14){($ID,$folder,$folderloc,$unfiltered,$tgs,$sp,$reg_trench,$caz,$guide,$xclen,$xclen2,$maxdepthplot,$tilted2,$tilted3)=@ARGV;}
if($reg_trench eq "na"){$reg_trench="library/forplotting/trenches_usgs_2017_depths.csv";}

print "$ID,$folder,$folderloc,$unfiltered,$reg_trench,$caz,$guide,$xclen";


`rm -r $folderloc/$ID\_slab2_xsecs_$folder`;
`mkdir $folderloc/$ID\_slab2_xsecs_$folder`;
`mkdir $folderloc/$ID\_slab2_xsecs_$folder/indata`;

####
# gmt must be on the PATH
$gmt="gmt";
## Set GMT output bounds
`$gmt gmtset PS_MEDIA=a2`;
$raw="$folderloc/$ID\_slab2_dep_$folder.grd";
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

## Cut down bathymetry data
$GEBCO="library/forplotting/world_lr.grd";
$R="-R$minlon/$maxlon/$minlat/$maxlat";

## Trench files
$outtrench="$folderloc/$ID\_slab2_xsecs_$folder/$ID\_finaltrench.dat";
$outtrench1="$folderloc/$ID\_slab2_xsecs_$folder/$ID\_trench.dat";
$out1="$outtrench1";
open(OUT1,">$out1");

## Order trench data from min lat to max lat
`sort -nk2 $reg_trench`;
############################
# Get trench data into an array
$tID="$ID";
open(DAT,"$reg_trench");
   @dta=<DAT>;close(DAT);
   $n=0;
   $nn=0;
   foreach(@dta) {
        ($lon1[$n],$lat1[$n],$str1[$n],$bound[$n],$slab[$n])=split(',',$_);
        chomp($slab[$n]);
        if($lon1[$n]<0){$lon1[$n]=$lon1[$n]+360;}
        if($slab[$n]eq$tID){print OUT1 "$lon1[$n] $lat1[$n] $str1[$n]\n";
                            $lon[$nn]=$lon1[$n];
                            $lat[$nn]=$lat1[$n];
                            $str[$nn]=$str1[$n];
                            $nn++;}
      $n++;
   }
close(OUT1);
`sort -nk1 $outtrench1`;
# Loop through the trench data we have and project to get more data points
$n=1;
$nmax=scalar @lon;
$count=0;
while ($n < $nmax) {
   if (($lon[$n-1]!=$lon[$n]) && ($lat[$n-1]!=$lat[$n])) {
      `$gmt project -C$lon[$n-1]/$lat[$n-1] -E$lon[$n]/$lat[$n] -G1 -Q > $folderloc/$ID\_slab2_xsecs_$folder/temp.dat`;
      $l=0;
      open(TEM,"$folderloc/$ID\_slab2_xsecs_$folder/temp.dat");
      @tdta=<TEM>;close(TEM);
      foreach(@tdta) {
         ($plon[$l],$plat[$l],$pdis[$l])=split;
         chomp($plon[$l],$plat[$l],$pdis[$l]);
         $l++;
      }
      $m=0;
      $refdis=0;
      $mmax=scalar @plon;
      foreach(@plon) {
        if(($str[$n-1] > 270) && ($str[$n] < 90)){
            $str[$n] = $str[$n]+360}
        if(($str[$n] > 270) && ($str[$n-1] < 90)){
            $str[$n-1] = $str[$n-1]+360}
         $slp=($str[$n-1]-$str[$n])/(-1*$pdis[$mmax-1]);
         $pstr[$m]=($slp*$pdis[$m])+$str[$n-1];
         $r=$m+$count;
         $trlon[$r]=$plon[$m];
         $trlat[$r]=$plat[$m];
         $trstr[$r]=$pstr[$m];
         $trdis[$r]=$pdis[$m];
         $m++;
      }
      $len=scalar @plon;
      $count=$count+$len;
      splice @plon;
      splice @plat;
      splice @pdist;
      splice @pstr;
      splice @totdis;
      unlink '$folderloc/$ID\_slab2_xsecs_$folder/temp.dat';
   }
   $n++;
}

# Calculate the distance of each trench location from the starting location
$m=0;
$totdis[0]=0;
$refdis=0;
foreach(@trlon) {
   if ($m>0 && $trdis[$m]==0) {
      $refdis=$refdis+$trdis[$m-1];
      $totdis[$m]=$totdis[$m-1];
   } else {
      $totdis[$m]=$refdis+$trdis[$m];
   }
   if ($trlon[$m]<0) {
      $trlon[$m]=$trlon[$m]+360;
   }
   $m++;
}

## Convert strike values to azimuth for projections
$n=0;
foreach(@trstr) {
   if ($trstr[$n]<90){
      $trstr[$n]=$trstr[$n]+360;
    }
    $traz[$n]=$trstr[$n]-90;
   $n++;
}

## Get rid of duplicate trench values
$n=0;
$l=0;
foreach(@trlon) {
   if (($trlon[$n]!=$trlon[$n-1]) && ($trlon[$n]!=$trlon[$n-1])) {
      $ftrlon[$l]=$trlon[$n];
      $ftrlat[$l]=$trlat[$n];
      $ftrstr[$l]=$trstr[$n]+180;
      $ftraz[$l]=$traz[$n]+180;
      $ftrdis[$l]=$totdis[$n];
      $l++;
   }
   $n++;
}

## Start at top or bottom of trench and place 50 km spaced trench data into arrays
$len=scalar @ftrlon;
$nsecs=int($ftrdis[$len-1]/$sp);
$count=0;
while ($count < $nsecs+1) {
   $loc=$count*$sp;
   $dis=10000;
   $ind=-1;
   $l=0;
   foreach(@ftrdis) {
      $temp=abs($ftrdis[$l]-$loc);
      if ($temp<$dis) {
         $ind=$l;
         $dis=$temp;
      }
      $l++;
   }
   $sftrlon[$count]=$ftrlon[$ind];
   $sftrlat[$count]=$ftrlat[$ind];
   $sftrstr[$count]=$ftrstr[$ind];
   $sftraz[$count]=$ftraz[$ind];
   $sftrdis[$count]=$ftrdis[$ind];
   $count++;
}

## Print final trench data from arrays into a file (spaced at 50 km)
unlink($outtrench);
open(OUTT,">$outtrench") or die "Could not open file '$outtrench' $!";

$n=0;
foreach(@sftrlon) {
   #if ($sftraz[$n]>190) {
    #  $sftraz[$n]=$sftraz[$n]-180;
   #}
   print OUTT "$sftrlon[$n] $sftrlat[$n] $sftraz[$n] $sftrdis[$n] $sftrstr[$n]\n";
   $n++;
}


## Data files:
$data="$folderloc/$ID\_slab2_dat_$folder.csv";
$nodes="$folderloc/$ID\_slab2_nod_$folder.csv";
$filler="$folderloc/$ID\_slab2_fil_$folder.csv";
$results="$folderloc/$ID\_slab2_res_$folder.csv";
$tilted = "$folderloc/$ID\_slab2_sup_$folder.csv";
$clipmask="$folderloc/$ID\_slab2_clp_$folder.csv";
$depcpt="$folderloc/$ID\_slab2_xsecs_$folder/dep.cpt";
`$gmt makecpt -Cseis -I -D -Fr -T$z0/$z1/10 -Z > $depcpt`;


$cont340="$folderloc/$ID\_slab2_340_$folder.csv";
$cont360="$folderloc/$ID\_slab2_360_$folder.csv";
$cont380="$folderloc/$ID\_slab2_380_$folder.csv";
$cont400="$folderloc/$ID\_slab2_400_$folder.csv";
$cont420="$folderloc/$ID\_slab2_420_$folder.csv";
$cont440="$folderloc/$ID\_slab2_440_$folder.csv";
$cont460="$folderloc/$ID\_slab2_460_$folder.csv";
$cont480="$folderloc/$ID\_slab2_480_$folder.csv";
$cont500="$folderloc/$ID\_slab2_500_$folder.csv";
$cont520="$folderloc/$ID\_slab2_520_$folder.csv";
$cont540="$folderloc/$ID\_slab2_540_$folder.csv";
$cont560="$folderloc/$ID\_slab2_560_$folder.csv";
$cont580="$folderloc/$ID\_slab2_580_$folder.csv";
$cont600="$folderloc/$ID\_slab2_600_$folder.csv";
$supp_izu="$folderloc/$ID\_slab2_sud_$folder.csv";

# Extracting input data and writing to individual files for plotting
$EQID = "EQ";
$ERID = "ER";
$ASID = "AS";
$TOID = "TO";
$AAID = "AA";
$BAID = "BA";
$RFID = "RF";
$CPID = "CP";
$indataEQ="$folderloc/$ID\_slab2_xsecs_$folder/indataEQ.dat";
open(OUTEQ,">$indataEQ");
$indataER="$folderloc/$ID\_slab2_xsecs_$folder/indataER.dat";
open(OUTER,">$indataER");
$indataAS="$folderloc/$ID\_slab2_xsecs_$folder/indataAS.dat";
open(OUTAS,">$indataAS");
$indataTO="$folderloc/$ID\_slab2_xsecs_$folder/indataTO.dat";
open(OUTTO,">$indataTO");
$indataAA="$folderloc/$ID\_slab2_xsecs_$folder/indataAA.dat";
open(OUTAA,">$indataAA");
$indataBA="$folderloc/$ID\_slab2_xsecs_$folder/indataBA.dat";
open(OUTBA,">$indataBA");
$indataRF="$folderloc/$ID\_slab2_xsecs_$folder/indataRF.dat";
open(OUTRF,">$indataRF");
$indataCP="$folderloc/$ID\_slab2_xsecs_$folder/indataCP.dat";
open(OUTCP,">$indataCP");
open(DAT,"$data");
   @dta=<DAT>;close(DAT);
   $n=0;
   foreach(@dta) {
        ($ilat[$n],$ilon[$n],$idep[$n],$iunc[$n],$itype[$n],$iid[$n],$imag[$n],$itime[$n],$iS1[$n],$iD1[$n],$iR1[$n],$iS2[$n],$iD2[$n],$iR2[$n],$isrc[$n])=split(',',$_);
        chomp($isrc[$n]);
        if($ilon[$n]<0){$ilon[$n]=$ilon[$n]+360;}
        if($itype[$n]eq$EQID){print OUTEQ "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        if($itype[$n]eq$ERID){print OUTER "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        if($itype[$n]eq$ASID){print OUTAS "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        if($itype[$n]eq$TOID){print OUTTO "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        if($itype[$n]eq$AAID){print OUTAA "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        if($itype[$n]eq$BAID){print OUTBA "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        if($itype[$n]eq$RFID){print OUTRF "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        if($itype[$n]eq$CPID){print OUTCP "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        $n++;
   }
close(OUTEQ);
close(OUTER);
close(OUTAS);
close(OUTTO);
close(OUTAA);
close(OUTBA);
close(OUTRF);
close(OUTCP);

#______________________ FOR SHIFTED DATA ______________________________________

$tiltsurf="$folderloc/$ID\_slab2_xsecs_$folder/tilted.dat";
open(TILTSURF,">$tiltsurf");
open(DAT,"$tilted");
   @dta=<DAT>;close(DAT);
   $n=0;
   #lat,lon,depth,unc,etype,ID,S1,D1,R1,S2,D2,R2,src,distance
   foreach(@dta) {
        ($ilon[$n],$ilat[$n],$idep[$n])=split(',',$_);
        chomp($idep[$n]);
        if($ilon[$n]<0){$ilon[$n]=$ilon[$n]+360;}
        print TILTSURF "$ilon[$n] $ilat[$n] $idep[$n] $idep[$n]\n";
        $n++;
   }
close(TILTSURF);

if ($tilted2 ne "na"){
    $tiltsurf2="$folderloc/$ID\_slab2_xsecs_$folder/tilted2.dat";
    open(TILTSURF2,">$tiltsurf2");
    open(DAT,"$tilted2");
       @dta=<DAT>;close(DAT);
       $n=0;
       #lat,lon,depth,unc,etype,ID,S1,D1,R1,S2,D2,R2,src,distance
       foreach(@dta) {
            ($ilon[$n],$ilat[$n],$idep[$n])=split(',',$_);
            chomp($idep[$n]);
            if($ilon[$n]<0){$ilon[$n]=$ilon[$n]+360;}
            print TILTSURF2 "$ilon[$n] $ilat[$n] $idep[$n] $idep[$n]\n";
            $n++;
       }
    close(TILTSURF2);
}

if ($tilted3 ne "na"){
    $tiltsurf3="$folderloc/$ID\_slab2_xsecs_$folder/tilted3.dat";
    open(TILTSURF3,">$tiltsurf3");
    open(DAT,"$tilted3");
       @dta=<DAT>;close(DAT);
       $n=0;
       #lat,lon,depth,unc,etype,ID,S1,D1,R1,S2,D2,R2,src,distance
       foreach(@dta) {
            ($ilon[$n],$ilat[$n],$idep[$n])=split(',',$_);
            chomp($idep[$n]);
            if($ilon[$n]<0){$ilon[$n]=$ilon[$n]+360;}
            print TILTSURF3 "$ilon[$n] $ilat[$n] $idep[$n] $idep[$n]\n";
            $n++;
       }
    close(TILTSURF3);
}

#______________________ FOR SHIFTED DATA ______________________________________


#______________________ FOR unfiltered DATA ______________________________________

$indataEQu="$folderloc/$ID\_slab2_xsecs_$folder/indataEQu.dat";
open(OUTEQu,">$indataEQu");
$indataERu="$folderloc/$ID\_slab2_xsecs_$folder/indataERu.dat";
open(OUTERu,">$indataERu");
$indataASu="$folderloc/$ID\_slab2_xsecs_$folder/indataASu.dat";
open(OUTASu,">$indataASu");
$indataTOu="$folderloc/$ID\_slab2_xsecs_$folder/indataTOu.dat";
open(OUTTOu,">$indataTOu");
$indataAAu="$folderloc/$ID\_slab2_xsecs_$folder/indataAAu.dat";
open(OUTAAu,">$indataAAu");
$indataBAu="$folderloc/$ID\_slab2_xsecs_$folder/indataBAu.dat";
open(OUTBAu,">$indataBAu");
$indataRFu="$folderloc/$ID\_slab2_xsecs_$folder/indataRFu.dat";
open(OUTRFu,">$indataRFu");
$indataCPu="$folderloc/$ID\_slab2_xsecs_$folder/indataCPu.dat";
open(OUTCPu,">$indataCPu");
open(DATu,"$unfiltered");
   @dtu=<DATu>;close(DATu);
   $n=0;
   foreach(@dtu) {
        ($ilat[$n],$ilon[$n],$idep[$n],$iunc[$n],$itype[$n],$imag[$n],$itime[$n],$iPaz[$n],$iPpl[$n],$iTaz[$n],$iTpl[$n],$iS1[$n],$iD1[$n],$iR1[$n],$iS2[$n],$iD2[$n],$iR2[$n],$imlon[$n],$imlat[$n],$imdep[$n],$iid[$n],$isrc[$n])=split(',',$_);
        chomp($isrc[$n]);
        if($ilon[$n]<0){$ilon[$n]=$ilon[$n]+360;}
        if($itype[$n]eq$EQID){print OUTEQu "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        if($itype[$n]eq$ERID){print OUTERu "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        if($itype[$n]eq$ASID){print OUTASu "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        if($itype[$n]eq$TOID){print OUTTOu "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        if($itype[$n]eq$AAID){print OUTAAu "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        if($itype[$n]eq$BAID){print OUTBAu "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        if($itype[$n]eq$RFID){print OUTRFu "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        if($itype[$n]eq$CPID){print OUTCPu "$ilon[$n] $ilat[$n] $idep[$n] $iid[$n]\n";}
        if(($idep[$n]*-1)<$z0){$z0=$idep[$n]*-1;}
        $n++;
   }
close(OUTEQu);
close(OUTERu);
close(OUTASu);
close(OUTTOu);
close(OUTAAu);
close(OUTBAu);
close(OUTRFu);
close(OUTCPu);

#______________________ FOR unfiltered DATA ______________________________________


# Extracting input nodes (pre and post shifted) and writing to individual files for plotting
$preshift="$folderloc/$ID\_slab2_xsecs_$folder/preShift.dat";
open(OUTPRE,">$preshift");
$posshift="$folderloc/$ID\_slab2_xsecs_$folder/posShift.dat";
open(OUTPOS,">$posshift");
$rezultz1="$folderloc/$ID\_slab2_xsecs_$folder/rezultz1.dat";
open(RESULTS1,">$rezultz1");
$rezultz2="$folderloc/$ID\_slab2_xsecs_$folder/rezultz2.dat";
open(RESULTS2,">$rezultz2");


open(DAT,"$nodes");
   @dta=<DAT>;close(DAT);
   $n=0;
   foreach(@dta) {
        ($nlon[$n],$nlat[$n],$ndep[$n],$nstdv[$n],$nsmag[$n],$nshiftstd[$n],$navstr[$n],$navdip[$n],$navrke[$n],$npsdep[$n],$nsstr[$n],$nsdip[$n],$nnID[$n],$npslon[$n],$npslat[$n],$bzlon[$n],$bzlat[$n],$ncentsurf[$n],$nthickness[$n],$nalen[$n],$nblen[$n],$nclen[$n],$nogstr[$n],$nogdip[$n],$hstdv[$n],$vstdv[$n])=split(',',$_);
        chomp($nogdip[$n]);
        if($nlon[$n]<0){$nlon[$n]=$nlon[$n]+360;}
        if($ndep[$n]!="nan"){
        print OUTPOS "$nlon[$n] $nlat[$n] $ndep[$n] $nnID[$n]\n";
        print OUTPRE "$bzlon[$n] $bzlat[$n] $npsdep[$n] $nnID[$n]\n";
        }
      $n++;
   }
close(OUTPOS);
close(OUTPRE);

# Extracting input nodes (pre and post shifted fillers) and writing to individual files for plotting
$preshift2="$folderloc/$ID\_slab2_xsecs_$folder/preShift2.dat";
open(OUTPRE2,">$preshift2");
$posshift2="$folderloc/$ID\_slab2_xsecs_$folder/posShift2.dat";
open(OUTPOS2,">$posshift2");

open(DAT,"$filler");
   @dta=<DAT>;close(DAT);
   $n=0;
   foreach(@dta) {
        ($nlon[$n],$nlat[$n],$ndep[$n],$nstdv[$n],$nsmag[$n],$nshiftstd[$n],$navstr[$n],$navdip[$n],$navrke[$n],$npsdep[$n],$nsstr[$n],$nsdip[$n],$nnID[$n],$npslon[$n],$npslat[$n],$bzlon[$n],$bzlat[$n],$ncentsurf[$n],$nthickness[$n],$nalen[$n],$nblen[$n],$nclen[$n],$nogstr[$n],$nogdip[$n],$hstdv[$n],$vstdv[$n])=split(',',$_);
        chomp($nogdip[$n]);
        if($nlon[$n]<0){$nlon[$n]=$nlon[$n]+360;}
        if($ndep[$n]!="nan"){
        print OUTPOS2 "$nlon[$n] $nlat[$n] $ndep[$n] $nnID[$n]\n";
        print OUTPRE2 "$bzlon[$n] $bzlat[$n] $npsdep[$n] $nnID[$n]\n";
        }
      $n++;
   }
close(OUTPOS2);
close(OUTPRE2);



############################

## Check basemaps to ensure trench locations are corret

# Basemap parameters
$J1="-JM15";
$R1="-R$minlon/$maxlon/$minlat/$maxlat";
$B1="-Bx5g -By5g -BSWen";
$C1="-Crelief";
$W1="-W1.5,red";

############################

## Create cross-sections

# Cross-section Variables
# Cross-section Variables
if($xclen eq "na"){$l=$z0*-1;}
else{$l=$xclen;}
if($xclen2 eq "na"){$sl=-20;}
else{$sl=$xclen2;}
if($maxdepthplot eq "na"){$ymax=$z0*-1;}
else{$ymax=$maxdepthplot;}
$xmin=$sl;
$xmax=$l;
$ymin=-5;
$xdist=$xmax-$xmin;
$ydist=$ymax-$ymin;
$xyratio=$xdist/$ydist;
$xplotmax=35;
$yplotmax=-1*($ydist*$xplotmax)/$xdist;
$shiftmap=-1*$yplotmax+2;
$R3="-R$sl/$l/$ymin/$ymax";
$J3="-JX$xplotmax/$yplotmax";
$B3='-B20g20:Distance:/20g10:Depth:WSen:."$n Lat=$sftrlat[$n] Lon=$sftrlon[$n] Az=$sftraz[$n]":';
$projectmin=$sl-20;
$projectmax=$l+20;
$pL="-L$projectmin/$projectmax";
$pW="-W-30/30";
$pWn="-W-10/10";
$pWt="-W-2/2";

# Project and plot data for each location along the trench on each xsec
$n=0;
foreach(@sftrlon) {
    
    if ($caz ne "na"){$sftraz[$n]=$caz;}
    ## Print statement to keep track of progress
    print "Generating cross-section [$n]: Lon = $sftrlon[$n] Lat = $sftrlat[$n] Str = $sftrstr[$n] Az = $sftraz[$n] \n";

    if($sftraz[$n]>360){$sftraz[$n]=$sftraz[$n]-360;}
    ## Create a basic basemap for each xsec
    `$gmt psbasemap $R3 $J3 -B50g25:Distance:/50g25:Depth:WSen:."$n Lat=$sftrlat[$n] Lon=$sftrlon[$n] Az=$sftraz[$n]": --MAP_GRID_PEN_PRIMARY=lightgrey -K > $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;

        # ___________________________________________ UNFILTERED _________________________________________________________________
        
    # PLOT UNFILTERED INPUT (ERS EQS and BA)
    `$gmt project $indataEQu -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pW $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inEQu.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inEQu.dat $R3 $J3 -Sc0.35 -W0.08,grey -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt project $indataERu -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pW $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inERu.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inERu.dat $R3 $J3 -Sc0.35 -W0.08,grey -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt project $indataBAu -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pW $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inBAu.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inBAu.dat $R3 $J3 -Sc0.35 -W0.08,grey -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;

    # PLOT AS and TO and AA
    `$gmt project $indataASu -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pW $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inASu.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inASu.dat $R3 $J3 -Sc0.35 -W0.08,grey -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt project $indataTOu -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pW $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inTOu.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inTOu.dat $R3 $J3 -Sc0.35 -W0.08,grey -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt project $indataAAu -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pW $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inAAu.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inAAu.dat $R3 $J3 -Sc0.35 -W0.08,grey -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt project $indataRFu -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pW $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inRFu.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inRFu.dat $R3 $J3 -Sc0.35 -W0.08,grey -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;

    # ___________________________________^^^_____ UNFILTERED _______^^^_______________________________________________________

    # PLOT USED INPUT (ERS EQS and BA)
    `$gmt project $indataEQ -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pW $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inEQ.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inEQ.dat $R3 $J3 -Sd0.25 -Gblack -W0.05,black -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt project $indataER -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pW $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inER.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inER.dat $R3 $J3 -Sd0.25 -Gorange -W0.05,black -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt project $indataBA -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pW $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inBA.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inBA.dat $R3 $J3 -Sd0.25 -Ggrey -W0.05,black -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    # PLOT EXTENDED DEPTHS
    `$gmt project $tiltsurf -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pWt $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tiltsurf.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tiltsurf.dat $R3 $J3 -Sc0.1 -W0.05,red -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    
    # PLOT SUPPLEMENTAL EXTENDED DEPTHS
    if ($tilted2 ne "na"){
    `$gmt project $tiltsurf2 -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pWt $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tiltsurf2.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tiltsurf2.dat $R3 $J3 -Sc0.1 -W0.05,green -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    }
    
    # PLOT SUPPLEMENTAL EXTENDED DEPTHS
    if ($tilted3 ne "na"){
    `$gmt project $tiltsurf3 -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pWt $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tiltsurf3.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tiltsurf3.dat $R3 $J3 -Sc0.1 -W0.05,cyan -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    }
    
    # PLOT fillers
    `$gmt project $posshift2 -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pWn $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_posshift2.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_posshift2.dat $R3 $J3 -Sc0.6 -Gturquoise -W0.05,black -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt project $preshift2 -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pWn $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_preshift2.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_preshift2.dat $R3 $J3 -Sc0.3 -Gmagenta -W0.05,black -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    
    # PLOT NODES
    `$gmt project $posshift -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pWn $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_posshift.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_posshift.dat $R3 $J3 -Sc0.5 -Gblue -W0.05,black -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt project $preshift -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pWn $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_preshift.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_preshift.dat $R3 $J3 -Sc0.25 -Gyellow -W0.05,black -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    # PLOT AS and TO and AA
    `$gmt project $indataAS -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pW $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inAS.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inAS.dat $R3 $J3 -Sd0.25 -Gred -W0.05,black -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt project $indataTO -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pW $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inTO.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inTO.dat $R3 $J3 -Sd0.25 -Ggreen -W0.05,black -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt project $indataAA -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pW $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inAA.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inAA.dat $R3 $J3 -Sd0.25 -Gpink -W0.05,black -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt project $indataRF -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pW $pL -Q -S -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inRF.dat`;
    `$gmt psxy $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_inRF.dat $R3 $J3 -Sd0.25 -Gcyan -W0.05,black -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    
    ## Project and plot Slab1.0
    if ($tgs ne "na"){
    `$gmt project -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pL -G10 -Q | awk '{if (\$1<0) {print 360+\$1, \$2} else {print \$1, \$2}}' | $gmt grdtrack -G$tgs > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tgs.dat`;
    `$gmt project $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tgs.dat -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pL $pW -Q -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tgs2.dat`;
    `awk '{print \$1,-\$2}' $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tgs2.dat | $gmt psxy -R -J -W1.75,green -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `awk '{print \$1,-\$2}' $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tgs2.dat | $gmt psxy -R -J -Sc0.15 -Ggreen -W0.25,black -O -K >>  $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    }
    
    ## Project and plot extra grid
    if ($guide ne "na"){
    `$gmt project -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pL -G10 -Q | awk '{if (\$1<0) {print 360+\$1, \$2} else {print \$1, \$2}}' | $gmt grdtrack -G$guide > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tguide.dat`;
    `$gmt project $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tguide.dat -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pL $pW -Q -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tguide2.dat`;
    `awk '{print \$1,-\$2}' $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tguide2.dat | $gmt psxy -R -J -W1.75,cyan -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `awk '{print \$1,-\$2}' $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_tguide2.dat | $gmt psxy -R -J -Sc0.15 -Gcyan -W0.25,black -O -K >>  $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    }

    ## Project and plot raw grid
    `$gmt project -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pL -G10 -Q | awk '{if (\$1<0) {print 360+\$1, \$2} else {print \$1, \$2}}' | $gmt grdtrack -G$raw > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_traw.dat`;
    `$gmt project $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_traw.dat -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pL $pW -Q -Fpz > $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_traw2.dat`;
    `awk '{print \$1,-\$2}' $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_traw2.dat | $gmt psxy -R -J -W1.75,red -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `awk '{print \$1,-\$2}' $folderloc/$ID\_slab2_xsecs_$folder/indata/$ID\_$n\_traw2.dat | $gmt psxy -R -J -Sc0.15 -Gred -W0.25,black -O -K >>  $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    
    ## Project and plot bathymetry data
    `$gmt project -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] -G10 $pL -Q | awk '{print \$1, \$2}' | $gmt grdtrack -G$GEBCO | $gmt project -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] $pL $pW -Q -Fpz | awk '{print \$1,-\$2/1000}' | $gmt psxy -R -J -W2.5,blue -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;

       ## Add a legend to overview map
    #`$gmt pslegend -Dx0.2/0.2/3.5c/BL -F+gwhite+p0.75,black -O -K << EOF >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps
    #N 1
    #S 0.2c d 0.25 black 0.25p 0.6c EQ
    #S 0.2c d 0.25 orange 0.25p 0.6c ER
    #S 0.2c d 0.25 red 0.25p 0.6c AS
    #S 0.2c d 0.25 green 0.25p 0.6c TO
    #S 0.2c c 0.3 yellow 0.25p 0.6c PREshift
    #S 0.2c c 0.5 blue 0.25p 0.6c POSTshift
    #S 0.2c c 0.25 blue 0.25p 0.6c BA
    #S 0.2c c 0.25 darkgreen 0.25p 0.6c BA
    #S 0.2c c 0.25 darkgrey 0.25p 0.6c BA
    #EOF`;

    ## Add an overview map
    `$gmt psbasemap -JM15 -R$minlon/$maxlon/$minlat/$maxlat -Bx5g -By5g -BSWen -X14.5i -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt grdimage $GEBCO --MAP_GRID_PEN_PRIMARY=lightgrey -JM15 -R$minlon/$maxlon/$minlat/$maxlat -O -K -Crelief >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    
    `$gmt psclip $clipmask -JM15 -R$minlon/$maxlon/$minlat/$maxlat -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt grdimage $raw -C$depcpt -JM15 -R$minlon/$maxlon/$minlat/$maxlat -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt grdcontour $raw -C20 -A40 -JM15 -R$minlon/$maxlon/$minlat/$maxlat -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt psclip -C -JM15 -R$minlon/$maxlon/$minlat/$maxlat -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt psxy -W2p $clipmask -P -JM15 -R$minlon/$maxlon/$minlat/$maxlat -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `awk '{print \$1,\$2}' $outtrench | $gmt psxy -R$minlon/$maxlon/$minlat/$maxlat -JM15 -W1.75,black -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;

    ## Testing azimuth direction
    `$gmt project -C$sftrlon[$n]/$sftrlat[$n] -A$sftraz[$n] -G10 $pL -Q -S | awk '{print \$1, \$2}' | $gmt psxy -R -J -Sc0.15 -Gwhite -O -K >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;
    `$gmt pscoast -J -R -O -K -Dl -W0.75p,100 -Bx10g5 -By10g5 >> $folderloc/$ID\_slab2_xsecs_$folder/$ID\_$n\_csec.ps`;

    ## Loop through to next lon,lat position
    $n++;
}
#`rm -r $folderloc/$ID\_slab2_xsecs_$folder/indata`;
unlink "$folderloc/$ID\_slab2_xsecs_$folder/$ID\_finaltrench.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/$ID\_trench.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataAA.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataAAu.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataAS.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataASu.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataBA.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataBAu.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataCP.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataCPu.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataEQ.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataEQu.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataER.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataERu.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataRF.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataRFu.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataTO.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/indataTOu.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/posShift.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/preShift.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/rezultz1.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/rezultz2.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/tilted.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/temp.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/posShift2.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/preShift2.dat";
unlink "$folderloc/$ID\_slab2_xsecs_$folder/dep.cpt";

#####################
sub round_to_nths {
    my ($num, $n) = @_;
    (int $num*$n)/$n
}
######################
