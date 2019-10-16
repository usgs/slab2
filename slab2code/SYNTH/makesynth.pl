#!/usr/bin/perl
#####
use Math::Trig;
#####
$outfile="exp_10-19_input.csv";
open(OUT,">$outfile");
#####
$unc=15;
$id0=100000;
$etype="EQ";
$mag=5.5;
$time="2000-01-01 01:00:00";
$string="nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan";
print OUT "lat,lon,depth,unc,ID,etype,mag,time,Paz,Ppl,Taz,Tpl,S1,D1,R1,S2,D2,R2,mlon,mlat,mdep,id_no,src\n";
#####
$lon0=126.4;
$lat0=1.3;
$az=15;
$dipdir=$az+90;
$len1=-300;
$len2=300;
$dx=20;
#####
# Global AS properties
$asdist=80;
$asdip=8;
$asdip=$asdip*(pi/180);
$tdep=7.5;
$as_end_dep=$tdep+($asdist*tan($asdip));
#####
# SZ properties
$szdist=100;
$szdip=20;
$szdip=$szdip*(pi/180);
$sz_end_dep=$as_end_dep+($szdist*tan($szdip));
$dsz=5;
$ddep1=$dsz*tan($szdip);
#####
# WBZ properties
$wbzdist=200;
$wbzdip=50;
$wbzdip=$wbzdip*(pi/180);
$ddep2=$dsz*tan($wbzdip);
#####
@tlocs=`gmt project -C$lon0/$lat0 -A$az -G$dx -Q -L$len1/$len2`;
#####
$n=0;
$id=$id0;
foreach(@tlocs){
    ($lon1,$lat1,$a)=split;
    print OUT "$lat1,$lon1,$tdep,$unc,$id,$etype,$mag,$time,$string\n";
    $id++;
####
    @tloc0=`gmt project -C$lon1/$lat1 -A$dipdir -L0/$asdist -G$asdist -Q`;
    ($lon2,$lat2,$a)=split(/\t/,$tloc0[1]);
    print OUT "$lat2,$lon2,$as_end_dep,$unc,$id,$etype,$mag,$time,$string\n";
    $id++;
####
    @szloc=`gmt project -C$lon2/$lat2 -A$dipdir -L0/$szdist -G$dsz -Q`;
    $cdep=$as_end_dep;
    $nn=0;
    foreach(@szloc){
	($lon3,$lat3,$a)=split;
	if($nn>0){
	    $cdep=$cdep+$ddep1;
	    print OUT "$lat3,$lon3,$cdep,$unc,$id,$etype,$mag,$time,$string\n";
	    $id++;
	}
	$nn++;
    }
####
    @wbzloc=`gmt project -C$lon3/$lat3 -A$dipdir -L0/$wbzdist -G$dsz -Q`;
    $nn=0;
    foreach(@wbzloc){
        ($lon4,$lat4,$a)=split;
        if($nn>0){
            $cdep=$cdep+$ddep2;
            print OUT "$lat4,$lon4,$cdep,$unc,$id,$etype,$mag,$time,$string\n";
            $id++;
        }
        $nn++;
    }
####
    $n++;
}
