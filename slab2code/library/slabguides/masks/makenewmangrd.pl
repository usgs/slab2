#!/usr/bin/perl -w
#####
# SLAB ID
###

$ref="/Users/ginevramoore/Documents/Slab2/Ginevra2017/gmtfiles/grdfiles/man_SG_man.grd";
$gmt="/Users/ginevramoore/GMT-5.3.1.app/Contents/Resources/bin/gmt";
$trenches="trenches_usgs_2016.csv";
#####
$minlon=999;
$maxlon=-999;
$uncmax=0;
$minlat=$minlon;$maxlat=$maxlon;


####
# CREATE SURFACE OF SMOOTHED SLAB2 DATA
###
$grd1ps="$folder/temp_slab_ps.grd";
$inc="-I0.5/0.5";
`$gmt blockmean $psmooth $bounds $inc | $gmt surface -G$grd1ps $inc $bounds -T0.0i -T0.0b -V`;
