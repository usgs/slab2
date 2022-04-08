#!/bin/bash
# Get slab2 contours and save as shapefile
# This makes contours from the grd files, then converts the contours into ESRI shapefiles
# Make sure to activate the slab2env (slab2 conda env) or have the GDAL package installed

if [ $# -eq 0 ]; then
echo " "
echo "How to use: bash MakeKML.sh Region Data Date ContourInterval"
echo "GDAL is Requied to make the KML file. The slab2 anaconda environemtn contains GDAL, so activate that enviornment before running this script."
echo "Run from directory where the slab2 file(s) is located."
echo "Region can be 'all' to run for all 27 slab2 regions"
echo "Region can also be the 3 letter acronym for a specific slab2 region (e.g., 'alu')"
echo "Data is the data you want to convert - depth (dep), dip, strike (str), uncertainty (unc), thickness (thk). Use the 3 letter acronym."
echo "Date is the date the slab model was generated (e.g., 04.08.22). This date is included in the name of the file (alu_slab2_dep_04.08.22.grd)."
echo "ContourInterval is the contour interval in km"
echo "Example: bash MakeKML.sh alu dep 04.08.22 20"
echo "copy and paste the above example into the terminal to genterate Shapefiles and a KML for the alu region"
echo " "
exit
fi

region=$1
data=$2
date=$3
ci=$4

if [ $region == "all" ]
then
    # Set slab array to run this for all 27 slab regions:
    declare -a slabArray=("alu" "cal" "cam" "car" "cas" "cot" "hal" "hel" "him" "hin" "izu" "ker" "kur" "mak" "man" "mue" "pam" "png" "phi" "puy" "ryu" "sam" "sco" "sol" "sul" "sum" "van")
else
    # If running for only 1 region
    declare -a slabArray=($region)
fi

for slab in "${slabArray[@]}"; do
    echo "Getting contours for ${slab} and making shapefile and KML file"
    
    # create region directory:
    mkdir ${slab}

    # use GMT to create a contour file:
    gmt grdcontour "${slab}"_slab2_"${data}"_"${date}".grd -C"${ci}" -D"${slab}"_slab2_"${data}"_"${date}"_cont.xyz
    
    # Make shapefile:
    gmt gmtconvert "${slab}"_slab2_"${data}"_"${date}"_cont.xyz -fg -aZ=contour:double+GLINE > "${slab}"_slab2_"${data}"_"${date}"_cont-01.xyz
    ogr2ogr -f "ESRI shapefile" "${slab}".shp "${slab}"_slab2_"${data}"_"${date}"_cont-01.xyz
    
    # Make kml file:
    gmt gmt2kml "${slab}"_slab2_"${data}"_"${date}"_cont.xyz -Fl -W2p -CGEdepth.cpt > "${slab}"_slab2_"${data}"_"${date}".kml
    
    # clean up files & move files to region directory:
    rm -rf "${slab}"_slab2_"${data}"_"${date}"_cont.xyz
    rm -rf "${slab}"_slab2_"${data}"_"${date}"_cont-01.xyz
    mv ${slab}_* ${slab}
    mv ${slab}.* ${slab}
done

echo "Finished! Files have been moved to the ${slab} folder in your current working directory"
