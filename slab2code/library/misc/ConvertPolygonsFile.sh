#!/bin/bash
# KLH 10/10/2019
# Run this script to use any changes made to slab_polygons2.txt in the slab2 code
# i.e., use to convert slab_polygons2.txt to slab_polygons.txt
# slab_polygons2.txt is used for plotting in GMT, whereas slab_polygons.txt is used in the slab2 code
# To run type: bash ConvertPolygonsFile.sh

# Move all rows into one row:
awk '{printf("%s%s",$0,NR?", ":"\n")}' slab_polygons2.txt > slab_polygons3.txt
# Remove the first instance of '>':
awk 'gsub(/^[>]+/,"")' slab_polygons3.txt > slab_polygons4.txt
# Replace all other instances of '>' with a new line:
sed -e $'s/>/\\\n/g' slab_polygons4.txt > slab_polygons5.txt
# Remove the last comma on each row:
sed 's/,$//' slab_polygons5.txt > slab_polygons6.txt
sed 's/, $//' slab_polygons6.txt > slab_polygons.txt
# Remove the files created with this script that are not needed:
rm -rf slab_polygons3.txt slab_polygons4.txt slab_polygons5.txt slab_polygons6.txt