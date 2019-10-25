# Converts polygon CSV file (columns should be lon, lat) into the proper format for slab_polygons.txt
# Change the CSV input file in line 8 for your polygon file
# Change the three-letter slab region code in line 10

import pandas as pd
import numpy as np

data = pd.read_csv('YourFileHere.csv')

polystr = 'exp, '
for index,row in data.iterrows():
	lon, lat = str(row['lon']), str(row['lat'])
	polystr += lon
	polystr += ', '
	polystr += lat
	polystr += ', '

print(polystr)
print('Now replace the lines printed above to slab_polygons.txt for the slab region being updated. Remove the last comma.')
