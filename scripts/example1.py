# -*- coding: utf-8 -*-

import geopandas as gpd

# make a test-model for an area of 10x10 km
extent = [140000, 150000, 450000, 460000]

# read the sources wihin this extent
bbox = [extent[0], extent[2], extent[1], extent[3]]
gdf = gpd.read_file('..\data\waterlopen_project.shp', bbox=bbox)

# generate a modflow model for this extent
delr = 1000
