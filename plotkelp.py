#!/usr/bin/env python

from datetime import datetime, timedelta
import numpy as np

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift
from kelp.kelpClass import PelagicPlanktonDrift
from pprint import pprint
from netCDF4 import Dataset, datetime, date2num, num2date
from scipy.ndimage.filters import gaussian_filter
import ogr
import osr
import matplotlib
import os

o = PelagicPlanktonDrift(loglevel=0)  # Set loglevel to 0 for debug information

#######################
# Preparing readers
#######################

base = 'results'
baseout = 'figures'
hexagon = False
startdate = '01052016'
enddate = '01082016'
experiment = 1

reader_basemap = reader_basemap_landmask.Reader(
    llcrnrlon=17.0, llcrnrlat=69.4,
    urcrnrlon=18.0, urcrnrlat=70.0,
    resolution='h', projection='merc')

shapefile = '/Users/trondkr/Dropbox/NIVA/KelpFloat/Kelp/Shapefile/ShapefilesHGU/kelpexpol_exp_grazed_combined.shp'

print "=> Using shapefile %s" % (shapefile)
s = ogr.Open(shapefile)
for layer in s:
    polygons = [x + 1 for x in xrange(layer.GetFeatureCount() - 1)]

    for ii, polygon in enumerate(polygons):
        polygonIndex = int(polygon) - 1
        feature = layer.GetFeature(polygonIndex)
        geom = feature.GetGeometryRef()
        points = geom.GetGeometryCount()
        ring = geom.GetGeometryRef(0)
        print ring.GetPointCount()

        if ring.GetPointCount() > 3:
            specie = 'Kelp'
            filename = 'results/%s_polygon_%s_experiment_%s_%s_to_%s.nc' % (
            specie, polygon, experiment, startdate, enddate)
            plotfilename = 'figures/%s_polygon_%s_experiment_%s_%s_to_%s.png' % (
            specie, polygon, experiment, startdate, enddate)
            plotfilenameColor = 'figures/%s_polygon_%s_experiment_%s_%s_to_%s_color.png' % (
            specie, polygon, experiment, startdate, enddate)
            plotfilenameVertical = 'figures/%s_polygon_%s_experiment_%s_%s_to_%s_vertical.png' % (
            specie, polygon, experiment, startdate, enddate)
            plotfilenameAnime = 'figures/%s_polygon_%s_experiment_%s_%s_to_%s.mp4' % (
            specie, polygon, experiment, startdate, enddate)
            print filename
            if os.path.exists(filename):
                print "=> Opening input file: %s" % (filename)

                o.io_import_file(filename)
                # o.add_reader([reader_basemap]) #Do not include basemap when stranding is deactivated

                # o.plot_vertical_distribution()
                # o.plot(linecolor='z',lvmin=-150, lvmax=0,filename=plotfilenameColor)
                # o.plot(filename=plotfilename)
                o.animation(filename=plotfilenameAnime)
