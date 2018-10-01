# coding=utf-8

import os, sys
import numpy as np
import glob
import string
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import cm 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from pylab import *
import datetime
from pprint import pprint
from netCDF4 import Dataset, datetime, date2num,num2date
from scipy.ndimage.filters import gaussian_filter
import ogr
import osr


__author__   = 'Trond Kristiansen'
__email__    = 'me (at) trondkristiansen.com'
__created__  = datetime(2016, 8, 10)
__modified__ = datetime(2016, 8, 10)
__version__  = "1.0"
__status__   = "Production"

# --------
# calculateaverages.py
#
# This script takes the output from opendrift and calculates area averages 
# as a function of time. The total area around North Sea is divided into
# bins of specific resolution and the number of particles within each bin 
# is summed for a specific time period (e.g. 1 month). The total output 
# is a heatmap of where the most particles reside for each time period.
# --------

def createBins(requiredResolution):

    print 'func: createBins() => Creating bins for averaging'
    xmin=17.2; xmax=19.3
    ymin=69.5; ymax=70.4
    
    deg2rad=np.pi/180.
    R = 6371  # radius of the earth in km
    # Distance from minimum to maximim longitude
    x = (xmax*deg2rad - xmin*deg2rad) * cos( 0.5*(ymax*deg2rad+ymax*deg2rad) )
    y =  ymax*deg2rad - ymax*deg2rad
    dx = R * sqrt( x*x + y*y )
    print "Distance from minimum to maximim longitude binned area is %s km"%(dx)

    # Distance from minimum to maximim latitude
    x = (xmax*deg2rad - xmax*deg2rad) * cos( 0.5*(ymax*deg2rad+ymin*deg2rad) )
    y =  ymax*deg2rad - ymin*deg2rad
    dy = R * sqrt( x*x + y*y )

    print "Distance from minimum to maximim latitude binned area is %s km"%(dy)


    ngridx = int(np.round(dx/requiredResolution,0))
    ngridy = int(np.round(dy/requiredResolution,0))
    
    xi = np.linspace(np.floor(xmin),np.ceil(xmax),ngridx)
    yi = np.linspace(np.floor(ymin),np.ceil(ymax),ngridy)

    print '=> created binned array of domain of size (%s,%s) with resolution %s'%(ngridx,ngridy,requiredResolution)

    return xi,yi

def calculateAreaAverages(xi,yi,cdf,first,weeksInYear):

    print 'func: calculateAreaAverages() => Calculating averages within bins'
    print '=> binned domain (%2.1f,%2.1f) to (%2.1f,%2.1f)'%(np.min(xi),np.min(yi),np.max(xi),np.max(yi))

    timesteps = cdf.variables['time'][:]
    timeunits = cdf.variables["time"].units
    
    print '=> found %s timesteps in input file'%(len(timesteps))
    newWeek=-9

    for tindex, t in enumerate(timesteps): 

        currentDate = num2date(t, units=timeunits, calendar="gregorian")
        year,weekNumber,DOW = currentDate.isocalendar()
       
        Xpos = cdf.variables['lon'][:,tindex]
        Ypos = cdf.variables['lat'][:,tindex]

        H, xedges, yedges = np.histogram2d(Xpos, Ypos, bins=(xi, yi), normed=False)

        if (tindex==0 and first is True):
            weeklyFrequency=np.zeros((weeksInYear,np.shape(H)[0],np.shape(H)[1]), dtype=float32)
            
        if weekNumber != newWeek:
            print "=> Adding data to week: %s (startdate: %s)"%(weekNumber,currentDate)
            newWeek=weekNumber
        weeklyFrequency[weekNumber,:,:]=weeklyFrequency[weekNumber,:,:] + H
        
    # Create log values and levels for frequencyplot
    weeklyFrequency=ma.log(weeklyFrequency)
    levels = np.arange(weeklyFrequency.min(),weeklyFrequency.max(),(weeklyFrequency.max()-weeklyFrequency.min())/10)
        
    sigma = 0.2 # this depends on how noisy your data is, play with it!
    first=False

    return gaussian_filter(weeklyFrequency, sigma), first

def plotDistribution(shapefile,kelpData,week,baseout,xii,yii,polygonIndex,distType):
    print "Plotting the distributions for week: %s"%(week)
    plt.clf()
    plt.figure(figsize=(10,10), frameon=False)
    ax = plt.subplot(111)

    mymap = Basemap(llcrnrlon=17.2, llcrnrlat=69.5,
                       urcrnrlon=19.3, urcrnrlat=70.4,
                       resolution='h', projection='merc',lon_0=18,lat_0=60,area_thresh=0.)

    x, y = mymap(xii,yii)
       

    levels=np.arange(np.min(kelpData),np.max(kelpData),0.5)
                              
    CS1 = mymap.contourf(x,y,np.fliplr(np.rot90(kelpData,3)),levels,cmap=cm.get_cmap('Spectral_r',len(levels)-1), extend='both',alpha=1.0)
    plt.colorbar(CS1,orientation='vertical',extend='both', shrink=0.5)

    mymap.drawcoastlines()
    mymap.fillcontinents(color='grey',zorder=2)
    mymap.drawcountries()
    mymap.drawmapboundary()

    plt.title('Kelp week: %s'%(week))
    print "Adding polygons to plot"

    mypatches=createPathsForPolygons(shapefile,mymap)
    p = PatchCollection(mypatches,alpha=1.0,facecolor='none',lw=1.0,edgecolor='purple',zorder=2)
    ax.add_collection(p)
  
    if polygonIndex=="all":
        plotfile=baseout+'/Kelp_distribution_polygon_'+str(polygonIndex)+'_fullperiod.png'
    else:
        plotfile=baseout+'/Kelp_distribution_polygon_'+str(polygonIndex)+'_week_'+str(week)+'.png'
    print "=> Creating plot %s"%(plotfile)             
    plt.savefig(plotfile,dpi=300)
                                                                

def getPathForPolygon(ring,mymap):
    codes=[]
    x = [ring.GetX(j) for j in range(ring.GetPointCount())]
    y = [ring.GetY(j) for j in range(ring.GetPointCount())]
    codes += [mpath.Path.MOVETO] + (len(x)-1)*[mpath.Path.LINETO]
    
    pathX,pathY=mymap(x,y)
    mymappath = mpath.Path(np.column_stack((pathX,pathY)), codes)

    return mymappath
   
def createPathsForPolygons(shapefile,mymap):

    mypatches=[]
    s = ogr.Open(shapefile)
    print shapefile

    for layer in s:
        # get projected spatial reference
        sr = layer.GetSpatialRef()
        # get geographic spatial reference
        geogr_sr = sr.CloneGeogCS()
        # define reprojection
        proj_to_geog = osr.CoordinateTransformation(sr, geogr_sr)

        polygons=[x+1 for x in xrange(layer.GetFeatureCount()-1)]
        for polygonIndex,polygon in enumerate(polygons):
            feature = layer.GetFeature(polygonIndex)
            geom = feature.GetGeometryRef()
            points = geom.GetGeometryCount()
            ring = geom.GetGeometryRef(0)
            geom.Transform(proj_to_geog)

            if ring.GetPointCount() > 3:
                polygonPath = getPathForPolygon(ring,mymap)
                path_patch = mpatches.PathPatch(polygonPath, lw=2, edgecolor="purple",facecolor='none')
                                
                mypatches.append(path_patch)

    return mypatches

def main():

    # EDIT --------------------------------------
    # Which species to calculate for
    
    # The timespan part of the filename
    startdate='01052016'
    enddate='01062016'

    # Results and storage folders
    base='results'
    baseout='distributionFigures'
    
    if not os.path.exists(baseout): os.makedirs(baseout)

    # The resolution of the output grid in kilometers
    requiredResolution = 0.200 # km between each binned box

    # END EDIT ----------------------------------

    # Create the grid you want to calculate frequency on
    xi,yi = createBins(requiredResolution)
    weeksInYear=52
    xii,yii=np.meshgrid(xi[:-1],yi[:-1])
    firstRead=True
    hexagon=False
    maxNumberOfPolygons=1

    shapefile='/Users/trondkr/Dropbox/NIVA/KelpFloat/Kelp/Shapefile/ShapefilesHGU/kelpexpol_exp_grazed_combined.shp'
        
    print "=> Using shapefile %s"%(shapefile)
    s = ogr.Open(shapefile)
    for layer in s:
        polygons=[x+1 for x in xrange(layer.GetFeatureCount()-1)]
    
 #   polygons=[2]

    if firstRead:
        allData=np.zeros((maxNumberOfPolygons,weeksInYear,len(xi)-1,len(yi)-1))
        print "=> Created final array for all data of size :",np.shape(allData)
        firstRead=False

    for polygonIndex,polygon in enumerate(polygons):
        first=True
            
        feature = layer.GetFeature(polygonIndex)
        geom = feature.GetGeometryRef()
        points = geom.GetGeometryCount()
        ring = geom.GetGeometryRef(0)
        if ring.GetPointCount() > 3:

            infile=base+'/Kelp_polygon_'+str(polygon)+'_kelp_opendrift_'+str(startdate)+'_to_'+str(enddate)+'_novertical.nc'
            print "=> Opening input file: %s"%(os.path.basename(infile))

            if os.path.exists(infile):
                cdf = Dataset(infile)
                filteredData, first = calculateAreaAverages(xi,yi,cdf,first,weeksInYear)
           
                allData[polygonIndex,:,:,:]=filteredData
            else:
                print "==>> Input file %s could not be opened"%(infile)
    
    for week in xrange(1,weeksInYear,1):
        # Calculate the cumulative distribution for each month and species
        first=True
        for polygonIndex,polygon in enumerate([x+1 for x in xrange(len(polygons))]):
            if first:
                kelpData=np.zeros((len(xi)-1,len(yi)-1))
                first=False
                print "==> Created array of data for week: ",week," with size: ",np.shape(kelpData)

            kelpData=kelpData + np.squeeze(allData[polygonIndex,week,:,:])
        levels=np.arange(np.min(kelpData),np.max(kelpData),0.5)
      
        if (len(levels)>2 and np.mean(kelpData)>0):
            plotDistribution(shapefile,kelpData,week,baseout,xii,yii,polygons[0],"weekly")

    # Plot the distribution for all weeks
    plotDistribution(shapefile,kelpData,week,baseout,xii,yii,"all","integrated")

if __name__ == "__main__":
    main()

