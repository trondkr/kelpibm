from netCDF4 import Dataset, datetime, date2num, num2date
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib.pyplot import cm
import etopo1
import pandas as pd
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import ogr
import osr
from math import *
from scipy.ndimage.filters import gaussian_filter
import mpl_util
import laplace_filter
import cmocean

__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@niva.no'
__created__ = datetime(2017, 5, 16)
__modified__ = datetime(2017, 5, 16)
__version__ = "1.0"
__status__ = "Development, modified on 16.05.2017"


def getPathForPolygon(ring, mymap):
    codes = []
    x = [ring.GetX(j) for j in range(ring.GetPointCount())]
    y = [ring.GetY(j) for j in range(ring.GetPointCount())]
    codes += [mpath.Path.MOVETO] + (len(x) - 1) * [mpath.Path.LINETO]

    pathX, pathY = mymap(x, y)
    mymappath = mpath.Path(np.column_stack((pathX, pathY)), codes)

    return mymappath


def createBins(requiredResolution):
    print 'func: createBins() => Creating bins for averaging'
    xmin = 15
    xmax = 21
    ymin = 69
    ymax = 72

    deg2rad = np.pi / 180.
    R = 6371  # radius of the earth in km
    # Distance from minimum to maximim longitude
    x = (xmax * deg2rad - xmin * deg2rad) * cos(0.5 * (ymax * deg2rad + ymax * deg2rad))
    y = ymax * deg2rad - ymax * deg2rad
    dx = R * sqrt(x * x + y * y)
    print "Distance from minimum to maximim longitude binned area is %s km" % (dx)

    # Distance from minimum to maximim latitude
    x = (xmax * deg2rad - xmax * deg2rad) * cos(0.5 * (ymax * deg2rad + ymin * deg2rad))
    y = ymax * deg2rad - ymin * deg2rad
    dy = R * sqrt(x * x + y * y)

    print "Distance from minimum to maximim latitude binned area is %s km" % (dy)

    ngridx = dx / requiredResolution
    ngridy = dy / requiredResolution

    xi = np.linspace(xmin, xmax, ngridx)
    yi = np.linspace(ymin, ymax, ngridy)

    print '=> created binned array of domain of size (%s,%s) with resolution %s' % (ngridx, ngridy, requiredResolution)

    return xi, yi


def calculateAreaAverages(xi, yi, Xpos, Ypos, area):
    print 'func: calculateAreaAverages() => Calculating averages within bins'
    print '=> binned domain (%2.1f,%2.1f) to (%2.1f,%2.1f)' % (np.min(xi), np.min(yi), np.max(xi), np.max(yi))
    print '=> drift domain (%2.1f,%2.1f) to (%2.1f,%2.1f)' % (np.min(Xpos), np.min(Ypos), np.max(Xpos), np.max(Ypos))

    H, xedges, yedges = np.histogram2d(np.asarray(Xpos), np.asarray(Ypos), bins=(xi, yi), normed=False)

    sigma = 0.2  # this depends on how noisy your data is, play with it!
    return gaussian_filter(H, sigma)


def createPathsForPolygons(shapefile, mymap):
    mypatches = []
    s = ogr.Open(shapefile)

    for layer in s:
        # get projected spatial reference
        sr = layer.GetSpatialRef()
        # get geographic spatial reference
        geogr_sr = sr.CloneGeogCS()
        # define reprojection
        proj_to_geog = osr.CoordinateTransformation(sr, geogr_sr)

        polygons = [x + 1 for x in xrange(layer.GetFeatureCount()-1)]

        for polygonIndex, polygon in enumerate(polygons):

            feature = layer.GetFeature(polygonIndex)
            geom = feature.GetGeometryRef()

            ring = geom.GetGeometryRef(0)
            geom.Transform(proj_to_geog)

            if ring.GetPointCount() > 3:
                #print "Looping over polygon index %s with %s points" % (polygonIndex, ring.GetPointCount())
                polygonPath = getPathForPolygon(ring, mymap)
                path_patch = mpatches.PathPatch(polygonPath, lw=0.2, edgecolor="purple", facecolor='none')

                mypatches.append(path_patch)
    return mypatches


def createScatterPlot(shapefile, sedLats, sedLons, sedDepths, sedCats, useEtopo1, etopo1name):
    plt.clf()
    # plt.figure(figsize=(10,10), frameon=False)
    ax = plt.subplot(111)

    if useEtopo1:
        """Get the etopo2 data"""
        e1 = Dataset(etopo1name, 'r')
        lons = e1.variables["lon"][:]
        lats = e1.variables["lat"][:]

        res = etopo1.findSubsetIndices(60, 63, 2, 7, lats, lons)

        lon, lat = np.meshgrid(lons[res[0]:res[1]], lats[res[2]:res[3]])
        print "Extracted data for area Kelp: (%s,%s) to (%s,%s)" % (lon.min(), lat.min(), lon.max(), lat.max())
        bathy = e1.variables["z"][int(res[2]):int(res[3]), int(res[0]):int(res[1])]
        bathySmoothed = laplace_filter.laplace_filter(bathy, M=None)

        etopo1levels = [-100, -75, -65, -50, -35, -25, -15, -10, -5, 0]

    mymap = Basemap(llcrnrlon=17.5, llcrnrlat=69.5,
                    urcrnrlon=21, urcrnrlat=70.5,
                    resolution='h', projection='merc', lon_0=np.mean(np.array(sedLons)),
                    lat_0=np.mean(np.array(sedLats)), area_thresh=0.)

    if useEtopo1:
        xe1, ye1 = mymap(lon, lat)

        CSE1 = mymap.contour(xe1, ye1, bathySmoothed, 10, alpha=1.0, linewidths=0.5)
        print "Depth levels", np.min(bathySmoothed), np.max(bathySmoothed)
    # Prepare and plot the kelp data
    x, y = mymap(sedLons, sedLats)
    levels = np.arange(np.min(sedDepths), np.max(sedDepths), 0.5)

    mypatches = createPathsForPolygons(shapefile, mymap)
    p = PatchCollection(mypatches, alpha=1.0, facecolor='none', lw=2.0, edgecolor='purple', zorder=2)
    ax.add_collection(p)

    sizes = np.array(sedCats) * 10 + 20
    CS1 = mymap.scatter(x, y, s=5, c=sedDepths, cmap=cm.get_cmap('Spectral_r', len(levels) - 1), alpha=0.5, lw=0)
    plt.colorbar(CS1, orientation='vertical', extend='both', shrink=0.5)

    mymap.drawcoastlines()
    mymap.fillcontinents(color='grey', zorder=2)
    mymap.drawcountries()
    mymap.drawmapboundary()

    plt.title('Tare sedimentering')
    #  plt.show()


def plotHistDistribution(shapefile, hist, xii, yii, polygon, experiment, startdate, enddate, plotType):
    plt.clf()

    # plt.figure(figsize=(10,10), frameon=False)
    ax = plt.subplot(111)

    mymap = Basemap(llcrnrlon=17.0, llcrnrlat=69.4,
                    urcrnrlon=19.0, urcrnrlat=70.0,
                    resolution='h', projection='merc', lon_0=np.mean(np.array(xii)),
                    lat_0=np.mean(np.array(yii)), area_thresh=0.)

    xiii, yiii = mymap(xii, yii)
    levels = np.arange(np.min(hist), np.max(hist), 1)
    # levels = np.arange(np.min(hist), 200, 1)
    levels = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]

    CS1 = mymap.contourf(xiii, yiii, np.fliplr(np.rot90(hist, 3)), levels,
                         cmap=mpl_util.LevelColormap(levels, cmap=cmocean.cm.haline),
                         origin='lower',
                         extend='both',
                         alpha=1.0)

    plt.colorbar(CS1, orientation='vertical', extend='both', shrink=0.5)
    mymap.drawcoastlines()
    mymap.fillcontinents(color='grey', zorder=2)
    mymap.drawcountries()
    mymap.drawmapboundary()

    mypatches = createPathsForPolygons(shapefile, mymap)
    # colors = 100 * np.random.rand(len(mypatches))
    # p = PatchCollection(mypatches, alpha=0.9, zorder=10)
    # p.set_array(np.array(colors))
    # ax.add_collection(p)

    p = PatchCollection(mypatches, alpha=0.8, facecolor='none', lw=0.5, edgecolor='red', zorder=10)
    ax.add_collection(p)

    if plotType == "individual":
        p = PatchCollection([mypatches[polygon]], alpha=0.9, facecolor='darkorange', lw=0.5, edgecolor='darkorange',
                            zorder=10)
        ax.add_collection(p)

    if plotType == "all":
        p = PatchCollection(mypatches, alpha=0.9, facecolor='darkorange', lw=0.5, edgecolor='darkorange', zorder=10)
        ax.add_collection(p)

    if plotType == "all":

        x = xii.flatten()
        y = yii.flatten()
        hh = np.fliplr(np.rot90(hist, 3))

        z = hh.flatten()

        np.savetxt("allDensityXYZ.csv", (x, y, z))


    # plt.title('Tare sedimentering')
    print "Adding %s release polygons to map" % (len(mypatches))
    # plt.show()

    if plotType == "individual":
        plotfilename = 'distributionFigures/Kelp_polygon_%s_experiment_%s_%s_to_%s.png' % (
            polygon+1, experiment, startdate, enddate)

    if plotType == "all":
        plotfilename = 'distributionFigures/Kelp_allPolygons_experiment_%s_%s_to_%s.png' % (
            experiment, startdate, enddate)

    print "=> Creating plot %s" % plotfilename
    plt.savefig(plotfilename, dpi=300)


def calculateLevelOfSedimentation(filename, densities, area):
    print "Filename: %s" % (filename)
    cdf = Dataset(filename)

    z = cdf.variables['z'][:]
    # Depth is positive
    h = - (cdf.variables['sea_floor_depth_below_sea_level'][:])

    time = cdf.variables['time'][:]
    trajectory = cdf.variables['trajectory'][:]
    density = cdf.variables['plantpart'][:]

    lat = cdf.variables['lat'][:]
    lon = cdf.variables['lon'][:]

    sedLats = []
    sedLons = []
    sedCats = []
    sedDepths = []

    for tr in xrange(len(trajectory[:])):

        dens = int(np.ma.mean(density[tr, :]))
        diffDepthtoBottom = np.abs(np.squeeze(z[tr, :]) - np.squeeze(h[tr, :]))
        ind = np.where(diffDepthtoBottom < 0.05)
        diff = 0.05  # Difference in movement in meters

        for pos in xrange(len(z[tr, :-2])):

            if ((z[tr, pos] - z[tr, pos + 1] < diff) and (z[tr, pos + 1] - z[tr, pos + 2] < diff) and diffDepthtoBottom[pos] < diff):

                #  print "Found index %s depth %s seafloor %s"%(pos,z[tr,pos],h[tr,pos])

                if (ind[0].size > 0 and dens in densities):
                    # currentIndex = ind[0][0]
                    sedLats.append(lat[tr, pos])
                    sedLons.append(lon[tr, pos])
                    sedCats.append(density)
                    sedDepths.append(z[tr, pos])
                    break

    sedimentrate = ((len(sedLons) - len(trajectory)) / (len(trajectory) * 1.0)) * 100.
    print "Found %s positions sedimented out of %s (%s percent)" % (len(sedLons), len(trajectory), sedimentrate)

    # plt.plot(time/3600.,z[tr,:],color=colors[ind])
    return sedLats, sedLons, sedDepths, sedCats, len(sedLons), len(trajectory)
    # plt.show()


    # plt.show()


useEtopo1 = True
densities = [0, 1, 2]
first = True

plotCumulative = True
plotScatter = False
etopo1name = '/Users/trondkr/Dropbox/Projects/arcwarm/maps/ETOPO1_Ice_g_gmt4.grd'
requiredResolution = 0.16# km between each binned box
base = 'results'
shapefile = '/Users/trondkr/Dropbox/NIVA/KelpFloat/Kelp/Shapefile/ShapefilesHGU/kelpexpol_exp_grazed_combined.shp'
experiment = 1

if experiment == 1:
    startdate = '01052016'
    enddate = '01082016'
if experiment == 2:
    startdate = '01032016'
    enddate = '15052016'
if experiment == 3:
    startdate = '01112015'
    enddate = '01042016'
if experiment == 4:
    startdate = '01052016'
    enddate = '01082016'
s = ogr.Open(shapefile)
for layer in s:
    polygons = [x + 1 for x in xrange(layer.GetFeatureCount()-1)]
   # polygons=[16,17]
    totalSedimentation = 0;
    totalParticles = 0

    for polygonIndex, polygon in enumerate(polygons):

        filename = 'results/Kelp_polygon_%s_experiment_%s_%s_to_%s.nc' % (polygon, experiment, startdate, enddate)

        # plotVerticalDistribution(filename)

        sedLats, sedLons, sedDepths, sedCats, ts, tp = calculateLevelOfSedimentation(filename, densities, polygonIndex)
        totalSedimentation = totalSedimentation + ts
        totalParticles = totalParticles + tp

        print "Total particles %s sedimented %s" % (totalParticles, totalSedimentation)
        print "Running:"
        print "=> Scatter on:%s Cumulative: %s" % (plotScatter, plotCumulative)

        if plotCumulative:
            xi, yi = createBins(requiredResolution)
            xii, yii = np.meshgrid(xi[:-1], yi[:-1])

            if first:
                allData = np.zeros((len(polygons), len(xi) - 1, len(yi) - 1))
                print "=> Created final array for all data of size :", np.shape(allData)
                first = False

            hist = calculateAreaAverages(xi, yi, sedLons, sedLats, polygonIndex)

            plotHistDistribution(shapefile, hist, xii, yii, polygonIndex, experiment, startdate, enddate, "individual")

            allData[polygonIndex, :, :] = hist

        if plotScatter:
            createScatterPlot(shapefile, sedLats, sedLons, sedDepths, sedCats, useEtopo1, etopo1name)

    if plotCumulative:

        sedimentrate = ((totalSedimentation - totalParticles) / (totalParticles * 1.0)) * 100.
        print "Found %s positions sedimented out of %s (%s percent)" % (
        totalSedimentation, totalParticles, sedimentrate)

        # Calculate the cumulative distribution for each month and species
        first = True
        for polygonIndex, polygon in enumerate([x + 1 for x in xrange(len(polygons))]):
            if first:
                kelpData = np.zeros((len(xi) - 1, len(yi) - 1))
                first = False
                print "==> Created array of data for polygon: ", polygonIndex, " with size: ", np.shape(kelpData)

            kelpData = kelpData + np.squeeze(allData[polygonIndex, :, :])
        levels = np.arange(np.min(kelpData), np.max(kelpData), 0.5)

        # Plot the distribution for all weeks
        plotHistDistribution(shapefile, kelpData, xii, yii, polygonIndex, experiment, startdate, enddate, "all")
