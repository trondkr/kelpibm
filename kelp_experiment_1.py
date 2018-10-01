#!/usr/bin/env python

from datetime import datetime, timedelta
import numpy as np
import time

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_ROMS_native
from kelp.kelpClass import PelagicPlanktonDrift
import os
from netCDF4 import Dataset, datetime, date2num, num2date
import random
import math
import glob
from random import randint

try:
    from osgeo import gdal, osr, ogr
except Exception as e:
    print e
    raise ValueError('OGR library is needed to read shapefiles.')


def commonKelpProperties():
    # Old blade / lamina = 0
    # New blade / lamina = 1
    # Stipe     = 2
    # Fragment  = 3

    weights = [0.373, 0.7122, 0.4684, 0.03749, 5.908e-6]
    areas = [0.189097, 0.1323, 0.000452, 0.0007549, 4.39e-6]
    diameters = [0.4887, 0.410, 0.0240, 0.09807, 0.002365]
    lengths = [0.002, 0.002, 0.984, 0.0031, 0.002365]
    volumes = (np.asarray(areas) * np.asarray(lengths)).tolist()
    densities = [1446.6, 1541.0, 1882, 234.99, 1035]
    SDdensities = [401, 668.71, 403, 51.09, 48.98]

    return weights, densities, SDdensities, areas, lengths, volumes, diameters


def commonDateProperties(experiment):
    if experiment == 1:
        startTime = datetime(2016, 5, 1, 0, 0, 0)
        endTime = datetime(2016, 8, 1, 0, 0, 0) #8
    if experiment == 2:
        startTime = datetime(2016, 3, 1, 0, 0, 0)
        endTime = datetime(2016, 5, 15, 0, 0, 0) # 5
    if experiment == 3:
        startTime = datetime(2015, 11, 20, 0, 0, 0)
        endTime = datetime(2016, 4, 1, 0, 0, 0) # 4
    if experiment == 4:
        startTime = datetime(2016, 5, 1, 0, 0, 0)
        endTime = datetime(2016, 8, 1, 0, 0, 0) #8
    if experiment == 5:
        startTime = datetime(2015, 8, 1, 0, 0, 0)
        endTime = datetime(2016, 8, 1, 0, 0, 0)

    return startTime, endTime


def kelpProperties(num, kelpTypes):
    # Get the options of weights and densities
    weights, densities, SDdensities, areas, lengths, volumes, diameters = commonKelpProperties()

    # Loop over num release dates and randomly select type from provided kelpTypes list
    # Kelptypes list contents indicate which of the indices/options you have selected in commonKelpProperties
    # If onlye old blades: kelpTypes=[0]

    kelpWeights = []
    kelpDensities = []
    kelpAreas = []
    kelpLengths = []
    kelpVolumes = []
    kelpDiameters = []
    kelpType = []

    for i in xrange(len(num)):
        ind = random.choice(kelpTypes)
        kelpWeights.append(weights[ind])
        kelpType.append(ind)
        # Calculate a random density based on the mean and std values

        # Divide STD by 3 to get teh sigma value require dby the function
        stds = np.random.normal(densities[ind], SDdensities[ind] / 3.0)

        kelpDensities.append(stds)
        kelpAreas.append(areas[ind])
        kelpVolumes.append(volumes[ind])
        kelpDiameters.append(diameters[ind])
        kelpLengths.append(lengths[ind])

    return kelpWeights, kelpDensities, kelpAreas, kelpVolumes, kelpDiameters, kelpLengths, kelpType


def createOutputFilenames(experiment, polygonIndex, shapefile, verticalBehavior):
    startTime, endTime = commonDateProperties(experiment)
    startDate = ''
    if startTime.day < 10:
        startDate += '0%s' % (startTime.day)
    else:
        startDate += '%s' % (startTime.day)

    if startTime.month < 10:
        startDate += '0%s' % (startTime.month)
    else:
        startDate += '%s' % (startTime.month)

    startDate += '%s' % (startTime.year)

    endDate = ''
    if endTime.day < 10:
        endDate += '0%s' % (endTime.day)
    else:
        endDate += '%s' % (endTime.day)

    if endTime.month < 10:
        endDate += '0%s' % (endTime.month)
    else:
        endDate += '%s' % (endTime.month)

    endDate += '%s' % (endTime.year)

    # Special file naming for KINO. Each layer has name 'species.shp' and we want teh species name only.
    head, tail = os.path.split(shapefile)

    specie = "Kelp"
    outputFilename = 'results/%s_polygon_%s_experiment_%s_%s_to_%s.nc' % (
        specie, polygonIndex, experiment, startDate, endDate)
    animationFilename = 'figures/%s_polygon_%s_experiment_%s_%s_to_%s.mp4' % (
        specie, polygonIndex, experiment, startDate, endDate)
    plotFilename = 'figures/%s_polygon_%s_experiment_%s_%s_to_%s.png' % (
        specie, polygonIndex, experiment, startDate, endDate)

    if not os.path.exists('figures'):
        os.makedirs('figures')
    if not os.path.exists('results'):
        os.makedirs('results')
    return outputFilename, animationFilename, plotFilename


def createAndRunSimulation(use_svim, experiment, mapResolution, interMethod, lowDepth, highDepth, layer, polygonIndex,
                           shapefile, outputFilename, animationFilename, plotFilename, kinoDirectory, pattern_kino,
                           svimfiles2015, svimfiles2016, verticalBehavior, allNum,
                           allReleaseTimes, allKelpWeights, allKelpDensities, allKelpAreas, allKelpVolumes,
                           allKelpDiameters, allKelpLengths, allKelpTypes):
    # Setup a new simulation
    o = PelagicPlanktonDrift(loglevel=1)  # Set loglevel to 0 for debug information
    startTime, endTime = commonDateProperties(experiment)

    allKelpWeights_flat = [item for sublist in allKelpWeights for item in sublist]
    allKelpDensities_flat = [item for sublist in allKelpDensities for item in sublist]
    allReleaseTimes_flat = [item for sublist in allReleaseTimes for item in sublist]
    allKelpAreas_flat = [item for sublist in allKelpAreas for item in sublist]
    allKelpVolumes_flat = [item for sublist in allKelpVolumes for item in sublist]
    allKelpDiameters_flat = [item for sublist in allKelpDiameters for item in sublist]
    allKelpLengths_flat = [item for sublist in allKelpLengths for item in sublist]
    allKelpTypes_flat = [item for sublist in allKelpTypes for item in sublist]

    allNum_flat = [item for sublist in allNum for item in sublist]

    print "=> Simulation will release a total of %s particles\n" % (np.sum(allNum_flat))

    # Randomly distribute the particles at depths varying between lowDepth and highDepth
    depths = [randint(lowDepth, highDepth) for i in range(len(allNum_flat))]

    #######################
    # Preparing readers
    #######################
    reader_basemap = reader_basemap_landmask.Reader(
        llcrnrlon=15, llcrnrlat=68,
        urcrnrlon=23, urcrnrlat=74,
        resolution=mapResolution, projection='merc')
    o.add_reader([reader_basemap])  # Do not include basemap when stranding is deactivated

    reader_kino = reader_ROMS_native.Reader([s for s in pattern_kino])
    reader_kino.interpolation = interMethod

    if use_svim:
        reader_svim2015 = reader_ROMS_native.Reader(svimfiles2015)
        reader_svim2015.interpolation = interMethod
        reader_svim2016 = reader_ROMS_native.Reader(svimfiles2016)
        reader_svim2016.interpolation = interMethod

        o.add_reader([reader_kino, reader_svim2016, reader_svim2015])
    else:
        o.add_reader([reader_kino])

    #######################
    # Adjusting configuration
    #######################
    o.set_config('processes:turbulentmixing', True)
    o.set_config('turbulentmixing:diffusivitymodel', 'windspeed_Sundby1983')
    o.set_config('turbulentmixing:timestep', 1)  # secondsS
    o.set_config('turbulentmixing:verticalresolution', 1)  # default is 1 meter, but since we have longer timestep we justify it
    o.set_config('processes:verticaladvection', True)
    o.set_config('turbulentmixing:TSprofiles', False)
    o.set_config('turbulentmixing:max_iterations', 100)

    o.set_config('drift:scheme', 'euler')
    o.set_config('general:coastline_action', 'previous')  # Prevent stranding, jump back to previous position
    print o

    # depths=[randint(lowDepth,highDepth) for i in xrange(len(allKelpWeights))]
    allNum_flat = map(int, allNum_flat)
    for i, nums in enumerate(allNum_flat):

        if nums <= 0:
            continue
        print "Running i=%s num=%s and polygon=%s" % (i, nums, polygonIndex)

        o.seed_from_shapefile(shapefile, allNum_flat[i], featurenum=[polygonIndex],
                              z="seafloor+1",  # depths[i],
                              weight=allKelpWeights_flat[i],
                              density=allKelpDensities_flat[i],
                              area=allKelpAreas_flat[i],
                              volume=allKelpVolumes_flat[i],
                              diameter=allKelpDiameters_flat[i],
                              length=allKelpLengths_flat[i],
                              time=allReleaseTimes_flat[i],
                              plantpart=allKelpTypes_flat[i])

    # reader_basemap.plot()

    #########################
    # Run the model
    #########################
   # o.plot()

    o.run(end_time=endTime, time_step=timedelta(hours=2),
          outfile=outputFilename)
    # export_variables=['lon', 'lat', 'z','temp','length','weight','survival'])
    # print o


def setupSeed(seedCount, intervalHours, startTime, endTime, startReleaseTime, endReleaseTime, releaseParticles):
    ##################################################
    # Create seed variation as function of day
    # Called multiple times  from setupSeedsForExperiment
    ##################################################

    difference = endTime - startTime
    hoursOfSimulation = divmod(difference.total_seconds(), 3600)

    difference = endReleaseTime - startReleaseTime
    hoursOfRelease = divmod(difference.total_seconds(), 3600)

    timeStepsSimulation = int(int(hoursOfSimulation[0]) / 3)

    # print "=>Release: Simulated Release will run for %s simulation hours\n initiated on %s and ending on %s"%(timeStepsSimulation,startReleaseTime,endReleaseTime)

    interval = timedelta(hours=intervalHours)
    hoursPerRelease = divmod(interval.total_seconds(), 3600)  # hours per Release event
    timeStepsRelease = int(int(hoursOfRelease[0]) / int(hoursPerRelease[0]))  # number of Release timesteps
    ReleaseTimes = [startReleaseTime + interval * n for n in range(timeStepsRelease)]  # times of Release

    # num=np.random.normal(releaseParticles,int(releaseParticles/2)-1, size=len(ReleaseTimes)).astype(int)
    num = [releaseParticles for n in range(timeStepsRelease)]
    # num=np.sort(num) #sort particles in increasing order

    print "=> Seed episode: %s => Release of %s kelp particles" % (seedCount, np.sum(num))

    return num, ReleaseTimes


def setupSeedsForExperiment(experiment, releaseParticles):
    print "\nSeed setup started --------"
    seedCount = 1
    allNum = []
    allReleaseTimes = []
    allKelpProps = []
    allKelpWeights = []
    allKelpDensities = []
    allKelpAreas = []
    allKelpVolumes = []
    allKelpDiameters = []
    allKelpLengths = []
    allKelpTypes = []

    # Batch one : Old lamina (77%) released evenly 6 times a day between 122 and 135. 
    startTime, endTime = commonDateProperties(experiment)
    startReleaseTime = startTime
    endReleaseTime = endTime

    intervalHours = 6
    print "=> Release: daily: %s to %s" % (startReleaseTime, endReleaseTime)

    num, ReleaseTimes = setupSeed(seedCount, intervalHours, startTime, endTime, startReleaseTime, endReleaseTime,
                                  releaseParticles)
    allNum.append(num)
    allReleaseTimes.append(ReleaseTimes)
    seedCount += 1

    # Release only old blades
    kelpTypes = [0]
    kelpWeights, kelpDensities, kelpAreas, kelpVolumes, kelpDiameters, kelpLengths, kelpTypes = kelpProperties(num,
                                                                                                               kelpTypes)
    allKelpDensities.append(kelpDensities)
    allKelpWeights.append(kelpWeights)
    allKelpAreas.append(kelpAreas)
    allKelpVolumes.append(kelpVolumes)
    allKelpDiameters.append(kelpDiameters)
    allKelpLengths.append(kelpLengths)
    allKelpTypes.append(kelpTypes)
    # Batch two : New lamina (23%) released evenly 6 times a day once a week 
    # Find the number of weeks between start and stop
    num_of_weeks = int(math.ceil((endTime - startTime).days / 7.0))

    # Release once per week
    for week in xrange(int(num_of_weeks) - 1):
        dayN = week * 7
        startReleaseTime = startTime + timedelta(days=dayN)
        endReleaseTime = startTime + timedelta(days=dayN + 1)
        print "=> Release: weekly: %s to %s" % (startReleaseTime, endReleaseTime)
        intervalHours = 6

        num, ReleaseTimes = setupSeed(seedCount, intervalHours, startTime, endTime, startReleaseTime, endReleaseTime,
                                      releaseParticles)
        allNum.append(num)
        allReleaseTimes.append(ReleaseTimes)

        seedCount += 1
        # Define the properties of kelp to be released in batch 2
        # new blades, stipes, and fragments (in equal counts - check with Eli)
        kelpTypes = [1, 2, 3, 4]
        kelpWeights, kelpDensities, kelpAreas, kelpVolumes, kelpDiameters, kelpLengths, kelpTypes = kelpProperties(num,
                                                                                                                   kelpTypes)
        allKelpDensities.append(kelpDensities)
        allKelpWeights.append(kelpWeights)
        allKelpAreas.append(kelpAreas)
        allKelpVolumes.append(kelpVolumes)
        allKelpDiameters.append(kelpDiameters)
        allKelpLengths.append(kelpLengths)
        allKelpTypes.append(kelpTypes)

    print "Seed setup done --------\n"
    # Return the total number of particles per release date
    return allNum, allReleaseTimes, allKelpWeights, allKelpDensities, allKelpAreas, allKelpVolumes, allKelpDiameters, allKelpLengths, allKelpTypes


#########################
# SETUP FOR KELP PROJECT
#########################

runLocally = True
experiments = [5]
releaseParticles = 5

for experiment in experiments:
    lowDepth, highDepth = -4, -2  # in negative meters
    verticalBehavior = False
    startTime, endTime = commonDateProperties(experiment)

    kinoDirectory = '/imr/vol1/NorFjords5/Malangen-160m_AUG2015-AUG2016/'
    if runLocally:
        kinoDirectory = "/Volumes/home/CloudStation/NorFjord/"

    #if (startTime.year<endTime.year):
    apattern = 'norfjords_160m_avg.nc4_*'
    #else:
    #    apattern = 'norfjords_160m_avg.nc4_%s*' % (startTime.year)

    argument = "%s%s" % (kinoDirectory, apattern)

    pattern_kino = glob.glob(argument)
    pattern_kino.sort()

    svim2015 = '/Volumes/home/CloudStation/SVIM/2015/ocean_avg_2015%02d.nc4'
    svimfiles2015 = [svim2015 % d for d in range(1, 13, 1)]

    svim2016 = '/Volumes/home/CloudStation/SVIM/2016/ocean_avg_2016%02d.nc4'
    svimfiles2016 = [svim2016 % d for d in range(1, 13, 1)]

    use_svim = True
    interMethod = 'linearNDFast'  # linearND

    mapResolution = "h"  # f, c, i, h

    print "============================================="
    print "Starting kelp simulations for experiment %s\n" % (experiment)
    print "    %s" % (time.strftime("%c"))
    print "============================================="
    allNum, allReleaseTimes, allKelpWeights, allKelpDensities, allKelpAreas, allKelpVolumes, allKelpDiameters, allKelpLengths, allKelpTypes = setupSeedsForExperiment(
        experiment, releaseParticles)

    shapefile = "/work/shared/nn9297k/Kelp/Shapefile/kelpexpol_exp_grazed_combined.shp"
    if runLocally:
        shapefile = '/Users/trondkr/Dropbox/NIVA/KelpFloat/Kelp/Shapefile/ShapefilesHGU/kelpexpol_exp_grazed_combined.shp'  # kelpexpol_exp_grazed.shp'
    # Remember to setup/edit kelp properties in function: commonKelpProperties unless
    # using the default values

    # Setup complete --------------------------------------------------------------
    # No edit below

    print "=> Using shapefile %s" % shapefile
    s = ogr.Open(shapefile)

    # Find all kelp polygons in Shapefile
    for layer in s:
        polygons = [x + 1 for x in xrange(layer.GetFeatureCount() - 1)]
        polygons=[13]

        print('Running for layer with %s features)' % (layer.GetFeatureCount()))
        # Loop over all kelp polygons, releasing kelp and tracking their drift, writing results to file
        for polygonIndex in polygons:

            feature = layer.GetFeature(polygonIndex - 1)

            print "Kelp area %s for polygon %s" % (feature.GetGeometryRef().GetArea(), polygonIndex)
            geom = feature.GetGeometryRef()
            points = geom.GetGeometryCount()
            ring = geom.GetGeometryRef(0)

            if ring.GetPointCount() > 3:
                outputFilename, animationFilename, plotFilename = createOutputFilenames(experiment, polygonIndex,
                                                                                        shapefile, verticalBehavior)

                print "Result files will be stored as:\nnetCDF=> %s\nmp4=> %s" % (outputFilename, animationFilename)

                print "Starting simulations...."
                createAndRunSimulation(use_svim, experiment, mapResolution, interMethod, lowDepth, highDepth,
                                       layer, polygonIndex, shapefile,
                                       outputFilename, animationFilename, plotFilename,
                                       kinoDirectory, pattern_kino, svimfiles2015, svimfiles2016,
                                       verticalBehavior, allNum, allReleaseTimes, allKelpWeights,
                                       allKelpDensities, allKelpAreas, allKelpVolumes, allKelpDiameters, allKelpLengths,
                                       allKelpTypes)
