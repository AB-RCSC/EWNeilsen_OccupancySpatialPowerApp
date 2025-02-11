

import json
import numpy as np
import statsmodels as sm
import math
from osgeo import gdal, ogr, gdalconst
import numpy as np
import os
import sys

import random
import inspect
import pandas as pd
import shutil

import matplotlib.pyplot as plt

rdriver = gdal.GetDriverByName('GTiff')
vdriver = ogr.GetDriverByName('ESRI Shapefile')
    
def normScaler (inArray):
    outArray = (inArray - np.min(inArray)) / (np.max(inArray) - np.min(inArray))
    return outArray

def askQuestion(respType,question):

    keepAsking = True
    while keepAsking:
        
        resp = input(question)
        
        if (respType == "float"):            
            try:
                resp = float(resp)                 
            except:
                print("Enter a number between 0 and 1.")
            if (0 <= resp <= 1):
                keepAsking = False 
            else:
                print("Enter a number between 0 and 1.")
                
        if (respType == "int"):            
            try:
                resp = int(resp)
                keepAsking = False
            except:
                print("Enter a whole number.")
        
        if (respType == "string"):            
            try:
                resp = str(resp)
                keepAsking = False
            except:
                print("Enter text.")
                
    return resp
    
def inputSpatial():

    global cwd
    cwd = os.getcwd() 
    print(cwd)    
    
    ## raster info ##
    global nodata;  global extent; global srs; global n_rows; global n_cols; global cell_size;
    nodata = -9999
    cell_size = 5000
    SAshp = vdriver.Open(cwd + '/Data/SpatialLayers/Test_SA.shp')
    rlayer = SAshp.GetLayer()
    extent = rlayer.GetExtent()
    srs = rlayer.GetSpatialRef()
    n_rows = int(math.ceil(abs(extent[1] - extent[0]) / cell_size))
    n_cols = int(math.ceil(abs(extent[3] - extent[2]) / cell_size))
    print("The study area has " + str(extent)  + ". \nIt has " + str(n_rows) + " rows and "  + str(n_cols) + " columns." );print('')
    
    global usePredictors
    usePredictors = {}
    
    ## elevation
    rasIn = gdal.Open(cwd + '/Data/SpatialLayers/DEM_5000_m.tif', gdalconst.GA_ReadOnly)
    band1 = rasIn.GetRasterBand(1)
    band1 = band1.ReadAsArray()
    elevationArray = normScaler(band1)
    usePredictors["elevation"] = elevationArray
    
    ## water proximity
    rasIn = gdal.Open(cwd + '/Data/SpatialLayers/WaterProx_dist_5000_m.tif', gdalconst.GA_ReadOnly)
    band1 = rasIn.GetRasterBand(1)
    band1 = band1.ReadAsArray()
    waterProxArray = normScaler(band1)
    usePredictors["water_proximity"] = waterProxArray
    
    print('------------------')   
    print("Finished inputing Use vars:", list(usePredictors.keys()));print('')

def simulateReponse():
   
    global usePredictors
    global nodata; global cell_size; global extent; global srs; global n_rows; global n_cols    
    global cwd
    
    global responseValues
    responseValues = {}  
   
    ###########################
    ## PROBABILITY OF USE AS PRODUCT OF INPUT LAYERS (for now)    
    varArray = []
    for i in usePredictors:
        varArray.append(usePredictors[i])      
    varArray = np.array(varArray)    
    responseValues["Use"] = normScaler(np.sum(varArray, axis=0))
    print('------------------')
    print("Used the inputted rasters to simulate the spatial probability of use across study area.");print('')

    ###########################
    ## CONVERT TO OCCUPANCY
    """True occupancy translates into a
    proportion of the landscape occupied.
    So, the prob of use, which is a function
    of covariates, needs to be discretized
    using trueOcc. So, if trueOcc is 0.2, the
    top twenty percent of prob use values are occupied."""
    
    global trueOcc
    print('------------------')
    trueOcc = askQuestion("float","Converting use into occupancy. What is the true proportion of the area that is occupied (number between 0 and 1)?")
    occThreshold = np.quantile(responseValues["Use"], 1 - trueOcc)        
    responseValues["Occupancy"] = np.zeros(shape=responseValues["Use"].shape)
    responseValues["Occupancy"][responseValues["Use"] > occThreshold] = 1     
    pxN = sum( [1 for line in responseValues["Occupancy"] for x in line if x ==1 ] )
    cellArea = cell_size ** 2
    saAreaKM = (cellArea/1000000) * pxN
    print("There are " + str(pxN) + " occupied pixels (" + str(saAreaKM) + " km occupied area). This leads to an instantaneous probability of detection in any cell for one, randomly moving, individual of " + str(round(1/pxN,4)));print('')

    ###########################
    ## SPATIAL PROBABILITY OF USE FOR A POPULATION 
    """ Density will apply to a number of
    animals over the area of the occupied cells """   
    print('------------------')
    dens = askQuestion("float","Simulate a population within the occupied cells using a population density. What is the density of individuals per km2 (0.001 - 1)?")
    global N
    N = float(dens) * saAreaKM
    global popPX
    popPX = N/pxN   
    if popPX > 1:
        popPX = 1
    print('') 
    print("With a density of " + str(dens) + " individuals per pixel across all occupied pixels, the total population is " + str(round(N,2)) + ". This gives an instantaneous probability of use of any occupied cell, of any randomly moving individual, of " + str(round(popPX,4)))
    
    ###########################
    ## SPATIAL PROBABILITY OF DETECTION    
    """ Use the product of occupancy and use to 
    extract the probability of use at occupied 
    sites only - we are currently assuming perfect
    detection at cameras. """
    detArray = []
    for i in responseValues:
        detArray.append(responseValues[i])    
    detArray = np.array(detArray)
    spatDet = normScaler(np.prod(detArray, axis=0))   
    spatDet = spatDet * popPX
    responseValues["Detection"] = spatDet  
    # print(np.mean(spatDet))
    global meanDetection
    meanDetection = np.mean(spatDet[spatDet != 0])
    print("The mean instantaneous probability of detection across occupied cells, for any randomly moveing individual, is " + str(round(meanDetection,4)));print('')
        
    ###########################
    ## RASTER OUTPUTS
    for res in responseValues:
        responsePath = cwd + "/Data/SpatialLayers/Simulated" + res + r'.tif'
        rasterizedDS = rdriver.Create(responsePath, n_rows, n_cols, 1, gdal.GetDataTypeByName("Float32"))
        rasterizedDS.SetGeoTransform([extent[0], cell_size, 0, extent[3], 0, -1 * cell_size])    
        rasterizedDS.SetProjection(srs.ExportToWkt());rasterizedDS.GetRasterBand(1).SetNoDataValue(nodata)
        rasterizedDS.GetRasterBand(1).Fill(nodata)
        rasterizedDS.GetRasterBand(1).WriteArray(responseValues[res])        
        rasterizedDS = None
        responsePath = None
    
    ###########################
    ##  PLOTTING
    print("Finished calculating distribution vars:", list(responseValues.keys()))
    for i in responseValues:
        plt.close()
        pltDat = responseValues[i]
        plt.title(i + str(np.amin(pltDat)) + "_" + str(np.amax(pltDat)))
        plt.imshow(pltDat)
        plt.show() 
  
def simulateOccupancyData():
    
    global nodata; global cell_size; global extent; global srs; global n_rows; global n_cols
    global responseValues
    global prevPos
    global responseVar
    global cwd
    
    ###########################
    ## SAVE DETECTION HISTORIES LOCATION
    # if os.path.isdir(cwd + '/Data/DetectionHistories'):
        # os.remove(cwd + '/Data/DetectionHistories')
    os.makedirs(cwd + '/Data/DetectionHistories', exist_ok=True)
    
    ###########################
    ## DEPLOYMENT SCENARIOS
    ## type of camera deployment
    camConfig = askQuestion("int","Enter the configuration of cameras \n (1) Systematic\n (2) Random\n (3) Stratigied Random \n")    
    
    ## deployment scenarios variables
    """The tool can assess multiple effort levels,
    defined here as the combination of number of
    sites and lenght of the deployment"""
    print('------------------')
    siteScenN = askQuestion("int","Enter the number of site scenarios.")
    maxCam = askQuestion("int","Enter the max number of cameras.")
    minCam = askQuestion("int","Enter the min number of cameras.")      
    durScenN = askQuestion("int","Enter the number of duration scenarios.")
    maxDur = askQuestion("int","Enter the max duration of deployments (weeks).")
    minDur = askQuestion("int","Enter the min duration of deployments (weeks).")
    
    sitesN = range(minCam,maxCam,round(maxCam/siteScenN)) 
    #sitesN = [40]
    dursN = range(minDur,maxDur,round(maxDur/durScenN))    
  
    ###########################
    ## SIMLUATE DETECTION HISTORIES
    for sn in sitesN:    
        sn = sn + 1        
        for dn in dursN:  

            ###########################
            ## ASSIGNING SITES
            
            ##SYSTEMATIC SITES
            if (camConfig == 1):            
                rArray = responseValues["Use"]
                rSize = rArray.size
                rSysN = int(round(rSize / sn))
                ##print("every " + str(rSysN))        
                #print("Setting pixels as sites.")    
                siteList = np.zeros(shape=responseValues["Use"].shape)
                siteList = np.ndarray.flatten(siteList)
                siteList[1::rSysN] = 1
                siteList = np.reshape(siteList,(responseValues["Use"].shape))
            
            ## RANDOM SITES
            elif(camConfig == 2):           
                rArray = responseValues["Use"]
                rSize = rArray.size
                siteInds = []
                for i in range(0, sn):
                     siteInds.append(random.randrange(rSize))
                siteList = np.zeros(shape=responseValues["Use"].shape)
                siteList = np.ndarray.flatten(siteList)
                siteList[siteInds] = 1
                siteList = np.reshape(siteList,(responseValues["Use"].shape))
                
            # ##  PLOTTING
            # plt.close()
            # plt.title(str(sn) + "cameras")
            # plt.imshow(siteList)
            # plt.show()     
            
            #########################
            ## GET TRUE DETECTION AND OCCUPANCY VALUES AT SITES
            siteOcc = []
            siteDet = []
            siteIndex = 0
            for rInd, row in enumerate(siteList):
                for cInd, value in enumerate(row):            
                    if value == 1:
                        OccValue = responseValues["Occupancy"][rInd,cInd]
                        DetValue = responseValues["Detection"][rInd,cInd]
                        siteOcc.append(OccValue)
                        siteDet.append(DetValue)        
            siteOcc = np.array(siteOcc)
            siteDet = np.array(siteDet)       
            
            #########################   
            ## POPULATED DECTECTION HISTORIES (NO MISSING DATA FUNCTIONS YET)
           
            ScenName = str(sn) + "_" + str(dn)
            os.makedirs(cwd + '/Data/DetectionHistories/' + ScenName, exist_ok=True)    

            simN = 10       
            for simn in range(simN):        
                DetectionHistory = []
                
                for sd in siteDet:    
                    siteDetHist = []
                    
                    for sv in range(dn):        
                        det = 0
                        samp = random.uniform(0, 1)
                        if samp < sd:
                            det = 1            
                        siteDetHist.append(det)    
                        
                    DetectionHistory.append(siteDetHist)

                ## export 
                dh = pd.DataFrame(DetectionHistory)
                dh.to_csv(cwd + '/Data/DetectionHistories/' + ScenName + "/" + str(simn + 1 ) + '_dh.csv', index=False)   


                