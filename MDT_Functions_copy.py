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
import ipywidgets as widgets # added by alif
from IPython.display import display # added by alif
import matplotlib.pyplot as plt


rdriver = gdal.GetDriverByName('GTiff')
vdriver = ogr.GetDriverByName('ESRI Shapefile')

def normScaler (inArray):
    outArray = (inArray - np.min(inArray)) / (np.max(inArray) - np.min(inArray))
    return outArray
    
#######################################################################################################################
# # For trueOcc value--------------------------------------------------------
# # Error message
# trueOcc_error_label = widgets.Label(value="❌ Error: Please enter a number between 0 and 1 then press submit again", style={'color': 'red'})

# # trueOcc widget
# trueOcc_resp_widget = widgets.FloatText(
#             description= "True Occupncy: ",
#             min=0.0, max=1.0, step=0.01, value=0.0  # Default value set to 0.5
#         )

# trueOcc_button = widgets.Button(description="Submit")

# # Create an Output widget to display the message
# trueOcc_output = widgets.Output()

# # For dense values---------------------------------------------------------
# # Error message
# dens_error_label = widgets.Label(value="❌ Error: Please enter a number between 0.001 and 1 then press submit again", style={'color': 'red'})

# # integer widget
# dens_resp_widget = widgets.FloatText(
#             description= "Dense; ",
#             min=0.0, max=1.0, step=0.01, value=0.0  # Default value set to 0.5
#         )

# # int Button widget
# dens_button = widgets.Button(description="Submit")

# # Create an Output widget to display the message
# dens_output = widgets.Output()

# # SPATIAL PROBABILITY OF DETECTION----------------------------------------------
# spatial_prob_button = widgets.Button(description="calculate")

# # Create an Output widget to display the message
# spatial_output = widgets.Output()

#######################################################################################################################

# Create two integer sliders A and B
slider_A = widgets.FloatSlider(min=0.01, max=1, step=0.1, value=0.5, description="trueOcc:")
slider_B = widgets.FloatSlider(min=0, max=0.1, step=0.001, value=0.05, description="dense:")

# Output widgets for displaying values and the bar chart
value_output = widgets.Output()
plot_output = widgets.Output()
out = widgets.Output()

# Function to update the bar chart dynamically based on equation 2A + 3B
def update_bar_chart(A, B):
    with plot_output:
        plot_output.clear_output(wait=True)  # Clear previous plot
        
        global usePredictors
        global nodata; global cell_size; global extent; global srs; global n_rows; global n_cols    
        global cwd
        
        global responseValues
        responseValues = {}  
        inputSpatial()
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
        trueOcc = A
        occThreshold = np.quantile(responseValues["Use"], 1 - trueOcc)        
        responseValues["Occupancy"] = np.zeros(shape=responseValues["Use"].shape)
        responseValues["Occupancy"][responseValues["Use"] > occThreshold] = 1 
        global pxN, popPX, saAreaKM, N, meanDetection
        pxN = sum( [1 for line in responseValues["Occupancy"] for x in line if x ==1 ] )
        cellArea = cell_size ** 2
        saAreaKM = (cellArea/1000000) * pxN
        print("There are " + str(pxN) + " occupied pixels (" + str(saAreaKM) + " km occupied area). This leads to an instantaneous probability of detection in any cell for one, randomly moving, individual of " + str(round(1/pxN,4)));print('')
    
        ###########################
        ## SPATIAL PROBABILITY OF USE FOR A POPULATION 
        """ Density will apply to a number of
        animals over the area of the occupied cells """   
        print('------------------')
        dens = B
    
        N = float(dens) * saAreaKM
        popPX = N/pxN   
        if popPX > 1:
            popPX = 1
        print('') 
        print("With a density of " + str(dens) + " individuals per pixel across all occupied pixels, the total population is " + str(round(N,2)) + ". This gives an instantaneous probability of use of any occupied cell, of any randomly moving individual, of " + str(round(popPX,4)))
        
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
        global meanDetection
        meanDetection = np.mean(spatDet[spatDet != 0])
        print("The mean instantaneous probability of detection across occupied cells, for any randomly moveing individual, is " +    str(round(meanDetection,4)));print('')

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
        ##  PLOTTING
        print("Finished calculating distribution vars:", list(responseValues.keys()))
        # Create an Output widget 
        # out = widgets.Output()

        # Number of images
        num_images = len(responseValues)
        
        # Create subplots with the required number of columns
        fig, axes = plt.subplots(1, num_images, figsize=(num_images * 4, 4))  # Dynamic width

        with out:
            out.clear_output(wait=True)  # Clear previous plot
            for ax, (i, pltDat) in zip(axes, responseValues.items()):
                # print(pltDat)
                ax.set_title(f"{i} {np.amin(pltDat):.2f}_{np.amax(pltDat):.2f}")
                ax.imshow(pltDat, cmap='viridis')  # Adjust colormap as needed
                ax.axis("off")  # Remove axes for a cleaner view
            
            plt.tight_layout()  # Adjust layout to prevent overlap
            plt.show()
                
        # Display the widget
        # display(out)  


# Function to print slider values dynamically
def update_values(change):
    with value_output:
        value_output.clear_output(wait=True)
        A, B = slider_A.value, slider_B.value
        # print(f"Slider A: {A}, Slider B: {B}, 2A + 3B: {2*A + 3*B}")
    update_bar_chart(A, B)  # Update chart when sliders change

# Attach the function to both sliders
slider_A.observe(update_values, names='value')
slider_B.observe(update_values, names='value')



def askQuestion(respType,question):

    keepAsking = True
    # while keepAsking:
    #     string_resp_widget = widgets.Text(
    #         description= question
    #     )

    #     if (respType == "trueOcc"):
    #         display(trueOcc_hbox)
    #         keepAsking = False
                
    #     if (respType == "dense"):            
    #         display(dens_hbox)
    #         keepAsking = False
        
    #     if (respType == "string"):   
    #         display(string_resp_widget)
    #         if string_resp_widget.value:
    #             keepAsking = False
    #             return string_resp_widget.value

                
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

    # Display widgets
    slider_Vbox = widgets.VBox([slider_A, slider_B], layout=widgets.Layout(width='30%'))
    graph_Vbox = widgets.VBox([plot_output, out], layout=widgets.Layout(width='70%'))
    main_Hbox = widgets.HBox([slider_Vbox, graph_Vbox])
    display(main_Hbox)
    # display(slider_A, slider_B, value_output, plot_output)
    
    # Initialize chart on first load
    update_values(None)
    # trueOcc = askQuestion("trueOcc","Converting use into occupancy. What is the true proportion of the area that is occupied (number between 0 and 1)?")
    # dens = askQuestion("dense","Simulate a population within the occupied cells using a population density. What is the density of individuals per km2 (0.001 - 1)?")
    # spatial_prob_button.on_click(lambda b: spatial_on_button_click(responseValues, b))






def func_simulateOccupancyData(camConfig, siteScenN, maxCam, minCam, durScenN, maxDur, minDur):
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






camConfig = widgets.IntSlider(min=1, max=3, step=1, value=0, description="camConfig:")
siteScenN = widgets.IntSlider(min=1, max=50, step=1, value=0, description="siteScenN:")
maxCam = widgets.IntSlider(min=1, max=200, step=1, value=0, description="maxCam:")
minCam = widgets.IntSlider(min=1, max=200, step=1, value=0, description="minCam:")
durScenN = widgets.IntSlider(min=1, max=50, step=1, value=0, description="durScenN:")
maxDur = widgets.IntSlider(min=1, max=50, step=1, value=0, description="maxDur:")
minDur = widgets.IntSlider(min=1, max=50, step=1, value=0, description="minDur:")








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
    
    # ###########################
    # ## DEPLOYMENT SCENARIOS
    # ## type of camera deployment
    # camConfig = askQuestion("int","Enter the configuration of cameras \n (1) Systematic\n (2) Random\n (3) Stratigied Random \n")    
    
    # ## deployment scenarios variables
    # """The tool can assess multiple effort levels,
    # defined here as the combination of number of
    # sites and lenght of the deployment"""
    # print('------------------')
    # siteScenN = askQuestion("int","Enter the number of site scenarios.")
    # maxCam = askQuestion("int","Enter the max number of cameras.")
    # minCam = askQuestion("int","Enter the min number of cameras.")      
    # durScenN = askQuestion("int","Enter the number of duration scenarios.")
    # maxDur = askQuestion("int","Enter the max duration of deployments (weeks).")
    # minDur = askQuestion("int","Enter the min duration of deployments (weeks).")
    simulateOccupancyData_Vbox = widgets.VBox([camConfig, siteScenN, maxCam, minCam, durScenN, maxDur, minDur], layout=widgets.Layout(width='30%'))
    display(simulateOccupancyData_Vbox)



