import plotly.express as px
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import math 

distanceList = []
l = 12; 
for i in range(0,100):
    if (i < 10):
        csvPath = 'C:\cygwin64\home\passa\mtcode\Code\out_folder\cond3\pos_i=00' + str(i) + '.csv'
    else:
        csvPath = 'C:\cygwin64\home\passa\mtcode\Code\out_folder\cond3\pos_i=0' + str(i) + '.csv'
    dataFrame = pd.read_csv(csvPath, names = ['x','y','z'])
    # grab the center of the first nanostar 
    centerNanostar1 = [dataFrame['x'][0], dataFrame['y'][0], dataFrame['z'][0]]
    centerNanostar2 = [dataFrame['x'][12], dataFrame['y'][12], dataFrame['z'][12]] # this index change depending on particles count
    distance = 0
    # should generalize to grab all the centers
    centers = [centerNanostar1, centerNanostar2]
    # iterate through shells 
    # start from (0,0)
    

    drList = np.linspace(1, 100, 50).tolist()
    countValues = [0]*len(drList)

    for dr in drList: 
        for center in centers: 
            distance = (center[0])**2 + (center[1])**2 + (center[2])**2
            if distance < dr: 
                countValues[drList.index(dr)] += 1 

    # normalize by the volume 
    normValues = []
    for dr in drList: 
        normValues.append(countValues[drList.index(dr)] / (4/3*math.pi*dr**3))

    dataDict = {'distance': dr, 'nv': normValues}
    outputDataframe = pd.DataFrame(dataDict)
    fileName2Write = str(i) + '_gr_data.csv'
    outputDataframe.to_csv(fileName2Write)