import plotly.express as px
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import math 

distanceList = []
stepList = []

for i in range(0,100):
    if (i < 10):
        csvPath = 'C:\cygwin64\home\passa\mtcode\Code\out_folder\cond6\pos_i=00' + str(i) + '.csv'
    else:
        csvPath = 'C:\cygwin64\home\passa\mtcode\Code\out_folder\cond6\pos_i=0' + str(i) + '.csv'
    dataFrame = pd.read_csv(csvPath, names = ['x','y','z'])
    # grab the center of the first nanostar 
    centerNanostar1 = [dataFrame['x'][0], dataFrame['y'][0], dataFrame['z'][0]]
    centerNanostar2 = [dataFrame['x'][13], dataFrame['y'][13], dataFrame['z'][13]] # this index change depending on particles count
    distance = (centerNanostar2[0] - centerNanostar1[0])**2 + (centerNanostar2[1] - centerNanostar1[1])**2 + (centerNanostar2[2] - centerNanostar1[2])**2
    # print(str(distance))
    distanceList.append(math.sqrt(distance))
    stepList.append(i)

# make dictionary 
tempDict = {'time_step': stepList, 'distance': distanceList}
distanceDf = pd.DataFrame(tempDict) 
fig = px.line(distanceDf, x="time_step", y="distance")
fig.update_xaxes(range=[0, 100])
fig.update_yaxes(range=[0, 15])
fig.show()
