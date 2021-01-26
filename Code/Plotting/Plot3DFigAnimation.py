import plotly.express as px
import pandas as pd
import numpy as np
import plotly.graph_objects as go

dfList = []

for i in range(0,100):
    if (i < 10):
        csvPath = 'C:\cygwin64\home\passa\mtcode\Code\pos_i=00' + str(i) + '.csv'
    else:
        csvPath = 'C:\cygwin64\home\passa\mtcode\Code\pos_i=0' + str(i) + '.csv'
    dataFrame = pd.read_csv(csvPath, names = ['x','y','z'])
    stepList = (np.zeros(13) + i*np.ones(13)).tolist()
    dataFrame['step_id'] = stepList
    dfList.append(dataFrame)

all = pd.concat(dfList)
figureToRender = px.scatter_3d(all, x='x', y='y',z='z',animation_frame='step_id', range_x=[0,11], range_y=[0,11], range_z=[0,11])
figureToRender.update_layout(
    scene = {
        'aspectmode': 'cube'
    }
)

figureToRender.show()
