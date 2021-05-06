import plotly.express as px
import pandas as pd

# read in the data as a dataframe
positionData = pd.read_csv('C:\cygwin64\home\passa\mtcode\Code\pos_i=099.csv', names=['x', 'y', 'z'])
figureToRender = px.scatter_3d(positionData, x='x', y='y',z='z',size_max=2)
figureToRender.show()
