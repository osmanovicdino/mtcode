#!/usr/bin/env python3

# import matplotlib.pyplot as plt
# import numpy as np
import csv
import sys
# from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
# import os
from PlotFrameFunction import graph



cmdargs = list(map(str,sys.argv))

filename =  cmdargs[1]


colarraytemp='red'

data=[]

if(len(cmdargs)> 2):
    x= False
    colarraystring = cmdargs[2]
    with open(colarraystring, newline='') as csvfile:
        data = list(csv.reader(csvfile,quoting=csv.QUOTE_NONNUMERIC))
else:
    x= True
    data.append(colarraytemp)

print(data[0])

graph(filename,data)

# outputfilename = os.path.splitext(filename)[0]+'.jpg'

# fig = plt.figure()

# ax = fig.add_subplot(111,projection='3d')

# def graph():
#     with open(filename, newline='') as csvfile:
#         spamreader = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
#         for eachLine in spamreader:
#             xpos = eachLine[0]
#             ypos = eachLine[1]
#             zpos = eachLine[2]

#             ax.scatter(xpos,ypos,zpos,color='blue')

#             ax.set_xlabel('x')
#             ax.set_ylabel('y')
#             ax.set_zlabel('z')

#     plt.savefig(outputfilename, format='jpg')
    #plt.show()


