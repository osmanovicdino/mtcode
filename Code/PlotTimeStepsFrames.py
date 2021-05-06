dfList = []

for i in range(0,100):
    if (i < 10):
        csvPath = 'C:\cygwin64\home\passa\mtcode\Code\pos_i=00' + str(i) + '.csv'
    else:
        csvPath = 'C:\cygwin64\home\passa\mtcode\Code\pos_i=0' + str(i) + '.csv'
    dataFrame = pd.read_csv(csvPath, names = ['x','y','z'])
    stepList = (np.zeros(13) + i*np.ones(13)).tolist()
    dataFrame['step_id'] = stepList
    # add "labels" to separate nanostars
    particlesPerNanostar = 13
    iter = 0
    labelList = []
    for row in dataFrame:
        if (iter < particlesPerNanostar - 1):
            labelList.append('Nanostar A')
        else:
            labelList.append('Nanostar B')
    dfList.append(dataFrame)

# convert frames to figures
# figures to png files
for d in range(0, 100):
    tempFig = px.scatter_3d(dfList[d], x='x', y='y',z='z', range_x=[0,11], range_y=[0,11], range_z=[0,11])
    tempFig.update_layout(
        scene = {
            'aspectmode': 'cube'
        }
    )
    savePath = 't='+str(d)+'.png' #C:\cygwin64\home\passa\mtcode\images\test\t=' + str(d) + '.png'
    tempFig.write_image(savePath)
