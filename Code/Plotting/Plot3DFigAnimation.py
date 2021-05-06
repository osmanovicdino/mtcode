import plotly.express as px
import pandas as pd
import numpy as np
import plotly.graph_objects as go

dfList = []

for i in range(0,100):
    if (i < 10):
        csvPath = 'C:\cygwin64\home\passa\mtcode\Code\out_folder\cond6\pos_i=00' + str(i) + '.csv'
    else:
        csvPath = 'C:\cygwin64\home\passa\mtcode\Code\out_folder\cond6\pos_i=0' + str(i) + '.csv'
    dataFrame = pd.read_csv(csvPath, names = ['x','y','z'])
    stepList = (np.zeros(26) + i*np.ones(26)).tolist()
    dataFrame['step_id'] = stepList # tagging the step number for the frame 


    # need to color the stars differently
    # each point on the star should have name (every 13 lines)
    nameList = []

    for d in range(0,13):
        nameList.append('Nanostar 1')

    for d in range(0,13):
        nameList.append('Nanostar 2')

    dataFrame['name'] = nameList
    dfList.append(dataFrame)

# for d in range(0, 100):
#     tempFig = px.scatter_3d(dfList[d], x='x', y='y',z='z',color='name', range_x=[0,11], range_y=[0,11], range_z=[0,11])
#     tempFig.update_layout(
#         scene = {
#             'aspectmode': 'cube'
#         }
#     )
#     savePath = 't='+str(d)+'.png' #C:\cygwin64\home\passa\mtcode\images\test\t=' + str(d) + '.png'
#     tempFig.write_image(savePath)


# # Create the frames
# frames = []
# imgs = glob.glob("*.png")
# for i in imgs:
#     new_frame = Image.open(i)
#     frames.append(new_frame)

# # Save into a GIF file that loops forever
# frames[0].save('cond1_animation.gif', format='GIF',
#                append_images=frames[1:],
#                save_all=True,
#                duration=125, loop=0)


all = pd.concat(dfList)
figureToRender = px.scatter_3d(all, x='x', y='y',z='z',animation_frame='step_id', color='name',range_x=[0,11], range_y=[0,11], range_z=[0,11])
figureToRender.update_layout(
    scene = {
        'aspectmode': 'cube'
    }
)

figureToRender.update_traces(marker=dict(size=6,
                              ),
                  selector=dict(mode='markers'))

figureToRender.show()
