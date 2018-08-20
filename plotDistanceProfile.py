import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import sys
from swarmPlot import swarm_plot

def main():

    fig_width = 12
    fig_height = 5 
    plt.figure(figsize = (fig_width, fig_height))
 
    panel_width = .4
    panel_height = .8

    distsPanel = plt.axes([0.05, 0.1, panel_width, panel_height])
    timesPanel = plt.axes([0.55, 0.1, panel_width, panel_height])

    inFile = open("distanceProfile", "r")
    timeIndex = int(inFile.readline().split()[-1])
    memIndex = inFile.readline().split()[-1]
    inFile.readline()
   
    myDists = []
    oldDists = []
    myTimes = []
    oldTimes = []

    for l in inFile.readlines():

        line = l.split()
        myDists.append(int(line[0]))
        oldDists.append(int(line[1]))
        myTimes.append(int(line[2]))
        oldTimes.append(int(line[3]))
  
    inFile.close()

    ########### Scatter plot panel

    distsPanel.plot(oldDists, myDists, marker = 'o', linewidth = 0, 
                    markerfacecolor = 'black', 
                    markeredgewidth = 0, markersize = 2) 
    maxVal = max(max(myDists), max(oldDists))
    distsPanel.plot([0,maxVal],[0, maxVal], linestyle="-", color = "red", linewidth = 0.5)
   
    distsPanel.set_xlabel("Paths based distance")
    distsPanel.set_ylabel("Snarl based distance")
    distsPanel.set_title("Predicted Distances")
    d = sorted(set(oldDists))[-2]
    distsPanel.set_xlim(min(oldDists) * .95,  d * 1.05)
    distsPanel.set_ylim(min(myDists) * .95,  max(myDists) * 1.05)



    ############# Time plot panel

    ''' 
    #Plot time by distances 
    timesPanel.plot(myDists, myTimes, marker = 'o', linewidth = 0, 
                    markerfacecolor = 'red', 
                    markeredgewidth = 0, markersize = 2, label = "Snarl based") 
    timesPanel.plot(oldDists, oldTimes, marker = 'o', linewidth = 0, 
                    markerfacecolor = 'blue', 
                    markeredgewidth = 0, markersize = 2, label = "Path based") 
   
    minTime = 0
    maxCount = max(max(oldTimes), max(myTimes))
    maxTime = max(max(oldDists), max(myDists))
    '''
    
    totalTimeOld = sum(oldTimes)
    totalTimeNew = sum(myTimes) + timeIndex

    timeDict = {} #maps time to (# oldTimes, # myTimes)
    for time in oldTimes:
        if time not in timeDict:
            timeDict[time] = (1, 0)
        else:
            (a, b) = timeDict[time]
            timeDict[time] = (a+1, b)
    for time in myTimes:
        if time not in timeDict:
            timeDict[time] = (0, 1)
        else :
            (a, b) = timeDict[time]
            timeDict[time] = (a, b + 1)

    minTime = min(timeDict.keys())
    maxTime = max(timeDict.keys())
 
 
    maxCount = 0
    for i in range(maxTime):
       xStart = i * 3
       (h1, h2) = timeDict.get(i+1, (0, 0))
       maxCount = max(maxCount, h1, h2)
       rect1 = mplpatches.Rectangle([xStart, 0], 1, h1, \
                     facecolor = "red")
       rect2 = mplpatches.Rectangle([xStart+1, 0], 1, h2, \
                     facecolor = "blue")

       timesPanel.add_patch(rect1)
       timesPanel.add_patch(rect2)
 
    extra = mplpatches.Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
    string = "Time for calculating \n index: " + str(timeIndex) + " " 
    timesPanel.legend([rect1, rect2],
                      ("Path based\n total time: " + str(totalTimeOld), 
                       "Snarl based\n total time: " + str(totalTimeNew)))       
    
    timesPanel.set_xlabel("Time")
    timesPanel.set_ylabel("Count")
    timesPanel.set_title("Run times")
    timesPanel.set_xlim(minTime - 1, maxTime * 1.1)
    timesPanel.set_ylim(0, maxCount *1.01)


    plt.savefig("distanceProfile", dpi=600)


main()
