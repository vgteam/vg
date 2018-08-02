import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import sys

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
    distsPanel.plot([0,maxVal],[0, maxVal], linestyle="-", color = "black")
   
    distsPanel.set_xlabel("Paths based distance")
    distsPanel.set_ylabel("Snarl based distance")
    distsPanel.set_title("Predicted Distances")
    distsPanel.set_xlim(min(oldDists) * .95,  max(oldDists) * 1.05)
    distsPanel.set_ylim(min(myDists) * .95,  max(myDists) * 1.05)



    ############# Time plot panel

    timesPanel.plot(oldDists, oldTimes, marker = 'o', linewidth = 0, 
                    markerfacecolor = 'blue', 
                    markeredgewidth = 0, markersize = 2, label = "Path based") 
    timesPanel.plot(myDists, myTimes, marker = 'o', linewidth = 0, 
                    markerfacecolor = 'red', 
                    markeredgewidth = 0, markersize = 2, label = "Snarl based") 
    timesPanel.plot([0,max(max(myDists), max(oldDists)) * 1.05],[timeIndex, timeIndex], linestyle="-", 
                                    color = "black",label = "Time for index calculation")
   
    timesPanel.set_xlabel("Predicted distance")
    timesPanel.set_ylabel("Time")
    timesPanel.set_title("Run times")
    timesPanel.set_xlim(min(min(myDists), min(oldDists)) * 0.95, 
                        max(max(myDists), max(oldDists)) * 1.05)
    timesPanel.set_ylim(min(min(myTimes), min(oldTimes), timeIndex) * .95,  
                        max(max(myTimes), max(oldTimes), timeIndex) * 1.05)

    timesPanel.legend()

    plt.savefig("distanceProfile", dpi=600)


main()
