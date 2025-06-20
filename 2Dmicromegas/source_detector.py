import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from vmm_tools import read_hit, read_cluster


'''takes data from ROOT file(s) and plots the number of hits and clusters vs time to check when the source was present'''

rootFolder = "VMMDataTestROOTFilesMay7UpdatedGeometry" #folder containing the ROOT files
rootFiles = sorted(glob.glob(os.path.join(rootFolder, "*.root"))) #using the sorted feature assuming the filenames have a meaning (e.g., chronological)

#given a folder full of chronologically sorted ROOT files (by filename), plot the number of hits vs time and number of clusters vs time
def plotVsTime(rootFiles):
    time = np.arange(1, len(rootFiles)+1)*0.5 #the 0.5 assumes you will be plotting time in hours, and each data file has 30 mins of data
    numHits = []
    numClusters = []
    for filePath in rootFiles:
        dfHits = read_hit(filePath)
        numHits = np.append(numHits, len(dfHits))
        
        dfClusters = read_cluster(filePath)
        numClusters = np.append(numClusters, len(dfClusters))
    
    fig = plt.figure()
    plt.plot(time, numHits, color='blue')
    plt.xlabel(f'Time (hours)')
    plt.ylabel('Number of Hits')
    plt.savefig(f'hits_vs_time_plot.png')
    plt.close()

    fig = plt.figure()
    plt.plot(time, numClusters, color='blue')
    plt.xlabel(f'Time (hours)')
    plt.ylabel('Number of Clusters')
    plt.savefig(f'clusters_vs_time_plot.png')
    plt.close()

plotVsTime(rootFiles)