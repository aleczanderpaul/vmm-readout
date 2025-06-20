import numpy as np
import matplotlib.pyplot as plt
from vmm_tools import read_hit


'''purpose of this script is to find noisy channels by identifying which ones produce significantly more events than the others (for no other apparent reason)'''

dataFrame = read_hit('SourceDetectorData/05-05-25_Source.root') #replace the filename with whichever ROOT file you are analyzing

#this function has two purposes: make a histogram of the number of events in each VMM channel, and print the noisy channels that need masking in the slow control
def channelPlotterPerVMM(dataFrame, vmmID, tooManyCounts):
    vmmNumber = dataFrame['vmm']
    channelNumber = dataFrame['ch']

    hitInVMMIndices = np.where(vmmNumber == vmmID)[0]
    channelNumber = channelNumber[hitInVMMIndices]

    #makes and saves a histogram of the events per channel for the VMM (that matches the number of vmmID)
    fig = plt.figure()

    binning = np.arange(-0.5, 64.5, 1)
    counts, bin_edges, _ = plt.hist(channelNumber, bins=binning, color='red', histtype='step')
    
    print(vmmID, (bin_edges[np.where(counts > tooManyCounts)[0]]+0.5).astype(int)) #prints the VMM and the channels that have counts higher than tooManyCounts (these are noisy channels)

    plt.xlabel(f'ch')
    plt.ylabel('Events')
    plt.title(f'VMM {vmmID}')
    plt.savefig(f'hits_per_channel_hist_VMM{vmmID}.png')
    plt.close()

#plots events per channel histogram for each of the VMMs in the micromegas
for i in range(0, 16):
    channelPlotterPerVMM(dataFrame, i, 1000)