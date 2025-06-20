import numpy as np
import matplotlib.pyplot as plt


'''purpose of this script is to find broken channels by identifying which ones have pedestals that are way too high'''

rawPedestalCSV = np.transpose(np.genfromtxt('pedestal.csv', delimiter=',')) #replace the filename with whichever pedestal csv file you are analyzing
rawThresholdMeasureCSV = np.transpose(np.genfromtxt('threshold.csv', delimiter=',')) #replace the filename with whichever threshold csv file you are analyzing

#get parallel arrays of the vmm being measured, the channel, and the measured pedestal value
vmm = rawPedestalCSV[2]
channel = rawPedestalCSV[3]
pedestal = rawPedestalCSV[4]
threshold = rawThresholdMeasureCSV[4]    

#make plots of pedestal vs channel for each vmm and CSV of the bad channels (>300 mV pedestal)
badChannelVMMNumber = np.array(())
badChannels = np.array(())
for i in range(0, 16):
    indices = np.where(vmm == i)[0]

    vmmChannels = channel[indices]
    vmmPedestal = pedestal[indices]
    vmmThreshold = threshold[indices]
    vmmTheoreticalThreshold = vmmPedestal + 100.0

    fig = plt.figure()
    linewidthValue = 0.75
    plt.plot(vmmChannels, vmmPedestal, color='blue', label='Pedestal', linewidth=linewidthValue)
    plt.plot(vmmChannels, vmmThreshold, color='red', label='Measured Threshold', linewidth=linewidthValue)
    plt.plot(vmmChannels, vmmTheoreticalThreshold, color='green', label='Theoretical Threshold', linewidth=linewidthValue)
    plt.legend(loc='upper center')
    plt.xlim(0, 63)
    plt.ylim(0, 1500)
    plt.xlabel(f'Channel Number')
    plt.ylabel('Voltage (mV)')
    plt.title(f'VMM {i}')
    plt.savefig(f'pedestal_threshold_vmm{i}.png')
    plt.close()

    badChannelIndices = np.where((vmmPedestal > 250) | (vmmPedestal > vmmThreshold))[0]
    test = vmmChannels[badChannelIndices]

    badChannelVMMNumber = np.append(badChannelVMMNumber, np.full(len(test), i))
    badChannels = np.append(badChannels, test)

columnHeaders = "VMM,Channel"
badChannelsCSV = np.transpose(np.array((badChannelVMMNumber, badChannels)))
np.savetxt("bad_channels.csv", badChannelsCSV, delimiter=",", header=columnHeaders, comments='')