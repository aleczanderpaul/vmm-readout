import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


'''purpose of this script is to find broken channels by identifying which ones have pedestals that are way too high'''

#fetches the vmm #, channel #, pedestal, and measured threshold from CSV files
def fetch_pedestal_and_threshold_info(pedestalCSV, thresholdCSV):
    rawPedestalCSV = np.genfromtxt(pedestalCSV, delimiter=',', unpack=True)
    rawThresholdMeasureCSV = np.genfromtxt(thresholdCSV, delimiter=',', unpack=True)

    vmm = rawPedestalCSV[2]
    channel = rawPedestalCSV[3]
    pedestal = rawPedestalCSV[4]
    threshold = rawThresholdMeasureCSV[4]

    #takes out the headers inside the arrays which got converted to NaN values in the np.genfromtxt step
    vmm = vmm[~np.isnan(vmm)]
    channel = channel[~np.isnan(channel)]
    pedestal = pedestal[~np.isnan(pedestal)]
    threshold = threshold[~np.isnan(threshold)]

    return vmm, channel, pedestal, threshold

#make plots of the pedestal, threshold, and the theoretical threshold for each VMM
def plot_pedestal_and_threshold(vmm, channel, pedestal, threshold, theoreticalThreshold):
    for i in range(0, 4):
        indices = np.where(vmm == i)[0]

        vmmChannels = channel[indices]
        vmmPedestal = pedestal[indices]
        vmmThreshold = threshold[indices]
        vmmTheoreticalThreshold = vmmPedestal + theoreticalThreshold

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

#find the bad channels and return their number and corresponding VMM number
def find_bad_channels(vmm, channel, pedestal, threshold, theoreticalThreshold):
    badChannelIndices = np.where((pedestal > 250) | (pedestal > threshold) | (np.abs(threshold - pedestal) < 70) | (pedestal < 140))[0]
    badChannelVMMNumber = vmm[badChannelIndices]
    badChannels = channel[badChannelIndices]

    return badChannelVMMNumber, badChannels

#save the bad channels and their corresponding VMM number to a CSV file
def save_bad_channels_to_csv(badChannelVMMNumber, badChannels):
    columnHeaders = "VMM,Channel"
    badChannelsCSV = np.transpose(np.array((badChannelVMMNumber, badChannels)))
    np.savetxt("bad_channels.csv", badChannelsCSV, delimiter=",", header=columnHeaders, comments='')

#main function to run the script to make user-friendly
def main(pedestalCSV, thresholdCSV, theoreticalThreshold, plot=True, savetoCSV=True):
    vmm, channel, pedestal, threshold = fetch_pedestal_and_threshold_info(pedestalCSV, thresholdCSV)

    if plot == True:
        plot_pedestal_and_threshold(vmm, channel, pedestal, threshold, theoreticalThreshold)
    
    if savetoCSV == True:
        badChannelVMMNumber, badChannels = find_bad_channels(vmm, channel, pedestal, threshold, theoreticalThreshold)
        save_bad_channels_to_csv(badChannelVMMNumber, badChannels)

#runs the script
if __name__ == "__main__":
    main('Micromegas/Pedestal_111154.csv', 'Micromegas/Threshold_111138.csv', 100, plot=True, savetoCSV=True) #just need to edit the filenames, theoretical threshold value, and whether or not to plot and save to CSV