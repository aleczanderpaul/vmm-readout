import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
matplotlib.use('Agg')
import pandas as pd
from vmm_tools import combineDataFrames, fiducializeArea
import os
import glob


'''Reconstruct Fe55 source data from VMMs and overlay data from different runs'''

strip_edges = np.arange(-0.5,499.5,1.0)
hist_colors = ['blue', 'orange', 'green', 'red', 'purple'] #colors for the histograms, if plotting more than 5 datasets overlaid, add more colors

#fetch the hits, clusters, and duration of data taking from multiple root folders, keeping them separated so they can be overlaid on plots later
def getHitsClustersAndDataDuration(rootFolders, single_file_duration_list):
    df_hits_list = [] #list to hold dataframes of hits
    df_clusters_list = [] #list to hold dataframes of clusters
    data_duration_list = [] #list to hold data durations for each root folder
    for i in range(0, len(rootFolders)):
        df_hits, df_clusters = combineDataFrames(rootFolders[i])
        data_duration = len(glob.glob(os.path.join(rootFolders[i], "*.root"))) * single_file_duration_list[i] * 60 #converts duration of all data files combined into seconds

        df_hits_list.append(df_hits)
        df_clusters_list.append(df_clusters)
        data_duration_list.append(data_duration)
    
    return df_hits_list, df_clusters_list, data_duration_list

'''-------------------------------------------------------------------'''
'''--------------------BEGIN HIT RELATED FUNCTIONS--------------------'''
'''-------------------------------------------------------------------'''
#run matplotlib histogram function for x hit rates
def histXHitRate(df_hits, strip_edges, data_duration, hist_color, hist_label):
    xHits = df_hits.loc[df_hits['plane'] == 0].reset_index()
    plt.hist(xHits['pos'], bins=strip_edges, histtype='step', color=hist_color, weights=np.full(len(xHits), 1/data_duration), label=hist_label)

#run matplotlib histogram function for y hit rates
def histYHitRate(df_hits, strip_edges, data_duration, hist_color, hist_label):
    yHits = df_hits.loc[df_hits['plane'] == 1].reset_index()
    plt.hist(yHits['pos'], bins=strip_edges, histtype='step', color=hist_color, weights=np.full(len(yHits), 1/data_duration), label=hist_label)

#plot overlaid x and y hit rates and save to png
def plotXYHitRates(df_hits_list, strip_edges, data_duration_list, hist_colors, hist_labels, logscale=True):
    #plot overlaid x strip hit rate
    fig = plt.figure()
    for i in range(0, len(rootFolders)):
        histXHitRate(df_hits_list[i], strip_edges, data_duration_list[i], hist_colors[i], hist_labels[i])
    plt.legend(loc='upper right')
    plt.xlabel("strips x")
    plt.ylabel("counts / s")
    if logscale == True:
        plt.yscale("log")
    plt.savefig(f'Micromegas/plots/x_hit_rate_overlaid.png', bbox_inches="tight")
    plt.close()

    #plot overlaid y strip hit rate
    fig = plt.figure()
    for i in range(0, len(rootFolders)):
        histYHitRate(df_hits_list[i], strip_edges, data_duration_list[i], hist_colors[i], hist_labels[i])
    plt.legend(loc='upper right')
    plt.xlabel("strips y")
    plt.ylabel("counts / s")
    if logscale == True:
        plt.yscale("log")
    plt.savefig(f'Micromegas/plots/y_hit_rate_overlaid.png', bbox_inches="tight")
    plt.close()
'''-------------------------------------------------------------------'''
'''---------------------END HIT RELATED FUNCTIONS---------------------'''
'''-------------------------------------------------------------------'''

'''-------------------------------------------------------------------'''
'''------------------BEGIN CLUSTER RELATED FUNCTIONS------------------'''
'''-------------------------------------------------------------------'''
#compute number of electrons in event to find and plot gain, then do a best fit to the data
# Per Lucian 1 ADC ~ 1 mV
#run matplotlib histogram function for avalanche gain
def histGain(df_clusters, gain_color, gain_label, x_gain, y_gain, fiducialize=False, fid_area='a', data_duration=None): #change area to desired section of the micromegas, see vmm_tools.py for options
    #6240 comes from 1 fC = 6240 electrons
    df_clusters["electrons_x"] = df_clusters['adc0'].apply(lambda x: 6240 * ( x  / x_gain ) ) # 9 mV/fC is VMM gain setting for x channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values
    df_clusters["electrons_y"] = df_clusters['adc1'].apply(lambda x: 6240 * ( x / y_gain ) ) # 4.5 mV/fC is VMM gain setting for y channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values
    df_clusters["electrons"] = df_clusters['electrons_x'] + df_clusters['electrons_y']
    df_clusters["gain"] = df_clusters["electrons"] / 167.5 #167.5 is the average number of primary electrons created by a 5.9 keV X-ray in Ar/CO2 70:30

    if fiducialize == True:
        areaName = fid_area
        df_clusters = fiducializeArea(df_clusters, area=areaName)
        pass
    elif fiducialize == False:
        areaName = 'all_areas'
        pass
    else:
        raise Exception("Pick a valid value for fiducialize, either True or False")
    
    gain = df_clusters['gain']
    #xmin, xmax = 0, gain.max()
    xmin, xmax = 2000, 15000
    nbins = 100
    #plt.hist(gain,nbins,(xmin,xmax), density = True, color=gain_color, histtype='step', label=gain_label) #if you want probability density as the y-axis, use this line
    plt.hist(gain,nbins,(xmin,xmax), color=gain_color, histtype='step', label=gain_label, weights=np.full(len(gain), 1/data_duration)) #if you want avg gain hit rate on the y-axis, use this line

#plot overlaid gains for each region, separating by pre-amp gain
def plotGainByRegion(df_clusters_list, data_duration_list, x_gains_list, y_gains_list, hist_colors, hist_labels):
    fid_areas = ['a', 'b', 'c', 'd'] #must be valid areas for fiducializeArea function in vmm_tools.py
    for l in range(0, len(fid_areas)):
        fig = plt.figure()
        for i in range(0, len(rootFolders)):
            histGain(df_clusters_list[i], hist_colors[i], hist_labels[i], x_gains_list[i], y_gains_list[i], fiducialize=True, fid_area=fid_areas[l], data_duration=data_duration_list[i])
        plt.legend(loc='upper right')
        plt.xlabel("Gain")
        plt.ylabel("Probability Density")
        plt.ylabel("counts / s") #if you want avg gain hit rate on the y-axis, use this line
        plt.ylim(0, 0.004) #to get all files to have the same y-axis range for easier comparison, change as necessary
        plt.xlim(2000, 15000) #to get all files to have the same x-axis range for easier comparison, change as necessary
        plt.savefig(f'Micromegas/plots/gain_hist_{fid_areas[l]}_overlaid.png', bbox_inches="tight")
        plt.close()

#plot overlaid gains for each pre-amp gain, separating for each region
def plotGainByPreAmpGain(df_clusters_list, data_duration_list, x_gains_list, y_gains_list, hist_colors):
    for l in range(0, len(rootFolders)):
        fig = plt.figure()
        fid_areas = ['a', 'b', 'c', 'd']
        area_labels = ['Region a', 'Region b', 'Region c', 'Region d']
        for i in range(0, len(fid_areas)):
                histGain(df_clusters_list[l], hist_colors[i], area_labels[i], x_gains_list[l], y_gains_list[l], fiducialize=True, fid_area=fid_areas[i], data_duration=data_duration_list[l])
        plt.legend(loc='upper right')
        plt.xlabel("Gain")
        plt.ylabel("Probability Density")
        plt.ylabel("counts / s") #if you want avg gain hit rate on the y-axis, use this line
        plt.ylim(0, 0.004) #to get all files to have the same y-axis range for easier comparison, change as necessary
        plt.xlim(2000, 15000) #to get all files to have the same x-axis range for easier comparison, change as necessary
        plt.savefig(f'Micromegas/plots/gain_hist_{x_gains_list[l]}_mVfC_overlaid.png', bbox_inches="tight")
        plt.close()

#calculate the charge sharing (defined as # electrons in x / # electrons in y)
def calculateChargeSharing(df_clusters, x_gain, y_gain, fiducialize=False, fid_area='a'):
    df_clusters["electrons_x"] = df_clusters['adc0'].apply(lambda x: 6240 * ( x  / x_gain ) ) # 9 mV/fC is VMM gain setting for x channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values
    df_clusters["electrons_y"] = df_clusters['adc1'].apply(lambda x: 6240 * ( x / y_gain ) ) # 4.5 mV/fC is VMM gain setting for y channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values

    if fiducialize == True:
        areaName = fid_area
        df_clusters = fiducializeArea(df_clusters, area=areaName)
        pass
    elif fiducialize == False:
        areaName = 'all_areas'
        pass
    else:
        raise Exception("Pick a valid value for fiducialize, either True or False")
    
    charge_sharing = 1.0*np.mean(df_clusters["electrons_x"] / df_clusters["electrons_y"])

    return charge_sharing

#print the calculated charge sharing values for each region and their pre-amp gain
def getChargeSharingPerRegion(df_clusters_list, x_gains_list, y_gains_list):
    for l in range(0, len(rootFolders)):
        fid_areas = ['a', 'b', 'c', 'd']
        for i in range(0, len(fid_areas)):
            charge_sharing = calculateChargeSharing(df_clusters_list[l], x_gains_list[l], y_gains_list[l], fiducialize=True, fid_area=fid_areas[i])
            print(f'Charge sharing for area {fid_areas[i]} with x={x_gains_list[l]} mV/fC and y={y_gains_list[l]} mV/fC preamp gain: {charge_sharing}')
'''-------------------------------------------------------------------'''
'''-------------------END CLUSTER RELATED FUNCTIONS-------------------'''
'''-------------------------------------------------------------------'''

'''Main execution'''
if __name__ == "__main__":
    rootFolders = ['Micromegas/16mV-fC_overnight_both_calib', 'Micromegas/DeprecatedData/12mV-fC_overnight_noise_removed'] #folder paths containing the ROOT files
    x_gains_list = [16.0, 12.0] #list of x preamp gains in mV/fC, must be parallel with rootFolders
    y_gains_list = [4.5, 4.5] #list of y preamp gains in mV/fC, must be parallel with rootFolders
    single_file_duration_list = [10, 10] #list of data-taking durations for a single file in minutes, must be parallel with rootFolders
    hist_labels = ['16mV-fC', '12mV-fC'] #labels for the histograms, must be parallel with rootFolders

    df_hits_list, df_clusters_list, data_duration_list = getHitsClustersAndDataDuration(rootFolders, single_file_duration_list)

    plotXYHitRates(df_hits_list, strip_edges, data_duration_list, hist_colors, hist_labels, logscale=True)

    plotGainByRegion(df_clusters_list, data_duration_list, x_gains_list, y_gains_list, hist_colors, hist_labels)
    plotGainByPreAmpGain(df_clusters_list, data_duration_list, x_gains_list, y_gains_list, hist_colors)
    getChargeSharingPerRegion(df_clusters_list, x_gains_list, y_gains_list)