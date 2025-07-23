import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from vmm_tools import combineDataFrames, fiducializeArea
import os
import glob


'''Reconstruct Fe55 source data from VMMs and overlay data from different runs'''

strip_edges = np.arange(-0.5,499.5,1.0)

def plotXHitRate(df_hits, strip_edges, data_duration, hist_color, hist_label):
    xHits = df_hits.loc[df_hits['plane'] == 0].reset_index()
    plt.hist(xHits['pos'], bins=strip_edges, histtype='step', color=hist_color, weights=np.full(len(xHits), 1/data_duration), label=hist_label)

def plotYHitRate(df_hits, strip_edges, data_duration, hist_color, hist_label):
    yHits = df_hits.loc[df_hits['plane'] == 1].reset_index()
    plt.hist(yHits['pos'], bins=strip_edges, histtype='step', color=hist_color, weights=np.full(len(yHits), 1/data_duration), label=hist_label)

#compute number of electrons in event to find and plot gain, then do a best fit to the data
# Per Lucian 1 ADC ~ 1 mV
def plotGain(df_clusters, gain_color, gain_label, x_gain, y_gain, fiducialize=False, fid_area='a'): #change area to desired section of the micromegas, see vmm_tools.py for options
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
    plt.hist(gain,nbins,(xmin,xmax), density = True, color=gain_color, histtype='step', label=gain_label)


'''Main execution'''
if __name__ == "__main__":
    rootFolders = ['Micromegas/16mV-fC', 'Micromegas/12mV-fC', 'Micromegas/9mV-fC'] #folder paths containing the ROOT files
    hist_labels = ['16 mV/fC', '12 mV/fC', '9 mV/fC'] #labels for the histograms, must be parallel with rootFolders
    hist_colors = ['blue', 'orange', 'green', 'red', 'purple'] #colors for the histograms, must be parallel with hist_labels or longer
    file_name_append = '_preamp_gain_test' #append to file names, change as needed

    df_hits_list = [] #list to hold dataframes of hits
    df_clusters_list = [] #list to hold dataframes of clusters

    single_file_duration_list = [10, 10, 10] #list of durations of a single file in minutes, must be parallel with rootFolders
    data_duration_list = [] #list to hold data durations for each root folder
    for i in range(0, len(rootFolders)):
        df_hits, df_clusters = combineDataFrames(rootFolders[i])
        data_duration = len(glob.glob(os.path.join(rootFolders[i], "*.root"))) * single_file_duration_list[i] * 60 #duration of data in seconds, change as needed

        df_hits_list.append(df_hits)
        df_clusters_list.append(df_clusters)
        data_duration_list.append(data_duration)
    
    #plot overlaid x strip hit rate
    fig = plt.figure()
    for i in range(0, len(rootFolders)):
        plotXHitRate(df_hits_list[i], strip_edges, data_duration_list[i], hist_colors[i], hist_labels[i])
    plt.legend(loc='upper right')
    plt.xlabel("strips x")
    plt.ylabel("counts / s")
    #plt.yscale("log")
    plt.savefig(f'Micromegas/plots/x_hit_rate_overlaid{file_name_append}.png', bbox_inches="tight")
    plt.close()

    #plot overlaid y strip hit rate
    fig = plt.figure()
    for i in range(0, len(rootFolders)):
        plotYHitRate(df_hits_list[i], strip_edges, data_duration_list[i], hist_colors[i], hist_labels[i])
    plt.legend(loc='upper right')
    plt.xlabel("strips y")
    plt.ylabel("counts / s")
    #plt.yscale("log")
    plt.savefig(f'Micromegas/plots/y_hit_rate_overlaid{file_name_append}.png', bbox_inches="tight")
    plt.close()

    #plot overlaid gains
    fig = plt.figure()
    x_gains_list = [16.0, 12.0, 9.0] #list of x preamp gains in mV/fC, must be parallel with rootFolders
    for i in range(0, len(rootFolders)):
        plotGain(df_clusters_list[i], hist_colors[i], hist_labels[i], x_gains_list[i], 4.5, fiducialize=False)
    plt.legend(loc='upper right')
    plt.xlabel("Gain")
    plt.ylabel("Probability Density")
    plt.savefig(f'Micromegas/plots/gain_hist_overlaid{file_name_append}.png', bbox_inches="tight")
    plt.close()