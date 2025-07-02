import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from vmm_tools import combineDataFrames, fiducializeArea


'''Reconstruct Fe55 source data from VMMs, overlay data with Fe55 present and without Fe55 present'''

rootFolder = "May7SourceData" #folder containing the ROOT files
df_hits, df_clusters = combineDataFrames(rootFolder)
data_duration = 46 * 30 * 60 #duration of data in seconds, change as needed
strip_edges = np.arange(-0.5,499.5,1.0)

rootFolderNoSource = "May7NoSourceData" #folder containing the ROOT files with no Fe55 source
df_hits_no_source, df_clusters_no_source = combineDataFrames(rootFolderNoSource)
data_duration_no_source = 32 * 30 * 60 #duration of data in seconds, change as needed


 # Isolate x and y hits and plot histograms
def plotXAndYHits(df_hits, df_hits_no_source, strip_edges, data_duration, data_duration_no_source):
    xHits = df_hits.loc[df_hits['plane'] == 0].reset_index()
    yHits = df_hits.loc[df_hits['plane'] == 1].reset_index()

    xHitsNoSource = df_hits_no_source.loc[df_hits_no_source['plane'] == 0].reset_index()
    yHitsNoSource = df_hits_no_source.loc[df_hits_no_source['plane'] == 1].reset_index()

    fig = plt.figure()
    plt.hist(xHits['pos'], bins=strip_edges, histtype='step', color='blue', weights=np.full(len(xHits), 1/data_duration), label='Fe55 Present')
    plt.hist(xHitsNoSource['pos'], bins=strip_edges, histtype='step', color='red', weights=np.full(len(xHitsNoSource), 1/data_duration_no_source), label='No Fe55 Present')
    plt.xlabel("strips x")
    plt.ylabel("counts / s")
    plt.yscale("log")
    plt.legend(loc='upper right')
    plt.savefig(f'x_hits_overlaid.png', bbox_inches="tight")
    plt.close()

    fig = plt.figure()
    plt.hist(yHits['pos'], bins=strip_edges, histtype='step', color='orange', weights=np.full(len(yHits), 1/data_duration), label='Fe55 Present')
    plt.hist(yHitsNoSource['pos'], bins=strip_edges, histtype='step', color='green', weights=np.full(len(yHitsNoSource), 1/data_duration_no_source), label='No Fe55 Present')
    plt.xlabel("strips x")
    plt.xlabel("strips y")
    plt.ylabel("counts / s")
    plt.yscale("log")
    plt.legend(loc='upper right')
    plt.savefig(f'y_hits_overlaid.png', bbox_inches="tight")
    plt.close()

#compute number of electrons in event to find and plot gain, then do a best fit to the data
# Per Lucian 1 ADC ~ 1 mV
def plotGain(df_clusters, df_clusters_no_source, fiducialize=False):
    #6240 comes from 1 fC = 6240 electrons
    df_clusters["electrons_x"] = df_clusters['adc0'].apply(lambda x: 6240 * ( x  / 9.0 ) ) # 9 mV/fC is VMM gain setting for x channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values
    df_clusters["electrons_y"] = df_clusters['adc1'].apply(lambda x: 6240 * ( x / 4.5 ) ) # 4.5 mV/fC is VMM gain setting for y channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values
    df_clusters["electrons"] = df_clusters['electrons_x'] + df_clusters['electrons_y']
    df_clusters["gain"] = df_clusters["electrons"] / 167.5 #167.5 is the average number of primary electrons created by a 5.9 keV X-ray in Ar/CO2 70:30

    df_clusters_no_source["electrons_x"] = df_clusters_no_source['adc0'].apply(lambda x: 6240 * ( x  / 9.0 ) ) # 9 mV/fC is VMM gain setting for x channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values
    df_clusters_no_source["electrons_y"] = df_clusters_no_source['adc1'].apply(lambda x: 6240 * ( x / 4.5 ) ) # 4.5 mV/fC is VMM gain setting for y channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values
    df_clusters_no_source["electrons"] = df_clusters_no_source['electrons_x'] + df_clusters_no_source['electrons_y']
    df_clusters_no_source["gain"] = df_clusters_no_source["electrons"] / 167.5 #167.5 is the average number of primary electrons created by a 5.9 keV X-ray in Ar/CO2 70:30

    if fiducialize == True:
        df_clusters = fiducializeArea(df_clusters, area='bottom right') #change area to desired section of the micromegas, see vmm_tools.py for options
        df_clusters_no_source = fiducializeArea(df_clusters_no_source, area='bottom right') #change area to desired section of the micromegas, see vmm_tools.py for options
        pass
    elif fiducialize == False:
        pass
    else:
        raise Exception("Pick a valid value for fiducialize, either True or False")
    
    fig = plt.figure()
    gain = df_clusters['gain']
    xmin, xmax = 0, gain.max()
    nbins = 100
    plt.hist(gain,nbins,(xmin,xmax), density = True, color='orange', histtype='step', label='Fe55 Present')

    gain_no_source = df_clusters_no_source['gain']
    xmin_no_source, xmax_no_source = 0, gain_no_source.max()
    nbins = 100
    plt.hist(gain_no_source,nbins,(xmin_no_source,xmax_no_source), density = True, color='green', histtype='step', label='No Fe55 Present')
    plt.legend(loc='upper right')
    plt.xlabel("Gain")
    plt.ylabel("Probability Density")
    plt.savefig('gain_hist_overlaid.png', bbox_inches="tight")
    plt.close()


'''Main execution'''
plotXAndYHits(df_hits, df_hits_no_source, strip_edges, data_duration, data_duration_no_source)
plotGain(df_clusters, df_clusters_no_source, fiducialize=False)