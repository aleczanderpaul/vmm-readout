import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
matplotlib.use('Agg')
import pandas as pd
from vmm_tools import combineDataFrames, fitCB , fiducializeArea, combineDataFramesMajd
import os
import glob


'''Reconstruct Fe55 source data from VMMs, single run. Has more features than overlaid version.'''

strip_edges = np.arange(-0.5,499.5,1.0)

'''-------------------------------------------------------------------'''
'''--------------------BEGIN HIT RELATED FUNCTIONS--------------------'''
'''-------------------------------------------------------------------'''
 # Isolate x and y hit rate and plot histograms
def plotXAndYHitRate(df_hits, strip_edges, data_duration, logscale=True):
    xHits = df_hits.loc[df_hits['plane'] == 0].reset_index()
    yHits = df_hits.loc[df_hits['plane'] == 1].reset_index()

    fig = plt.figure()
    plt.hist(xHits['pos'], bins=strip_edges, histtype='bar', color='blue', weights=np.full(len(xHits), 1/data_duration))
    plt.xlabel("strips x")
    plt.ylabel("counts / s")
    if logscale == True:
        plt.yscale("log")
    plt.savefig(f'Micromegas/plots/x_hit_rate.png', bbox_inches="tight")
    plt.close()

    fig = plt.figure()
    plt.hist(yHits['pos'], bins=strip_edges, histtype='bar', color='orange', weights=np.full(len(yHits), 1/data_duration))
    plt.xlabel("strips y")
    plt.ylabel("counts / s")
    if logscale == True:
        plt.yscale("log")
    plt.savefig(f'Micromegas/plots/y_hit_rate.png', bbox_inches="tight")
    plt.close()
'''-------------------------------------------------------------------'''
'''---------------------END HIT RELATED FUNCTIONS---------------------'''
'''-------------------------------------------------------------------'''

'''-------------------------------------------------------------------'''
'''------------------BEGIN CLUSTER RELATED FUNCTIONS------------------'''
'''-------------------------------------------------------------------'''
#plot cluster positions in a 2D histogram
def plotClusterLocations2D(df_clusters, strip_edges, data_duration, logscale=False):
    pos0Data = df_clusters['pos0']
    pos1Data = df_clusters['pos1']

    fig = plt.figure()
    if logscale == True:
        hist, xedges, yedges, im = plt.hist2d(pos0Data, pos1Data, bins=[strip_edges,strip_edges], cmap=plt.cm.jet, norm=LogNorm())
    elif logscale == False:
        hist, xedges, yedges, im = plt.hist2d(pos0Data, pos1Data, bins=[strip_edges,strip_edges], cmap=plt.cm.jet, cmin=1) #cmin=1 is to mask 0 count bins when linearscale)
    im.set_array(hist.T / data_duration)  # Convert to counts per second
    im.autoscale()  # Rescale color limits
    fiducialize_line_color = 'black'
    plt.axvline(156,linestyle = "--", color=fiducialize_line_color)
    plt.axvline(217,linestyle = "--", color=fiducialize_line_color)
    plt.axvline(280,linestyle = "--", color=fiducialize_line_color)
    plt.axvline(342,linestyle = "--", color=fiducialize_line_color)
    plt.axhline(156,linestyle = "--", color=fiducialize_line_color)
    plt.axhline(217,linestyle = "--", color=fiducialize_line_color)
    plt.axhline(280,linestyle = "--", color=fiducialize_line_color)
    plt.axhline(342,linestyle = "--", color=fiducialize_line_color)
    plt.text(170.5,170.5,"(c)",color=fiducialize_line_color, fontsize="large", fontweight="bold")
    plt.text(295,295,"(b)",color=fiducialize_line_color, fontsize="large", fontweight="bold")
    plt.text(170.5,295,"(d)",color=fiducialize_line_color, fontsize="large", fontweight="bold")
    plt.text(295,170.5,"(a)",color=fiducialize_line_color, fontsize="large", fontweight="bold")
    cbar1 = plt.colorbar()
    cbar1.set_label("Counts / s", rotation=270, labelpad=15) 
    plt.xlabel(f'x strip #')
    plt.ylabel(f'y strip #')
    plt.savefig(f'Micromegas/plots/cluster_locations.png', bbox_inches="tight")
    plt.close()

#compute number of electrons in event to find and plot gain, then do a best fit to the data
# Per Lucian 1 ADC ~ 1 mV
def plotGainAndFits(df_clusters, x_y_gain, fiducialize=False, fid_area='a', fit=True): #change area to desired section of the micromegas, see vmm_tools.py for options
    #6240 comes from 1 fC = 6240 electrons
    df_clusters["electrons_x"] = df_clusters['adc0'].apply(lambda x: 6240 * ( x  / x_y_gain[0] ) ) # 9 mV/fC is VMM gain setting for x channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values
    df_clusters["electrons_y"] = df_clusters['adc1'].apply(lambda x: 6240 * ( x / x_y_gain[1] ) ) # 4.5 mV/fC is VMM gain setting for y channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values
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
    
    if fit == True:
        fig = plt.figure()
        fitCB(df_clusters, plot=True, saveFig=False)
        plt.legend(loc='upper right')
        plt.savefig(f'Micromegas/plots/gain_fits_{areaName}.png', bbox_inches="tight")
        plt.close()
    elif fit == False:
        fig = plt.figure()
        gain = df_clusters['gain']
        xmin, xmax = 0, gain.max()
        nbins = 100
        plt.hist(gain,nbins,(xmin,xmax), density = True, color='g',alpha=0.6, label='Source Present')
        plt.legend(loc='upper right')
        plt.xlabel("Gain")
        plt.ylabel("Probability Density")
        plt.savefig(f'Micromegas/plots/gain_hist_{areaName}.png', bbox_inches="tight")
        plt.close()
    else:
        raise Exception("Pick a valid value for fit, either True or False")
'''-------------------------------------------------------------------'''
'''-------------------END CLUSTER RELATED FUNCTIONS-------------------'''
'''-------------------------------------------------------------------'''


'''Main execution'''
if __name__ == '__main__':
    rootFolder = "Micromegas/16mV-fC_overnight_both_calib" #folder containing the ROOT files
    df_hits, df_clusters = combineDataFrames(rootFolder)
    data_duration = len(glob.glob(os.path.join(rootFolder, "*.root"))) * 10 * 60 #duration of data in seconds, change as needed
    x_y_gain = (16.0, 4.5) #in mV/fC, change as needed

    plotClusterLocations2D(df_clusters, strip_edges, data_duration, logscale=True)
    plotXAndYHitRate(df_hits, strip_edges, data_duration, logscale=True)
    plotGainAndFits(df_clusters, x_y_gain, fiducialize=True, fid_area='a', fit=True)
    plotGainAndFits(df_clusters, x_y_gain, fiducialize=True, fid_area='b', fit=True)
    plotGainAndFits(df_clusters, x_y_gain, fiducialize=True, fid_area='c', fit=True)
    plotGainAndFits(df_clusters, x_y_gain, fiducialize=True, fid_area='d', fit=True)