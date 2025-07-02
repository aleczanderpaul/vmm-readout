import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from vmm_tools import combineDataFrames, fitCB, fitLorentzian, fitGaussian, fiducializeArea


'''Reconstruct Fe55 source data from VMMs'''

rootFolder = "May7SourceData" #folder containing the ROOT files
df_hits, df_clusters = combineDataFrames(rootFolder)
data_duration = 46 * 30 * 60 #duration of data in seconds, change as needed
strip_edges = np.arange(-0.5,499.5,1.0)


#plot cluster positions in a 2D histogram
def plotClusterLocations2D(df_clusters, strip_edges):
    pos0Data = df_clusters['pos0']
    pos1Data = df_clusters['pos1']

    fig = plt.figure()
    plt.hist2d(pos0Data, pos1Data, bins=[strip_edges,strip_edges], cmap=plt.cm.jet)
    plt.axvline(156,linestyle = "--", color="w")
    plt.axvline(217,linestyle = "--", color="w")
    plt.axvline(280,linestyle = "--", color="w")
    plt.axvline(342,linestyle = "--", color="w")
    plt.axhline(156,linestyle = "--", color="w")
    plt.axhline(217,linestyle = "--", color="w")
    plt.axhline(280,linestyle = "--", color="w")
    plt.axhline(342,linestyle = "--", color="w")
    plt.text(170.5,170.5,"(c)",color="w", fontsize="large", fontweight="bold")
    plt.text(295,295,"(b)",color="w", fontsize="large", fontweight="bold")
    plt.text(170.5,295,"(d)",color="w", fontsize="large", fontweight="bold")
    plt.text(295,170.5,"(a)",color="w", fontsize="large", fontweight="bold")
    cbar1 = plt.colorbar()
    cbar1.set_label("Count", rotation=270, labelpad=15) 
    plt.xlabel(f'pos0')
    plt.ylabel(f'pos1')
    plt.savefig(f'cluster_locations.png', bbox_inches="tight")
    plt.close()


 # Isolate x and y hits and plot histograms
def plotXAndYHits(df_hits, strip_edges, data_duration):
    xHits = df_hits.loc[df_hits['plane'] == 0].reset_index()
    yHits = df_hits.loc[df_hits['plane'] == 1].reset_index()

    fig = plt.figure()
    plt.hist(xHits['pos'], bins=strip_edges, histtype='bar', color='blue', weights=np.full(len(xHits), 1/data_duration))
    plt.xlabel("strips x")
    plt.ylabel("counts / s")
    plt.yscale("log")
    plt.savefig(f'x_hits.png', bbox_inches="tight")
    plt.close()

    fig = plt.figure()
    plt.hist(yHits['pos'], bins=strip_edges, histtype='bar', color='orange', weights=np.full(len(yHits), 1/data_duration))
    plt.xlabel("strips y")
    plt.ylabel("counts / s")
    plt.yscale("log")
    plt.savefig(f'y_hits.png', bbox_inches="tight")
    plt.close()

#compute number of electrons in event to find and plot gain, then do a best fit to the data
# Per Lucian 1 ADC ~ 1 mV
def plotGainAndFits(df_clusters, fiducialize=False, fit=True):
    #6240 comes from 1 fC = 6240 electrons
    df_clusters["electrons_x"] = df_clusters['adc0'].apply(lambda x: 6240 * ( x  / 9.0 ) ) # 9 mV/fC is VMM gain setting for x channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values
    df_clusters["electrons_y"] = df_clusters['adc1'].apply(lambda x: 6240 * ( x / 4.5 ) ) # 4.5 mV/fC is VMM gain setting for y channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values
    df_clusters["electrons"] = df_clusters['electrons_x'] + df_clusters['electrons_y']
    df_clusters["gain"] = df_clusters["electrons"] / 167.5 #167.5 is the average number of primary electrons created by a 5.9 keV X-ray in Ar/CO2 70:30

    if fiducialize == True:
        df_clusters = fiducializeArea(df_clusters, area='bottom right') #change area to desired section of the micromegas, see vmm_tools.py for options
        pass
    elif fiducialize == False:
        pass
    else:
        raise Exception("Pick a valid value for fiducialize, either True or False")
    
    if fit == True:
        fitCB(df_clusters, plot=True, saveFig=True)
        fitLorentzian(df_clusters, plot=True, saveFig=True)
        fitGaussian(df_clusters, plot=True, saveFig=True)

        fig = plt.figure()
        fitCB(df_clusters, plot=True, saveFig=False)
        fitLorentzian(df_clusters, plot=True, saveFig=False)
        fitGaussian(df_clusters, plot=True, saveFig=False)
        plt.legend(loc='upper right')
        plt.savefig('gain_fits.png', bbox_inches="tight")
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
        plt.savefig('gain_hist.png', bbox_inches="tight")
        plt.close()
    else:
        raise Exception("Pick a valid value for fit, either True or False")


'''Main execution'''
plotClusterLocations2D(df_clusters, strip_edges)
plotXAndYHits(df_hits, strip_edges, data_duration)
plotGainAndFits(df_clusters, fiducialize=False, fit=False)
plotGainAndFits(df_clusters, fiducialize=False, fit=True)