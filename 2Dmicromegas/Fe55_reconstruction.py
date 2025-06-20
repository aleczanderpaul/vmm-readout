import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from vmm_tools import combineDataFrames, fitCB, fitLorentzian, fitGaussian, fiducializeArea


'''Reconstruct Fe55 source data from VMMs'''

rootFolder = "VMMDataTestROOTFilesMay7UpdatedGeometry" #folder containing the ROOT files
df_hits, df_clusters = combineDataFrames(rootFolder)
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
    plt.savefig(f'cluster_locations.png')
    plt.close()


 # Isolate x and y hits and plot histograms
def plotXAndYHits(df_hits, strip_edges):
    xHits = df_hits.loc[df_hits['plane'] == 0].reset_index()
    yHits = df_hits.loc[df_hits['plane'] == 1].reset_index()

    fig = plt.figure()
    plt.hist(xHits['pos'], bins=strip_edges, histtype='bar', color='blue')
    plt.xlabel("strips x")
    plt.ylabel("count")
    plt.yscale("log")
    plt.savefig(f'x_hits.png')
    plt.close()

    fig = plt.figure()
    plt.hist(yHits['pos'], bins=strip_edges, histtype='bar', color='orange')
    plt.xlabel("strips y")
    plt.ylabel("count")
    plt.yscale("log")
    plt.savefig(f'y_hits.png')
    plt.close()

#compute number of electrons in event to find and plot gain, then do a best fit to the data
# Per Lucian 1 ADC ~ 1 mV
def plotGainAndFits(df_clusters, fiducialize=False):
    df_clusters["electrons_x"] = df_clusters['adc0'].apply(lambda x: 6240 * ( x  / 9.0 ) ) # 9 mV/fC is VMM gain setting for x channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values
    df_clusters["electrons_y"] = df_clusters['adc1'].apply(lambda x: 6240 * ( x / 4.5 ) ) # 4.5 mV/fC is VMM gain setting for y channels, 170mV is the pedestal, 1200mV is th operating voltage 1024 is the number of possible ADC values
    df_clusters["electrons"] = df_clusters['electrons_x'] + df_clusters['electrons_y']
    df_clusters["gain"] = df_clusters["electrons"] / 167.5

    if fiducialize == True:
        df_clusters = fiducializeArea(df_clusters, area='bottom right') #change area to desired section of the micromegas, see vmm_tools.py for options
        pass
    elif fiducialize == False:
        pass
    else:
        raise Exception("Pick a valid value for fiducialize, either True or False")
    
    fitCB(df_clusters)
    fitLorentzian(df_clusters)
    fitGaussian(df_clusters)

    fig = plt.figure()
    fitCB(df_clusters, plot=True, saveFig=False)
    fitLorentzian(df_clusters, plot=True, saveFig=False)
    fitGaussian(df_clusters, plot=True, saveFig=False)
    plt.legend(loc='upper right')
    plt.savefig('gain_fits.png')
    plt.close()


'''Main execution'''
#plotClusterLocations2D(df_clusters, strip_edges)
#plotXAndYHits(df_hits, strip_edges)
plotGainAndFits(df_clusters, fiducialize=True)