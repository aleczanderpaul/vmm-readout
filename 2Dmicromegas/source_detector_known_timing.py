import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from vmm_tools import read_hit


def plotKeywordHist(keyword, numBins):
    beforeSourceData = beforeSourceHitDataframe[keyword]
    sourceData = sourceHitDataframe[keyword]
    afterSourceData = afterSourceHitDataframe[keyword]

    combinedData = pd.concat([beforeSourceData, sourceData, afterSourceData])
    minimum = combinedData.min()
    maximum = combinedData.max()
    histRange = (minimum, maximum)

    fig = plt.figure()
    plt.hist(beforeSourceData, bins=numBins, color='blue', label='Before Source', histtype='step', range=histRange)
    plt.hist(sourceData, bins=numBins, color='red', label='Source', histtype='step', range=histRange)
    plt.hist(afterSourceData, bins=numBins, color='green', label='After Source', histtype='step', range=histRange)
    plt.yscale('log')
    plt.xlabel(f'{keyword}')
    plt.ylabel('Events')
    plt.legend(loc = 'upper left')
    plt.savefig(f'SourceDetectorPlots/{keyword}_hist.png')
    plt.close()

beforeSourceHitDataframe = read_hit('SourceDetectorData/05-05-25_BeforeSource.root')
sourceHitDataframe = read_hit('SourceDetectorData/05-05-25_Source.root')
afterSourceHitDataframe = read_hit('SourceDetectorData/05-05-25_AfterSource.root')

plotKeywordHist('id', 100)
plotKeywordHist('det', 100)
plotKeywordHist('plane', 100)
plotKeywordHist('fec', 100)
plotKeywordHist('vmm', 100)
plotKeywordHist('readout_time', 100)
plotKeywordHist('time', 100)
plotKeywordHist('ch', 100)
plotKeywordHist('pos', 100)
plotKeywordHist('bcid', 100)
plotKeywordHist('tdc', 100)
plotKeywordHist('adc', 100)
plotKeywordHist('over_threshold', 100)
plotKeywordHist('chip_time', 100)