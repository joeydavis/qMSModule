# -*- coding: utf-8 -*-
"""
Created on Mon May  5 18:23:21 2014

@author: jhdavis
"""

import numpy
import qMS
import vizLib
import matplotlib.pyplot as plt
import pylab
from matplotlib import gridspec

def getAllOccupancyMRM(dataFrame, fileName, listOfProteins=None,
                       num=['light'], den=['heavy'], total=False,
                        normProtein=None, normValue=1.0, offset=0.0):
    if listOfProteins is None:
        listOfProteins = list(dataFrame['ProteinName'].unique())
    pDict = {}
    for p in listOfProteins:
        pDict[p] = getOccupancyMRM(dataFrame, fileName, p, num=num, den=den, 
                                    total=total, normValue=normValue, 
                                    offset=offset, allData=False).values
    if not (normProtein is None):
        normValue = 1/numpy.median(pDict[normProtein])
        for p in listOfProteins:
            pDict[p] = getOccupancyMRM(dataFrame, fileName, p, num=num, den=den, 
                                        total=total, normValue=normValue, 
                                        offset=offset, allData=False).values
    return pDict
   
def getAllOccupancyFileListMRM(dataFrame, fileList, listOfProteins=None, 
                               num=['light'], den=['heavy'], total=False, 
                                normProtein=None, normValue=1.0, offset=0.0):
    filePDict = {}
    for i in fileList:
        filePDict[i] = getAllOccupancyMRM(dataFrame, i, listOfProteins=listOfProteins, 
                                            num=num, den=den, total=total, normProtein=normProtein, 
                                            normValue=normValue, offset=offset)
    return filePDict

def getOccupancyMRM(dataFrame, fileName, proteinName, num=['light'], den=['heavy'], normValue=1.0, offset=0.0, total=False, allData=False):
    """getOccupancy takes a pandas dataframe generated from reading an MRM CSV file
        as well as a fileName (eg.wiff) and a a proteinName. It calculates a 
        protein occupancy using the fields specified in the numerator and 
        denominator lists.

    :param dataFrame: the dataFrame to work on. Should be generated from a 
        Skyline MRM export results CSV (e.g. dataFrame = pandas.read_csv(csvFile, na_values='#N/A'))
        Must bear the columns "FileName, ProteinName, Area, TotalArea, TotalBackground
    :type dataFrame: a pandas dataframe
    :param fileName: the file to consider (the wiff file)
    :type fileName: a string
    :param proteinName: the protein to consider
    :type proteinName: a string
    :param proteinName: the protein to consider
    :type proteinName: a string
    :param num: species to calc in the numerator
    :type num: a list of strings (light or heavy)
    :param den: species to calc in the denominator
    :type den: a list of strings (light or heavy)
    :param total: boolean to calculate the total or on a product by product basis
    :type total: boolean (defaults to false)
    :param allData: a boolean to return the full dataframe or just the calculated value
    :type allData: boolean (defaults to false)
    
    :returns:  a pandas dataframe with the calculated value (either appended or on its own)

    """
    allProducts = dataFrame[(dataFrame['FileName']==fileName) & (dataFrame['ProteinName']==proteinName)]
    
    if total:
        stringAppend = ' TotalArea'
    else:
        stringAppend = ' Area'
    
    allProducts['num'] = 0.0
    allProducts['den'] = 0.0
    
    for i in num:
        allProducts['num'] = allProducts['num'] + allProducts[i + stringAppend]
        if total:
            allProducts['num'] = allProducts['num'] - allProducts[i + ' TotalBackground']
    for i in den:
        allProducts['den'] = allProducts['den'] + allProducts[i + stringAppend]
        if total:
            allProducts['den'] = allProducts['den'] - allProducts[i + ' TotalBackground']
    
    allProducts['calcValue'] = (allProducts['num']/allProducts['den'])*normValue+offset

    if allData:
        return allProducts
    else:
        if total:
            return qMS.dropDuplicatesPandas(allProducts)['calcValue'].dropna()
        else:
            return allProducts['calcValue'].dropna()
            
def getInfoMRM(dataFrame, index):
    """getInfoMRM takes a pandas dataframe (read from the tsv output from skyline),
        and an index into that dataframe. Returns info related to that index 
        in the order [filename, mod. pep. seq, precharge, fragIon, prodCharge, isotope]

    :param dataFrame: pandas dataFrame - generated using tsv from skyline and the command
        dataFrame = pandas.read_csv(file, na_values='#N/A', sep='\t')
    :type dataFrame: pandas dataframe
    :param index: a int of what index to inspect
    :type index: int
    :returns:  a list of data, [filename, mod. pep. seq, precharge, fragIon, prodCharge, isotope]

    """
    fn = dataFrame.ix[index]['FileName']
    pepSeq = dataFrame.ix[index]['PeptideModifiedSequence']
    precursorCharge = dataFrame.ix[index]['PrecursorCharge']
    fragIon = dataFrame.ix[index]['FragmentIon']
    prodCharge = dataFrame.ix[index]['ProductCharge']
    isotope = dataFrame.ix[index]['IsotopeLabelType']
    return [fn, pepSeq, precursorCharge, fragIon, prodCharge, isotope]

def getTotalChromatographMRM(dataFrame, index, fragIon=None):
    [fn, pms, preC, fi, proC, isotope] = getInfoMRM(dataFrame, index)
    subFrame = dataFrame[(dataFrame['FileName'] == fn) &
                         (dataFrame['PeptideModifiedSequence'] == pms) &
                         (dataFrame['PrecursorCharge'] == preC) &
                         (dataFrame['IsotopeLabelType'] == isotope)]
    if fragIon is None:
        first = subFrame['Intensities'].values[0].split(',')
        totArray = numpy.array([float(i) for i in first])
        for a in subFrame['Intensities'].values[1:]:
            b = numpy.array([float(x) for x in a.split(',')])
            totArray = numpy.add(totArray,b)
    else:
        print "NEED TO WRITE THIS CODE TO DEAL WITH SUBSETS OF FRAGMENT IONS"
    
    return [numpy.array([float(x) for x in subFrame['Times'].values[0].split(',')]), totArray]

def getPairedIonMRM(dataFrame, index):
    [fn, pms, preC, fi, proC, isotope] = getInfoMRM(dataFrame, index)
    if isotope == 'light':
        isotope = 'heavy'
    else:
        isotope = 'light'
    return dataFrame[(dataFrame['FileName'] == fn) &
                         (dataFrame['PeptideModifiedSequence'] == pms) &
                         (dataFrame['PrecursorCharge'] == preC) &
                         (dataFrame['FragmentIon'] == fi) &
                         (dataFrame['ProductCharge'] == proC) &
                         (dataFrame['IsotopeLabelType'] == isotope)].index.values[0]

def getIsotopPairTotalsMRM(dataFrame, index):
    [fn, pms, preC, fi, proC, isotope] = getInfoMRM(dataFrame, index)
    indexOpp = getPairedIonMRM(dataFrame, index)
    if isotope == 'light':
        [lX, lY] = getTotalChromatographMRM(dataFrame, index)
        [hX, hY] = getTotalChromatographMRM(dataFrame, indexOpp)
    
    else:
        [lX, lY] = getTotalChromatographMRM(dataFrame, indexOpp)
        [hX, hY] = getTotalChromatographMRM(dataFrame, index)

    return [lX, lY, hX, hY]

def getRelatedIndeciesMRM(dataFrame, index):
    [fn, ps, preC, fi, proC, i] = getInfoMRM(dataFrame, index)
    return dataFrame[(dataFrame['IsotopeLabelType'] == i) &
                        (dataFrame['FileName'] == fn) &
                        (dataFrame['PeptideModifiedSequence'] == ps) &
                        (dataFrame['PrecursorCharge'] == preC)].index.values


def plotTotalChromPairsMRM(dataFrame, toPlotIndex, axis, 
                            colors=['blue', 'red'], smooth=0, zoom=False):
    [lightX, lightY, heavyX, heavyY] = getIsotopPairTotalsMRM(dataFrame, toPlotIndex)
    if smooth > 0:
        lightY = vizLib.smoothListGaussian([float(i) for i in lightY], degree=smooth)
        heavyY = vizLib.smoothListGaussian([float(i) for i in heavyY], degree=smooth)
    axis.plot(lightX, lightY, colors[0], label='light')
    axis.plot(heavyX, heavyY, colors[1], label='heavy')
    if zoom:
        maxIndex = numpy.argmax(numpy.array(lightY))
        axis.set_xlim(lightX[maxIndex]-2, lightX[maxIndex]+2)
    return axis

                        
def plotAllTransitionsMRM(dataFrame, index, a, colors=None, smooth=0, zoom=False):
    allIndicies = getRelatedIndeciesMRM(dataFrame, index)
    colors = vizLib.getBCs('q', min(len(allIndicies), 9))
    colors.append('grey')
    for i in range(len(allIndicies)):
        [fn, ps, preC, fi, proC, iso] = getInfoMRM(dataFrame, allIndicies[i])
        a = plotMRM_index(dataFrame, allIndicies[i], a, color=colors[i], 
                          smooth=smooth, zoom=zoom)
    return a

def plotMRM(dataFrame, fileName, pepSeq, precursorCharge, fragIon, prodCharge, 
            isotopeLabel, ax, color='grey', smooth=0, zoom=False):
    subDF = dataFrame[(dataFrame['PeptideModifiedSequence'] == pepSeq) & 
                        (dataFrame['FragmentIon'] == fragIon) & 
                        (dataFrame['PrecursorCharge'] == precursorCharge) & 
                        (dataFrame['IsotopeLabelType'] == isotopeLabel) & 
                        (dataFrame['ProductCharge'] == prodCharge) & 
                        (dataFrame['FileName'] == fileName)]
    X = subDF['Times'].values[0].split(',')
    Y = subDF['Intensities'].values[0].split(',')
    
    if smooth > 0:
        Y = vizLib.smoothListGaussian([float(i) for i in Y], degree=smooth)
    ax.plot(X, Y, color=color)
    if zoom:
        if numpy.max(Y) > 100:
            maxIndex = numpy.argmax(numpy.array(Y))
            ax.set_xlim(float(X[maxIndex])-2, float(X[maxIndex])+2)

    return ax
    
def plotMRM_index(dataFrame, index, ax, color='grey', smooth=0, zoom=False):
    [fn, pms, preC, fi, proC, isotope] = getInfoMRM(dataFrame, index)
    plotMRM(dataFrame, fn, pms, preC, fi, proC, isotope, ax, color=color, smooth=smooth, zoom=zoom)
    return ax

def plotLHMRM_index(dataFrame, index, ax, smooth=0):
    [fn, pms, preC, fi, proC, isotope] = getInfoMRM(dataFrame, index)
    indexOpp = getPairedIonMRM(dataFrame, index)
    if isotope == 'heavy':
        [index, indexOpp] = [indexOpp, index]
    ax = plotMRM_index(dataFrame, index, ax, color='blue', smooth=smooth)
    ax = plotMRM_index(dataFrame, indexOpp, ax, color='red', smooth=smooth)
    return ax

def plotLHMRM(dataFrame, fileName, pepSeq, precursorCharge, prodCharge, fragIon, ax, smooth=0):
    ax = plotMRM(dataFrame, fileName, pepSeq, precursorCharge, prodCharge, fragIon, 'light', ax, color='blue', smooth=smooth)
    ax = plotMRM(dataFrame, fileName, pepSeq, precursorCharge, prodCharge, fragIon, 'heavy', ax, color='red', smooth=smooth)
    return ax

def prettyPlot3TransMRM(dataFrame, toPlotIndex, figsize=(33,8.5)):
    axisArray = []
    pylab.figure(figsize=figsize)
    for i in range(len(toPlotIndex)):
        a = plt.subplot2grid((2,9), (0, i*3), colspan=1, rowspan=2)
        a = plotTotalChromPairsMRM(dataFrame, toPlotIndex[i], a, smooth=5, zoom=True)
        info = getInfoMRM(dataFrame, toPlotIndex[i])
        a.set_title(info[1])
        vizLib.tickNum(a, xAxis=4, yAxis=4)
        vizLib.cleanAxis(a, ticksOnly=True)
        axisArray.append(a)
    
        b = plt.subplot2grid((2,9), (0,i*3+1), colspan=2, rowspan=1)
        b = plotAllTransitionsMRM(dataFrame, toPlotIndex[i], b, smooth=5, zoom=True)
        info = getInfoMRM(dataFrame, toPlotIndex[i])
        b.set_title('light')
        vizLib.tickNum(b, xAxis=4, yAxis=4)
        vizLib.cleanAxis(b, ticksOnly=True)
        axisArray.append(b)
    
        c = plt.subplot2grid((2,9), (1,i*3+1), colspan=2, rowspan=1)
        opp = getPairedIonMRM(dataFrame, toPlotIndex[i])
        c = plotAllTransitionsMRM(dataFrame, opp, c, smooth=5, zoom=True)
        info = getInfoMRM(dataFrame, toPlotIndex[i])
        c.set_title('heavy')
        vizLib.tickNum(c, xAxis=4, yAxis=4)
        vizLib.cleanAxis(c, ticksOnly=True)
        axisArray.append(c)
    pylab.tight_layout()
    return axisArray