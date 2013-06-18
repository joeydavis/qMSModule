"""
.. module:: qMS 
   :platform: Any
   :synopsis: A collection of functions for qMS data processing and analysis

.. moduleauthor:: Josh Silverman <josh.silverman@gmail.com>; Joey Davis <joeydavis@gmail.com>

"""

import os
import csv
import sys
import urllib2
import numpy
import pylab
import qMSDefs
import pandas as pd

#qMS utilities

path = os.getcwd()
sys.setrecursionlimit(10000000)
[startindex, pep_seq_index, exp_mr_index, csvopen] = [0, 0, 0, ''];
referencePath = '/home/jhdavis/scripts/python/modules/qMSmodule/'

def readIsoCSV(filename, columns=None):
    if columns is None:
        columns=['isofile', 'isoprotein', 'isopep', 'mw', 'isoz_charge', 'tim', 
                 'chisq', 'symb', 'mz', 'B', 'OFF', 'GW', 'AMP_U', 'AMP_L', 
                 'rt_n14', 'rt_n15', 'mz_n14', 'mz_n15',
                 'ppm_n14', 'ppm_n15', 'n14mass', 'n15mass', 'protein', 'startres',
                 'endres', 'charge', 'missed', 'seq', 'mod', 'seqmod', 'file', 'currentCalc', 'resid']
        r = csv.reader(open(filename))
        header = r.next()
        pulse = 'AMP_S' in header
        varLab = 'FRC_NX' in header
        if pulse:
            columns.append('AMP_S')
        if varLab:
            columns.append('FRC_NX')
    data = pd.read_csv(filename, usecols=columns)
    if not pulse:
        data = data.rename(columns={'AMP_L': 'AMP_S'})

    data['70Spos']=qMSDefs.positionLookup70S[data['protein']].values
    data['50Spos']=qMSDefs.positionLookup50S[data['protein']].values
    data['30Spos']=qMSDefs.positionLookup30S[data['protein']].values
    data['currentPos']=data['70Spos']
    data['ppmDiff']=data['ppm_n14'] - data['ppm_n15']
    data['rtDiff']=data['rt_n14'] - data['rt_n15']
    #data['50Spos'] = qMSDefs.positionLookup50S[data['protein']]
    #data['30Spos'] = qMSDefs.positionLookup30S[data['protein']]
    return data
    

def calcPercent(f, sigfig=2):
    """calcPercent takes a floating point number, rounds it to the number of sigfigs (default 2) and
        returns a string of that number multiplied by 100 and with a ' %'

    :param f: a floating point number to be converted to percentage
    :type f: float
    :param sigfig: the number of sigfigs to round to (default 2)
    :type sigfig: int
    :returns:  a string of the float rounded + ' %'

    """
    return str(round(f, sigfig)*100)+" %"

def maxLabFunc(k,t):
    """maxLabFunc is a function to calculate the max labeling based on the equation from Stephen Chen's paper

    :param k: the growth rate (calculated as ln(2)/doubling time)
    :type k: float
    :param t: the time (if a single number, you get back teh max lab at that point), 
        can accept an array and will give back the curve
    :type t: float/int/array
    :returns:  the max lab (array or single value)
    
    """
    return 1 - numpy.exp(-k*t)

def poolFunc(k,t,P):
    """poolFunc is a function to calculate the labeling kinetics for a protein with a given pool size
        Equation from Stephen Chen's paper

    :param k: the growth rate (calculated as ln(2)/doubling time)
    :type k: float
    :param t: the time (if a single number, you get back teh max lab at that point), 
        can accept an array and will give back the curve    
    :type t: float/int/array
    :param P: The pool size (expressed as precursorPool/completedRibosomePool)
    :type t: float
    :returns:  the expected labeling (array or single value)
    
    """
    
    return 1.0 + (P*numpy.exp((0.0-k)*(1.0+(1.0/P))*t)) - ((1.0 + P)*numpy.exp((0.0-k)*t))

def overLabelingFunc(k,t,d):
    """overLabelingFunc is a function to calculate the labeling kinetics for a protein with a given turnover rate
        Equation from Stephen Chen's paper

    :param k: the growth rate (calculated as ln(2)/doubling time)
    :type k: float
    :param t: the time (if a single number, you get back teh max lab at that point), 
        can accept an array and will give back the curve    
    :type t: float/int/array
    :param d: the turnover rate
    :type d: float
    :returns:  the expected labeling lab (array or single value)
    
    """
    
    return 1.0 - (numpy.exp(((0.0-k)+d)*t))

def growthRate(doublingTime):
    """growthRate calculates a growth rate k given a doubling time (units of k can be used in equations)
        listed in Stephen Chen's paper

    :param doublingTime: the doubling time of the cell in units of time (mins, secs)
    :type doublingTime: floatz
    :returns:  the growth rate k (units of inverse time) - calculated at ln(2)/k
    
    """
    
    return numpy.log(2)/float(doublingTime)

def getMedian(dataByProtein, proteinName):
    """getMediantakes the outputs of calculateStatsFile (a dataByProteinDictionary structure), and a protein name,
        it returns the median value of that protein

    :param dataByProtein: a dictionary datastructure bearing keys of protein names, with values of dictionaries bearing
        keys of vals, nvals, loc, flab
    :type dataByProtein: a calculateStatsFiles dictionary of dictionaries
    :param proteinName: a string of the protein to get the median value of 
    :type proteinName: string
    :returns:  a float of the median value in the datastructure for that protein

    """
    return float(numpy.median(dataByProtein[proteinName]["vals"]))

def calculateStatsFileDF(dataFrame, protein_set, numerator, denominator, normalization=1.0):
    """calculateStatsFile takes a dataFrame and a protein_set, a numerator, a denominator,
        and a normalization factor. It returns the information necessary for a stats file (vals, nvals, loc, flab)

    :param data_dictionary: a dictionary of dictionaries, must contain uid keys for each peptide fit
        each element in the dict is a dict with the keys protein, ampu, ampl, amps
    :type filename: dictionary or dictionaries
    :param protein_set: a set of protein names (must be the same names as given in the key protein in the subdict listed above)
    :type protein_set: list
    :param numerator: a list of strings identifying the numerator (must be ampu, ampl, and/or amps)
    :type numerator: list of strings
    :param denominator: a list of strings identifying the numerator (must be ampu, ampl, and/or amps)
    :type denominator: list of strings
    :param normalization: a float normalization factor if you want to scale all of the values uniformly
    :type pnormalization: float
    :returns:  a dictionary of dictionaries. First key is the protein name (one of those given in protein_set).
        This leads to a dictionary bearing keys "vals", "nvals", "loc", "flab" which are a list of values,
        the size of that list, the std dev of that list, the mean of that list

    """
    dataByProtein = {}
    for p in protein_set:
        dataByProtein[p] = {"flab" : -1, "loc" : -1, "nvals" : 0, "vals" : []}
        
    for uid in data_dictionary:
        parentProtein = data_dictionary[uid]["protein"]

        valNum = 0
        for i in numerator:
            valNum = valNum + data_dictionary[uid][i]
        valDen = 0
        for i in denominator:
            valDen = valDen + data_dictionary[uid][i]

        value = float(valNum)/float(valDen)/normalization
        dataByProtein[parentProtein]["vals"].append(value)
        dataByProtein[parentProtein]["nvals"] = dataByProtein[parentProtein]["nvals"] + 1
        
    for p in protein_set:
        dataByProtein[p]["loc"] = numpy.std(dataByProtein[p]["vals"])
        dataByProtein[p]["flab"] = numpy.mean(dataByProtein[p]["vals"])
        
    return dataByProtein

def calculateStatsFile(data_dictionary, protein_set, numerator, denominator, normalization=1.0):
    """calculateStatsFile takes the outputs of readCSV (data_dictionary, protein_set), a numerator, a denominator,
        and a normalization factor. It returns the information necessary for a stats file (vals, nvals, loc, flab)

    :param data_dictionary: a dictionary of dictionaries, must contain uid keys for each peptide fit
        each element in the dict is a dict with the keys protein, ampu, ampl, amps
    :type filename: dictionary or dictionaries
    :param protein_set: a set of protein names (must be the same names as given in the key protein in the subdict listed above)
    :type protein_set: list
    :param numerator: a list of strings identifying the numerator (must be ampu, ampl, and/or amps)
    :type numerator: list of strings
    :param denominator: a list of strings identifying the numerator (must be ampu, ampl, and/or amps)
    :type denominator: list of strings
    :param normalization: a float normalization factor if you want to scale all of the values uniformly
    :type pnormalization: float
    :returns:  a dictionary of dictionaries. First key is the protein name (one of those given in protein_set).
        This leads to a dictionary bearing keys "vals", "nvals", "loc", "flab" which are a list of values,
        the size of that list, the std dev of that list, the mean of that list

    """
    dataByProtein = {}
    for p in protein_set:
        dataByProtein[p] = {"flab" : -1, "loc" : -1, "nvals" : 0, "vals" : []}
        
    for uid in data_dictionary:
        parentProtein = data_dictionary[uid]["protein"]

        valNum = 0
        for i in numerator:
            valNum = valNum + data_dictionary[uid][i]
        valDen = 0
        for i in denominator:
            valDen = valDen + data_dictionary[uid][i]

        value = float(valNum)/float(valDen)/normalization
        dataByProtein[parentProtein]["vals"].append(value)
        dataByProtein[parentProtein]["nvals"] = dataByProtein[parentProtein]["nvals"] + 1
        
    for p in protein_set:
        dataByProtein[p]["loc"] = numpy.std(dataByProtein[p]["vals"])
        dataByProtein[p]["flab"] = numpy.mean(dataByProtein[p]["vals"])
        
    return dataByProtein

def outputStatsFile(data_dictionary, protein_set, filePrefix, numerator, denominator):
    """outputStatsFile takes the outputs of readCSV (data_dictionary, protein_set), a numerator, a denominator,
        and a normalization factor. It returns the information necessary for a stats file (vals, nvals, loc, flab)

    :param data_dictionary: a dictionary of dictionaries, must contain uid keys for each peptide fit
        each element in the dict is a dict with the keys protein, ampu, ampl, amps
    :type filename: dictionary or dictionaries
    :param protein_set: a set of protein names (must be the same names as given in the key protein in the subdict listed above)
    :type protein_set: list
    :param numerator: a list of strings identifying the numerator (must be ampu, ampl, and/or amps)
    :type numerator: list of strings
    :param denominator: a list of strings identifying the numerator (must be ampu, ampl, and/or amps)
    :type denominator: list of strings
    :param normalization: a float normalization factor if you want to scale all of the values uniformly
    :type normalization: float
    :returns:  a dictionary of dictionaries. First key is the protein name (one of those given in protein_set).
        This leads to a dictionary bearing keys "vals", "nvals", "loc", "flab" which are a list of values,
        the size of that list, the std dev of that list, the mean of that list
   
   """
    
    fileSuffix = "-"
    for i in numerator:
        fileSuffix = fileSuffix+str(i)+"_p_"
    fileSuffix = fileSuffix[:-3]+"_o_"
    for i in denominator:
        fileSuffix = fileSuffix+str(i)+"_p_"
    fileName = filePrefix+fileSuffix[:-3]+".stats"
    outfile = open(fileName, 'w')
    
    outfile.write("Protein flab +\\- loc nval vals\n")

    dataByProtein = {}
    for p in protein_set:
        dataByProtein[p] = {"flab" : -1, "loc" : -1, "nvals" : 0, "vals" : []}
        
    for uid in data_dictionary:
        parentProtein = data_dictionary[uid]['protein']
        value = calcValue(data_dictionary[uid], numerator, denominator)
        dataByProtein[parentProtein]["vals"].append(value)
        dataByProtein[parentProtein]["nvals"] = dataByProtein[parentProtein]["nvals"] + 1
            
    ordered = list(protein_set)
    ordered.sort()
    locIter = 1
    for p in ordered:
        outString = str(p) + " " + str(numpy.mean(dataByProtein[p]["vals"]))[:6] + " " + \
                    str(numpy.std(dataByProtein[p]["vals"]))[:6] + " " + str(locIter) + " " + \
                    str(dataByProtein[p]["nvals"]) + " "
        for v in dataByProtein[p]["vals"]:
            outString = outString + str(v)[:6] + " "    
        outString = outString[:-1] + "\n"
        outfile.write(outString)
    outfile.close()

def calcValue(data, num, den, offset=0.0):
    """calcVale takes a data dictionary bearing the keys ampu, ampl, and (optionally) amps and calcultes
        the ratio of the num/den (specified as lists of the ampu, ampl, amps above)

    :param data: a dictionary with minimal keys for the ampu, ampl, amps listed in num, den
    :type data: a dict
    :param num: a list of the strings in the numerator (ampu, ampl, amps)
    :type den: list
    :param den: a list of strings identifying the numerator (must be ampu, ampl, and/or amps)
    :type den: list of strings
    :param offset: an optional offset float to move the point up or down
    :type offset: float
    :returns:  a float of the [(sum of the nums) divided by the (sum of the dens)] + the offset
    
    """
    ns = [data[x] for x in num]
    valNum = reduce(lambda x,y:x+y, ns)
    ds = [data[x] for x in den]
    valDen = reduce(lambda x,y:x+y, ds)
    return float(valNum)/float(valDen) + offset

def calcValueDF(df, num, den, offset=0.0):
    nsDF = df[num[0]]
    dsDF = df[den[0]]
    for x in num[1:]:
        nsDF = nsDF + df[x]
    for x in den[1:]:
        dsDF = dsDF + df[x]
    return nsDF/dsDF + offset

def boolParse(s):
    return s.upper()=='TRUE'
    
def preProcessIsoCSV(isoPath, genPlots):
    dataFrame = readIsoCSV(isoPath)
    dataFrame['currentCalc'] = calcValueDF(dataFrame, ['AMP_U'], ['AMP_U', 'AMP_S'])
    rootPath = '/'.join(isoPath.split('/')[:-1])+'/'
    dataFrame = calcResidual(rootPath, dataFrame, genPlots=genPlots)

    fileName = isoPath.split('/')[-1:][0].replace('.', '_res.')
    dataFrame.to_csv(rootPath+fileName, index=False)
    return dataFrame

def calcResidual(datapath, dataFrame, genPlots=False):
    dataFrame['resid']=0
    for iso in dataFrame['isofile'].values:
        datFileName = datapath+iso+'.dat'
        try:
            datPd = pd.read_csv(datFileName, names=['offset', 'dat', 'resid'], header=None)
            del datPd['offset']
            datPd['residAdj'] = datPd['resid']+(datPd['dat'].median()-datPd['resid'].median())
            datPd['fit'] = datPd['residAdj']+datPd['dat']
            calcResid = datPd['resid'].abs().sum()/min([datPd['fit'].max(), datPd['dat'].max()])
        except IOError:
            print "Had trouble reading the .dat file for " + datFileName
            datPd['fit'] = numpy.nan
            calcResid = numpy.nan
        row = dataFrame[dataFrame['isofile']==iso]
        row['resid'] = calcResid
        dataFrame.update(row)
        if genPlots:
            datPd.to_csv(datapath+iso+'.plots', index=False)
    return dataFrame

def getProtMedValue(statsFileDict, proteinName):
    """getProtMedValue gets the median value of proteins values from a StatsFile-style data structure

    :param statsFileDict: a statsFile-style dictionary (can be generated by calculateStatsFile 
        (which takes a data_dictionary from readCSV))
    :type statsFileDict: dictionary (statsFilte type (needs keys of protein names (strings) and second dict of 'vals'))
    :param proteinName: name of the protein to get the med from
    :type proteinName: string
    :returns:  the median value in the 'vals' field for that protein
    
    """
    return numpy.median(numpy.array(statsFileDict[proteinName]['vals']))

def mergeFiles(fileList, pulse, numerator, denominator, proteinToNormalizeTo=None):
    """mergeFiles takes a list of csv files and returns a merged statsFile data dictionary structure

    :param fileList: a list of strings with the full path for reach .csv file
    :type fileList: list (of strings)
    :param pulse: a boolean of whether the datasetes are pulsed
    :type pulse: boolean
    :param numerator: a list of strings for what elements in the numerator (ampu, ampl, amps)
    :type numerator: list of strings
    :param denominator: a list of strings for what elements in the numerator (ampu, ampl, amps)
    :type denominator: list of strings
    :param proteinToNormalizeTo: string of what protein to normalize to (defaults to None)
    :type proteinToNormalizeTo: strings
    :returns:  the median value in the 'vals' field for that protein
    
    """
    normalizedDataByProtein = getInfoToPlotStats(fileList[0], pulse, numerator, denominator, proteinToNormalizeTo)
    compositeDataSet = normalizedDataByProtein
    for datapath in fileList[1:]:
        normalizedDataByProtein = getInfoToPlotStats(datapath, pulse, numerator, denominator, proteinToNormalizeTo)
        compositeDataSet = composeDataSets(compositeDataSet, normalizedDataByProtein)
    return compositeDataSet

def makePlotWithDataSets(listOfDataSets, proteinsToPlot, names, colors=None, yMax=1.5):
    """makePlotWithDataSets is a helper function to make nice "massage-like" plots of sets of datasets

    :param listOfDataSets: a list of datasets to be plotted (each list must be a stats file-like data dictionary)
    :type listOfDataSets: a list of stats file -like datasets
    :param proteinsToPlot: a list of the proteins to be plotted (list of strings that are keys to the datasets)
    :type proteinsToPlot: list of strings
    :param name: the names of the datasets
    :type name: list of strings for the dataset names (same order as listOfDataSets)
    :param yMax: float of maximum y value
    :type yMax: float
    :returns: an axis with the plotted data
    
    """
    if colors is None:
        colors = ['#ae2221', '#d72c2b', '#e78180', '#25557d', '#3170a4', '#5696cc']
    set1 = listOfDataSets[0]
    offsets = float(len(listOfDataSets)+1)
    ax = plotStatsDataStruct(proteinsToPlot, set1, names[0], offset=1.0/offsets, markerSize = 10/float(len(listOfDataSets))+4, yMax=yMax)
    i = 1
    for dataSet in listOfDataSets[1:]:
        ax = addStatsDataStructToPlot(proteinsToPlot, dataSet, ax, names[i], offset=(1.0/offsets)*(i+1.25), markerSize = 10/float(len(listOfDataSets))+4, color=colors[i])
        i = i +1
    return ax

def getInfoToPlotStats(datapath, pulse, numerator, denominator, proteinToNormalizeTo=None):
    """getInfoToPlotStats is a helper funciton that reads a csv file and returns a normalized
        statsfiledictionary

    :param datapath: the full path to the .csv file
    :type datapath: string
    :param pulse: a bool if the dataset has been pulsed (does the csv have all of ampu, ampl, amps or just ampu, ampl)
    :type pulse: boolean
    :param numerator: a list of the items in the numerator (can be ampu, ampl, amps or any mix)
    :type numerator: list of strings
    :param denominator: a list of the items in the numerator (can be ampu, ampl, amps or any mix)
    :type denominator: list of strings
    :param proteinToNormalizeTo: the protein to normalize to (defaults to None)
    :type proteinToNormalizeTo: string
    :returns: a stats dictionary data structure that can then be easily plotted/manipulated
    
    """
    
    [data_dictionary, protein_set, peptide_list] = readCSV(datapath, pulse)
    dataByProtein = calculateStatsFile(data_dictionary, protein_set, numerator, denominator, normalization=1.0)
    if proteinToNormalizeTo is None:
        return dataByProtein
    else:
        normMed = getProtMedValue(dataByProtein, proteinToNormalizeTo)
        normalizedDataByProtein = calculateStatsFile(data_dictionary, protein_set, numerator, denominator, normalization=normMed)
        return normalizedDataByProtein

def composeDataSets(comp, toAdd):
    """composeDataSets takes two StatsFileDatasets and composes them to return one dataset

    :param comp: the StatsFileDatasturct to add to
    :type comp: statsFile dictionary of dictionaries (needs the 'vals' keys at least)
    :param toAdd: the StatsFileDatasturct to add 
    :type toAdd: statsFile dictionary of dictionaries (needs the 'vals' keys at least)
    :returns: the composed statsfileDataStruct
    
    """
    proteinsToAdd = toAdd.keys()
    for i in proteinsToAdd:
        if i in comp.keys():
            comp[i]["vals"] = comp[i]["vals"] + toAdd[i]["vals"]
            comp[i]["nvals"] = len(comp[i]["vals"])
            comp[i]["loc"] = numpy.std(comp[i]["vals"])
            comp[i]["flab"] = numpy.mean(comp[i]["vals"])
        else:
            comp[i] = toAdd[i]
    return comp

def plotStatsDataStruct(proteins, dataByProtein, name, offset=0.0, markerSize=12, color='#e31a1c', yMax = 1.5):
    """plotStatsDataStruct plots the contents of a stats file-like datastructure. proteins to be plotted are 
        listed in the non-redundent list, proteins. The data is in dataByProtein, the name is in in name.
        Decent colors are red (['#ae2221', '#d72c2b', '#e78180']) and blue (['#25557d', '#3170a4', '#5696cc'])

    :param proteins: a non-redudent list of the protein names to use (should be the IDs in dataByProtein)
    :type proteins: list
    :param dataByProtein: a dictionary (easily created by calculateStatsFile), must contain 'vals' for each p in proteins
    :type dataByProtein: dictionary
    :param name: the name of the dataset
    :type name: string
    :param offset: where to center the point (x axis), scales 0 (left edge; default) to 1.0 (right edge)
    :type offset: float
    :param markerSize: size of the dots (default = 12)
    :type markerSize: int
    :param color: color for the dataset (default #e31a1c)
    :type color: color
    :param yMax: the max value for the y axis
    :type yMax: float
    :returns:  a pyplot axis with the data plotted
        
    """

    xAxis = range(1,len(proteins)+1)
    fig = pylab.figure(figsize=(22,5))
    ax = fig.add_subplot(111)
    xs = []
    ys = []
    for x in xAxis:
        p = proteins[x-1]
        if p in dataByProtein.keys():
            for v in dataByProtein[p]["vals"]:
                xs.append(x+offset)
                ys.append(v)
    pylab.grid(b=True, which='major', color='grey', linestyle='--', axis='y', linewidth=1.5, alpha=0.5)
    pylab.grid(b=True, which='major', color='grey', linestyle='-', axis='x', linewidth=1.5, alpha=0.75)
    ax.plot(xs, ys, 'o', color=color, markersize=markerSize, label=name)
    pylab.xticks(xAxis, [item for item in proteins], rotation=45, size=15)
    pylab.xlim(1, len(proteins)+1)
    ####################################
    pylab.yticks([0,yMax/4, yMax/2, 3*yMax/4, yMax], size=15)
    ####################################
    pylab.ylim(0, yMax)
    return ax
    
def addStatsDataStructToPlot(proteins, dataByProtein, ax, name, offset=0.0, markerSize=12, color='#377db8'):
    """addStatsDataStructToPlot adds the contents of a stats file-like datastructure to an existing plot. proteins to be plotted are 
        listed in the non-redundent list, proteins. THESE MUST BE THE SAME PROTEINS GIVEN IN THE FIRST CALL TO PLOTSTATSDATASTRUCT
        The data is in dataByProtein, the axis to add to is in ax, the name is in in name.
        Decent colors are red (['#ae2221', '#d72c2b', '#e78180']) and blue (['#25557d', '#3170a4', '#5696cc'])

    :param proteins: a non-redudent list of the protein names to use (should be the IDs in dataByProtein)
    :type proteins: list
    :param dataByProtein: a dictionary (easily created by calculateStatsFile), must contain 'vals' for each p in proteins
    :type dataByProtein: dictionary
    :param ax: the pyplot axis to modify
    :type ax: pyplot axis
    :param name: the name of the dataset
    :type name: string
    :param offset: where to center the point (x axis), scales 0 (left edge; default) to 1.0 (right edge)
    :type offset: float
    :param markerSize: size of the dots (default = 12)
    :type markerSize: int
    :param color: color for the dataset (default #e31a1c)
    :type color: color
    :returns:  a pyplot axis with the data plotted
        
    """

    xAxis = range(1,len(proteins)+1)
    xs = []
    ys = []
    for x in xAxis:
        p = proteins[x-1]
        if p in dataByProtein.keys():
            for v in dataByProtein[p]["vals"]:
                xs.append(x+offset)
                ys.append(v)
    ax.plot(xs, ys, 'o', color=color, markersize=markerSize, label=name)
    return ax

def getRPSeqData(ID):
    """getRPSeqData is a helper function to get the AA and cDNA seqs for ribosomal proteins

    :param ID: string with the geneID
    :type ID: string
    :returns:  a list of strings, first is by cDNA, second by AA;
        
    """
    address = 'http://ribosome.med.miyazaki-u.ac.jp/rpg.cgi?id='+ID+'&mode=seq'
    website = urllib2.urlopen(address)
    website_html = website.read()
    cDNAi = website_html.find('>cDNA Sequence')
    cDNAj = website_html.find('</textarea></td></tr>\n<tr><td align="center" bgcolor="#FFFF80" width="150">Amino Acids Sequence</td><td>')
    cDNA = website_html[cDNAi+74:cDNAj-1]

    AAj = website_html.find('</textarea></td></tr>\n</table>\n<div class="footer">')
    AA = website_html[cDNAj+155:AAj]
    
    return [cDNA, AA] 

def getRPInfo(numberGenes, baseURL='http://ribosome.med.miyazaki-u.ac.jp/rpg.cgi?mode=gene&id=ECO100'):
    """getRPInfo generates two dictionaries with relevant ribosomal protein information from a database.
        It prefers the url listed above, but can be used with other organisms so long as the
        find commands still work

    :param numberGenes: the number of genes to look for
    :type numberGenes: int
    :param baseURL: a string pointing to the base url, defaults to the japanese database
    :type baseURL: string
    :returns:  a list of dictionaries, first is by geneNames, second by geneProduct;
        each dictionary is a dict of dicts with subkeys size, cDNA, AA and either GP or GN (opposite
        of the base key)
        
    """
    gs = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
    genes = range(10, numberGenes)
    genes = gs + genes
    base = baseURL
    addys = [base + str(i) for i in genes]
    rpdictGN = {}
    rpdictGP = {}
    for address in addys:
        website = urllib2.urlopen(address)
        website_html = website.read()
        gni = website_html.find('>Gene Name</td><td>')
        gn = website_html[gni+19:gni+23]
        gpi = website_html.find('ibosomal protein ')
        gp = website_html[gpi+17:gpi+20]
        gp = gp.upper()
        if gp[-1] == '<':
            gp = gp[:-1]
        if gp == 'ITL':
            gp = 'L9'
        if gp == 'L7/':
            gp = 'L7/L12'
        si = website_html.find('Gene Size [bp]</td><td>')
        s = website_html[si+23:si+27]
        if s[-1] == '<':
            s = s[:-1]
        cDNA, AA = getRPSeqData(address[-8:])
        rpdictGN[gn] = {'GP':gp, 'size':s, 'cDNA':cDNA, 'AA':AA}
        rpdictGP[gp] = {'GN':gn, 'size':s, 'cDNA':cDNA, 'AA':AA}
    return [rpdictGN, rpdictGP]
    

#fetches protein sequence from uniprot ID
def getsequence(uniprot):
	try:
		urlbase = "http://www.uniprot.org/uniprot/"
		req = urllib2.Request(urlbase + uniprot + ".fasta")
		response = urllib2.urlopen(req)
		lines = [x for i, x in enumerate(response) if i > 0]
		return "".join(map(lambda x: x.rstrip("\n"), lines))
	except urllib2.URLError:
		print "No internet connection or bad Uniprot id"
		
#Returns the MAD of a list
def MAD(list):
	temp_median = numpy.median(list)
	temp_absdevs = map(lambda item: abs(item - temp_median), list)
	return numpy.median(temp_absdevs)

'''DEPRECATED FUNCTIONS'''
def readCSV(datapath, pulse):
    """readCSV takes a filename (fully specified with datapath) and a boolean if this was a pulsed dataset
        it return s a list of datastructures encompassing relavent information from the .csv file

    :param datapath: a string specifiy the full path/name of the file to be read
    :type datapath: string
    :param pulse: a boolean indicating if the dataset is a pulsed dataset
    :type pulse: boolean
    :returns:  a list ([data_dictionary, protein_set, peptide_list])
        data_dictionary has keys as UIDs, and values as dicts with
            keys for protein, ampu, ampl, amps, isofile
        protein_set is a set (list) of each protein in the dataset
        peptide_list is a list of all peptides observed 
        
    """
    data = list( csv.reader( open(datapath, 'rU') ) )
    header = data[0]

# 1 - Find the indices for the quantities of interest using list comprehensions
    protein_index = [index for index, item in enumerate(header) if item == "protein"][0]
    ampu_index = [index for index, item in enumerate(header) if item == "AMP_U"][0]
    ampl_index = [index for index, item in enumerate(header) if item == "AMP_L"][0]
    UID_index = [index for index, item in enumerate(header) if item == "isofile"][0]
    if pulse:
        amps_index = [index for index, item in enumerate(header) if item == "AMP_S"][0]

# 2 - Declare the peptide list and the protein set
    peptide_list = []
    protein_set  = set()

    data_dictionary = {}
    csv_hold_dict = {}
    UID_list = []
    
# 3 - Loop over data set, collect amplitudes, charge state, peptide sequence, protein id into protein_set
    for line in data[1:]:
        protein = line[protein_index]
        ampu = float(line[ampu_index])
        ampl = float(line[ampl_index])
        if pulse:
            amps = float(line[amps_index])
        else:
            amps = float(0)
        UID = line[UID_index]
        csv_hold_dict[UID] = line
        UID_list.append(UID)
        identifier = {"protein": protein, "ampu": ampu, "ampl": ampl, "amps": amps, "uid": UID}
        data_dictionary[UID] = {"ampu": ampu, "ampl": ampl, "amps": amps, "protein": protein}
        protein_set.add(protein)
        peptide_list.append(identifier)
    return [data_dictionary, protein_set, peptide_list]