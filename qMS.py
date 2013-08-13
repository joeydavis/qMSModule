"""
.. module:: qMS 
   :platform: Any
   :synopsis: A collection of functions for qMS data processing and analysis

.. moduleauthor:: Joey Davis <joeydavis@gmail.com>; Josh Silverman <josh.silverman@gmail.com>

"""

import os
import csv
import sys
import urllib2
import numpy
import qMSDefs
import pandas as pd
import re

#qMS utilities

path = os.getcwd()
sys.setrecursionlimit(10000000)
[startindex, pep_seq_index, exp_mr_index, csvopen] = [0, 0, 0, ''];
referencePath = '/home/jhdavis/scripts/python/modules/qMSmodule/'

#######Nice sorting and printing utilities######
def tryint(s):
    try:
        return int(s)
    except:
        return s
     
def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)
    return l


def readIsoCSV(filename, columns=None):
    """readIsoCSV takes a filename pointing to a _iso.csv file. It returns the calculated
        pandas dataFrame. Optional argument columns can be used to specify specific column headers

    :param filename: full path to the _iso.csv file
    :type filename: string
    :param columns: optional list of strings with the columns to incorporate into the dataFrame
    :type columns: list of strings
    :returns:  a pandas dataFrame with the relevant contents of the _iso.csv. Function
        automatically determines if the dataset is a pulse (bears AMP_S) or variable labeling (bears FRC_NX)

    """

    if columns is None:
        columns=['isofile', 'isoprotein', 'isopep', 'mw', 'isoz_charge', 'tim', 
                 'chisq', 'symb', 'mz', 'B', 'OFF', 'GW', 'AMP_U', 'AMP_L', 
                 'rt_n14', 'rt_n15', 'mz_n14', 'mz_n15',
                 'ppm_n14', 'ppm_n15', 'n14mass', 'n15mass', 'protein', 'startres',
                 'endres', 'charge', 'missed', 'seq', 'mod', 'seqmod', 'file', 'currentCalc', 'resid', 'minIntensity', 'handDelete', 'handSave']
        r = csv.reader(open(filename))
        header = r.next()
        pulse = 'AMP_S' in header
        varLab = 'FRC_NX' in header
        if pulse:
            columns.append('AMP_S')
        if varLab:
            columns.append('FRC_NX')
    data = pd.read_csv(filename, usecols=columns)
    data = data[pd.notnull(data['AMP_U'])]
    data = data[pd.notnull(data['AMP_L'])]
    if not pulse:
        data = data.rename(columns={'AMP_L': 'AMP_S'})

    positionOtherDict = {key:int(value)+1 for value, key in enumerate(sort_nicely(sorted(set(data['protein'].values))))}
    positionLookupOther = pd.Series(positionOtherDict)
    data['70Spos']=qMSDefs.positionLookup70S[data['protein']].values
    data['50Spos']=qMSDefs.positionLookup50S[data['protein']].values
    data['30Spos']=qMSDefs.positionLookup30S[data['protein']].values
    data['otherpos']=positionLookupOther[data['protein']].values
    
    data['currentPos']=data['otherpos']
    data['ppmDiff']=data['ppm_n14'] - data['ppm_n15']
    data['rtDiff']=data['rt_n14'] - data['rt_n15']
    data['handDelete'] = False
    data['handSave'] = False
    return data

def calcStatsDict(dataFrame, numerator, denominator, normalization=1.0, offset=0.0):
    """calcStatsDict takes a dataFrame and , a numerator, a denominator, an offset (applied to value first)
        and a normalization factor (scaling factor applied last). It returns a dictionary with keys of 
        protein names and values as a numpy array of calculated values based on numerator and denominator keys.

    :param dataFrame: a pandas dataFrame as calculated from an _iso.csv file (qMS.readIsoCSV)
    :type dataFrame: pandas dataFrame
    :param numerator: a list of strings identifying the numerator (must be AMP_U, AMP_L, and/or AMP_S)
    :type numerator: list of strings
    :param denominator: a list of strings identifying the numerator (must be AMP_U, AMP_L, and/or AMP_S)
    :type denominator: list of strings
    :param normalization: a float normalization factor if you want to scale all of the values uniformly (applied last)
    :type normalization: float
    :param offset: a float offset factor if you want to offset the values (applied first)
    :type offset: float
    :returns:  a dictionary of numpy arrays. First key is the protein name (one of those given in protein_set).
        This leads to a numpy array with the list of values

    """
    ps = list(set(dataFrame['protein'].values))
    ps.sort()
    return {p:calcValue(dataFrame[dataFrame['protein']==p], numerator, denominator, offset=offset).values*normalization for p in ps}

def multiStatsDict(isoFileList, num, den, normalization=1.0, offset=0.0, normProtein=None):
    """multiStatsDict takes a list of _iso.csv files and a list of nums and dens.
        It returns a dict of dict (first key is the file name) - this leads to a statsDict
        that is keyed by protein names. All of the statsDicts contain a full compliment of keys 
        (based on the file list), with empty numpy arrays if there were no values in the original
        dataset

    :param isoFileList: a list of file paths (full paths with entire file name)
    :type isoFileList: a list of strings
    :param num: a list of strings identifying the numerator (must be AMP_U, AMP_L, and/or AMP_S)
    :type num: list of strings
    :param den: a list of strings identifying the numerator (must be AMP_U, AMP_L, and/or AMP_S)
    :type den: list of strings
    :param normalization: a float normalization factor if you want to scale all of the values uniformly
    :type normalization: float
    :param offset: a float offset factor if you want to alter the values uniformly
    :type offset: float
    :param normProtein: string of the protein to normalize to (will be to the median)
    :type normProtein: string
    :returns:  a dictionary of of dicationaries of numpy arrays. First key is the file name, this leads
        to a dictionary where the first key is the protein name (one of those given in protein_set).
        This leads to a numpy array with the list of values

    """
    allPs = calcStatsDict(readIsoCSV(isoFileList[0]), num, den)
    fileStatsDict = dict()
    for f in isoFileList:
        fileStatsDict[f] = calcStatsDict(readIsoCSV(f), num, den, normalization=normalization, offset=offset)
        if not normProtein is None:
            normValue = 1/numpy.median(fileStatsDict[f][normProtein])
            fileStatsDict[f] = calcStatsDict(readIsoCSV(f), num, den, normalization=normValue, offset=offset)
        
    for f in isoFileList[1:]:
        allPs = appendKeys(allPs, fileStatsDict[f])
    
    for f in isoFileList:
        fileStatsDict[f] = appendKeys(fileStatsDict[f], allPs)
    
    return fileStatsDict

def mergeFiles(fileList, numerator, denominator, normProtein=None):
    """mergeFiles takes a list of _iso.csv files and returns a merged statsFile data dictionary structure

    :param fileList: a list of strings with the full path for reach .csv file
    :type fileList: list (of strings)
    :param numerator: a list of strings for what elements in the numerator (ampu, ampl, amps)
    :type numerator: list of strings
    :param denominator: a list of strings for what elements in the numerator (ampu, ampl, amps)
    :type denominator: list of strings
    :param proteinToNormalizeTo: string of what protein to normalize to (defaults to None)
    :type proteinToNormalizeTo: strings
    :returns:  the median value in the 'vals' field for that protein
    
    """

    splitDict = multiStatsDict(fileList, numerator,denominator, normProtein=normProtein)
    mergedDict = splitDict[sorted(splitDict.keys())[0]]
    for k in sorted(splitDict.keys())[1:]:
        cd = splitDict[k]
        for p in sorted(cd.keys()):
            mergedDict[p] = numpy.concatenate((mergedDict[p],cd[p]))
    return mergedDict
    


def getProtMedValue(statsFileDict, proteinName):
    """getProtMedValue gets the median value of proteins values from a StatsFile-style dictionary

    :param statsFileDict: a statsFile-style dictionary (can be generated by calcStatsDict 
        (which takes a pandas dataFrame from readIsoCSV))
    :type statsFileDict: dictionary (statsFilte type (needs keys of protein names (strings) and values of numpy arrays))
    :param proteinName: name of the protein to get the med from
    :type proteinName: string
    :returns:  the median value
    
    """
    return numpy.median(statsFileDict[proteinName])

def calcPercent(f, sigfig=2):
    """calcPercent takes a floating point number, rounds it to the number of sigfigs (default 2) and
        returns a string of that number multiplied by 100 and with a ' %'

    :param f: a floating point number to be converted to percentage
    :type f: float
    :param sigfig: the number of sigfigs to round to (default 2)
    :type sigfig: int
    :returns:  a string of the float rounded + ' %'

    """
    s = str(round(f, sigfig)*100)
    return s+" %"

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
        Equation from Stephen Chen's paper. USE THIS TO FIT THE LABELING OF TERMINAL (70S) pool.

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

def poolInterFunc(k,t,P):
    """poolInterFunc is a function to calculate the labeling kinetics for a protein using the 
        overlabeling of an intermediate. Derived from the differential equation in Stephen Chen's paper

    :param k: the growth rate (calculated as ln(2)/doubling time)
    :type k: float
    :param t: the time (if a single number, you get back teh max lab at that point), 
        can accept an array and will give back the curve    
    :type t: float/int/array
    :param d: the turnover rate
    :type d: float
    :returns:  the expected labeling lab (array or single value)
    
    """
    return 1.0 - (numpy.exp((0.0-k)*(1.0+(1.0/P))*t))
    
def poolInterFracXFunc(k,t,P,X=0.52):
    """poolInterFunc is a function to calculate the labeling kinetics for a protein using the 
        overlabeling of an intermediate. Derived from the differential equation in Stephen Chen's paper

    :param k: the growth rate (calculated as ln(2)/doubling time)
    :type k: float
    :param t: the time (if a single number, you get back teh max lab at that point), 
        can accept an array and will give back the curve    
    :type t: float/int/array
    :param d: the turnover rate
    :type d: float
    :returns:  the expected labeling lab (array or single value)
    
    """
    return X*poolInterFunc(k,t,P) + (1.0-X)*poolFunc(k,t,P)

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
    
    return 1.0 - (numpy.exp(-(d+k)*t))

def growthRate(doublingTime):
    """growthRate calculates a growth rate k given a doubling time (units of k can be used in equations)
        listed in Stephen Chen's paper

    :param doublingTime: the doubling time of the cell in units of time (mins, secs)
    :type doublingTime: floatz
    :returns:  the growth rate k (units of inverse time) - calculated at ln(2)/k
    
    """
    
    return numpy.log(2)/float(doublingTime)

def addBlankKey(d, k):
    """addBlankKey takes a stats dict and a key. It checks if the key in in the stats file,
        if it is not, it make a new dict entry with the key k and the value as an empty numpy array

    :param d: statsDictionary
    :type d: a dict of numpy arrays - keyed by protein names
    :param k: a string name for the key to check
    :type k: string
    :returns:  statsDictionary (a dictionary of numpy arrays). First key is the protein name (one of those given in protein_set).
        This leads to a numpy array with the list of values

    """
    if not d.has_key(k):
        d[k]=numpy.array([])
    return d
    
def appendKeys(d1, d2):
    """appendKeys a helper function that adds all the keys in d2 to d1 (making empty entries into d1)

    :param d1: statsDictionary
    :type d1: a dict of numpy arrays - keyed by protein names
    :param d2: statsDictionary
    :type d2: a dict of numpy arrays - keyed by protein names
    :returns:  statsDictionary (a dictionary of numpy arrays). First key is the protein name (one of those given in protein_set).
        This leads to a numpy array with the list of values

    """

    for k in d2:
        d1 = addBlankKey(d1, k)
    return d1

def calcValue(df, num, den, offset=0.0):
    """calcVale takes a pandas dataFrame bearing the keys to be used and calculates
        the ratio of the num/den (specified as lists of the keys - AMP_U, AMP_L, AMP_S)

    :param df: pandas dataFrame with the information from an _iso.csv
    :type df: pandas dataFrame
    :param num: a list of the strings in the numerator (AMP_U, AMP_L, AMP_S)
    :type den: list of strings
    :param den: a list of strings identifying the denominator (AMP_U, AMP_L, AMP_S)
    :type den: list of strings
    :param offset: an optional offset float to move the point up or down
    :type offset: float
    :returns:  a float of the [(sum of the nums) divided by the (sum of the dens)] + the offset
    
    """
    nsDF = df[num[0]]
    dsDF = df[den[0]]
    for x in num[1:]:
        nsDF = nsDF + df[x]
    for x in den[1:]:
        dsDF = dsDF + df[x]
    value = nsDF/dsDF + offset
    return value

def boolParse(s):
    """boolParse takes a string and returns a bool. Any capitilization of "true" results in True
        all other strings result in False

    :param s: string to process
    :type s: string
    :returns:  bool (see description)
    
    """
    return s.upper()=='TRUE'

def preProcessIsoCSV(isoPath, genPlots):
    """preProcessIsoCSV processes an _iso.csv and _plots directory to append the 
        residual column as well as a "currentCalc" column showing protein levels.
        Function is a helper that calls calcResidual

    :param isoPath: full path to the _iso.csv file
    :type isoPath: string
    :param genPlots: bool indicating if .plots files should be generated
    :type genPlots: list of strings
    :returns: a complteted dataFrame with resids and currentCalc columns. Has
        externality of generating .plots files in the _peaks directory if genPlots is true
    
    """
    dataFrame = readIsoCSV(isoPath)
    dataFrame['currentCalc'] = calcValue(dataFrame, ['AMP_U'], ['AMP_U', 'AMP_S'])
    rootPath = '/'.join(isoPath.split('/')[:-1])+'/'
    dataFrame = calcResidual(rootPath, dataFrame, genPlots=genPlots)

    fileName = isoPath.split('/')[-1:][0].replace('.', '_res.')
    dataFrame.to_csv(rootPath+fileName, index=False)
    return dataFrame

def cleanPlotsDir(path, extensions=['.newcon', '.fit', '.png']):
    """cleanPlotsDir executes a system command to remove files from a _peaks directory

    :param path: full path to the _peaks directory (no trailing /)
    :type path: path
    :param extensions: list of strings with extension types
    :type extensions: list of strings
    :returns: no return, outputs the command issued and any output from the command.
        Has externalities of deleting all *.ext files in the path directory
    
    """
    for ext in extensions:
        command = 'rm -rf ' + path + '/*' + ext
        print command
        output = os.popen(command).read()
        print output

def calcResidual(datapath, dataFrame, genPlots=False):
    """calcResidual takes a path to a _peaks directory and a pandas dataFrame containing the contents
        of the _iso.csv file. It appends a column to the dataFrame with the calculated residual
        (the difference between the fit and data in a .dat file) for each peptide. It also appends a
        column to the dataFrame with the max fit intensity
        Optional paramter genPlots will also generate .plots files that can be used to plot the datasets

    :param datapath: A string with the full path to the _plots directory
    :type datapath: string
    :param dataFrame: A pandas dataframe with ['isofile'] at the minimum (must point to the .dat files in 
        the _plots directory)
    :type dataFrame: pandas dataframe
    :param genPlots: Optional boolean telling function if it should write .plots files
    :type genPlots: boolean

    :returns:  the dataFrame modified to include the 'resid' and 'minIntensity' columns. .Dats that cause any errors
        are given fits with constant 666, resids with constant 666, and minIntensity with constant -666.
    
    """
    dataFrame['resid']=0
    dataFrame['minIntensity']=0
    for iso in dataFrame['isofile'].values:
        datFileName = datapath+iso+'.dat'
        try:
            datPd = pd.read_csv(datFileName, names=['offset', 'dat', 'resid'], header=None)
            del datPd['offset']
            datPd['residAdj'] = datPd['resid']+(datPd['dat'].median()-datPd['resid'].median())
            datPd['fit'] = datPd['residAdj']+datPd['dat']
            calcResid = datPd['residAdj'].abs().sum()/min([datPd['fit'].max(), datPd['dat'].max()])
            calcMinIntensity = min([datPd['fit'].max(), datPd['dat'].max()])
        except (IOError, TypeError, Exception) as e:
            print "Error " + e.message + " in " + datFileName
            sys.stdout.flush()
            datPd = pd.DataFrame()
            datPd['fit'] = 666
            datPd['residAdj'] = 666
            calcResid = 666
            calcMinIntensity = -666
        row = dataFrame[dataFrame['isofile']==iso]
        row['resid'] = calcResid
        row['minIntensity'] = calcMinIntensity
        dataFrame.update(row)
        if genPlots:
            datPd.to_csv(datapath+iso+'.plots', index=False)
    return dataFrame

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

def getRPInfo(numberGenes, start=10, baseURL='http://ribosome.med.miyazaki-u.ac.jp/rpg.cgi?mode=gene&id=ECO100'):
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
    #gs = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
    genes = range(start, start+numberGenes)
    #genes = gs + genes
    base = baseURL
    addys = [base + str(i) for i in genes]
    print addys
    rpdictGN = {}
    rpdictGP = {}
    for address in addys:
        website = urllib2.urlopen(address)
        website_html = website.read()
        gni = website_html.find('>Gene Name</td><td>')
        gn = website_html[gni+19:gni+26]
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

def listReplace(l, to, rv):
    """listReplace is a helper function replaces all occurances of to with rv in a list l

    :param l: list to be replaced
    :type l: list
    :param to: item to be replaced
    :type to: string
    :param rv: item to replace with
    :type rv: string
    :returns: a list with all occurances of to replaced with rv
    
    """
    tr = []
    for i in l:
        if i == to:
            tr.append(rv)
        else:
            tr.append(i)
    return tr

def printSortedDict(d):
    """printSortedDict is a helper function to print a dictionary

    :param d: a dictionary to be printed
    :type d: dict
    :returns: a string of the dictionary
    
    """
    k = d.keys()
    k.sort()
    tp = ''
    for i in k:
        tp = tp + str(i) + ":" + str(d[str(i)]) + ", "
    return tp

def calcMW(seq):
    mw = 0
    for aa in seq:
        mw = mw+qMSDefs.aaweights[aa.upper()]
    return mw
