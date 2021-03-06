# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 21:29:59 2013

@author: jhdavis
"""

'''DEPRECATED FUNCTIONS'''
def _calculateStatsFile(dataFrame, protein_set, numerator, denominator, normalization=1.0):
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

def _getProtMedValue(statsFileDict, proteinName):
    """getProtMedValue gets the median value of proteins values from a StatsFile-style data structure

    :param statsFileDict: a statsFile-style dictionary (can be generated by calculateStatsFile 
        (which takes a data_dictionary from readCSV))
    :type statsFileDict: dictionary (statsFilte type (needs keys of protein names (strings) and second dict of 'vals'))
    :param proteinName: name of the protein to get the med from
    :type proteinName: string
    :returns:  the median value in the 'vals' field for that protein
    
    """
    return numpy.median(numpy.array(statsFileDict[proteinName]['vals']))

def _mergeFiles(fileList, pulse, numerator, denominator, proteinToNormalizeTo=None):
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

def _getInfoToPlotStats(datapath, pulse, numerator, denominator, proteinToNormalizeTo=None):
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

def _composeDataSets(comp, toAdd):
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
    
def _readCSV(datapath, pulse):
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

def _getMedian(dataByProtein, proteinName):
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

def _calculateStatsFile(data_dictionary, protein_set, numerator, denominator, normalization=1.0):
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

def _outputStatsFile(data_dictionary, protein_set, filePrefix, numerator, denominator):
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


def _calcValue(data, num, den, offset=0.0):
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
