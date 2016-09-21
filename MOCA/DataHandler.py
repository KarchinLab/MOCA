'''
Get python modules
'''
import os
from random import shuffle
from itertools import chain
import cPickle
from collections import Counter
import csv

'''
Get third-party modules
'''
from numpy import mean, std, arange

'''
Get MOCA modules
'''
from Statistics import binary_matrix, z_matrix

class dot_dict(dict):
    '''
    For nested dictionaries, dot notation sometimes looks cleaner looking than 
    the standard python dictionary notation 
    '''
    def __getattr__(self, name):

        return self[name]
    
def get_labels(File, Row=0, Start=0):
    '''
    Gets desired "Row" (default=0). Initially written to get
    header (i.e., Row=0), but can get any row. "Start" (default=
    0) is the column to start at. Useful if the 0,0 coordinate
    of a matrix is not a column label, for instance. 
    '''

    if File[-4:] == ".csv":
        return filter(lambda x: x != "", [Line for Line in csv.reader(open(File,"rb"))][Row][Start:])
    elif File[-4:] == ".txt":
        return filter(lambda x: x != "", [Line for Line in csv.reader(open(File,"rU"), delimiter="\t")][Row][Start:])
    else:
        return open(File).readlines()[Row].split()[Start:]

def get_features(File, Column=0, Start=1):
    '''
    Gets desired 'Column' (default=0). Initially written to get
    row names (e.g., gene or drug names), but can get any column.
    'Start' (default=1) is the row to start at (defaults to '1' 
    because it is assumned there is a header (e.g., sample or 
    patient IDs (see get_labels above)).
    '''
    
    if File[-4:] == ".csv":
        return [Line[Column] for Line in csv.reader(open(File,"rb"))][Start:]
    elif File[-4:] == ".txt":
        return [Line[Column] for Line in csv.reader(open(File,"rU"), delimiter="\t")][Start:]
    else:
        return [Line.split()[Column] for Line in file(File)][Start:]

def get_variates(File, Row=1, Column=1, Type=str, NA="NA"):
    '''
    Returns the coordinates of an nxm matrix. Assumes there are header
    and row names (i.e., zeroth row and column are NOT part of the matrix
    and therefore they both default to '1'). 'Type' defaults to 'str',
    but can also be 'float' or 'int'. 
    '''

    if File[-4:] == ".csv":
        Matrix = [filter(lambda x: x != "", Line[Column:]) for Line in csv.reader(open(File,"rb"))][Row:]
    elif File[-4:] == ".txt":
        Matrix = [filter(lambda x: x != "", Line[Column:]) \
                  for Line in csv.reader(open(File,"rU"), delimiter="\t")][Row:]
    else:
        Matrix = [Line.split()[Column:] for Line in file(File)][Row:]

    return [[NA if Element == NA else Type(Element) for Element in Row] for Row in Matrix]

def get_ordered_vector(TemplateLabels, VectorLabels, Vector):
    '''
    Often we get multiple data types where the column labels 
    (patients or samples, for instance) are not in the same order 
    across data types. This function returns a vector ordered to 
    match a template. Additionally, the output will be reduced to the 
    intersection of VectorLabels with TemplateLabels. 
    '''
    
    return [Vector[VectorLabels.index(Label)] for Label in TemplateLabels]

def get_ordered_matrix(TemplateLabels, VectorLabels, Matrix):
    '''
    Often we get multiple data types where the column labels 
    (patients or samples, for instance) are not in the same order 
    across data types. This function returns a matrix ordered to 
    match a template vector. Additionally, the output will be reduced 
    to the intersection of VectorLabels with TemplateLabels. 
    '''
    
    return [get_ordered_vector(TemplateLabels, VectorLabels, Vector) \
                for Vector in Matrix]   

def get_path(DataType):
    '''
    For a particular data type (e.g., drug or expression) find the 
    file containing said data using moca paths file. Absolute paths
    can start with ${HOME} to indicate "home/user"
    '''

    Path = [Line.split()[2] for Line in open("Paths").readlines() \
            if Line.strip() and Line.split()[0] == DataType][0]
    if Path[:7] == "${HOME}":
        Path = Path.replace("${HOME}", os.environ["HOME"])
    
    return Path

def decompose_matrix(Data):
    '''
    Provided a path to a data matrix, return the Labels, Features, and Variates
    '''

    Data = get_path(Data)
    
    return get_labels(Data), get_features(Data), get_variates(Data)

def make_report(Data, Labels, Features, UntransformedFeatures, Variates, Arguments):
    '''
    Records information about the data being processed and any variables that were used during processing.
    These records are extracted from every matrix in the MOCA.data directory, and reported back to the
    user, if the 'Reports = Data' is set.
    '''

    Report = {}
    #Use the 'Entries' key to set the order of output when 'Reports = Results' is run.
    Report["Entries"] = ["Source data", "Label count", "Type      ", "Normalized", "Threshold",
                         "Feature count", "Missing data"]
    Report["Source data"] = get_path(Data)
    Report["Label count"] = len(Labels)
    Report["Feature count"] = len(Features)
    MissingData = Counter(chain(*Variates))[Arguments.NA]
    TotalData = float(len(list(chain(*Variates))))
    Report["Missing data"] = [MissingData, "(" + str(100*(MissingData/TotalData)) + "%)"]
    Report["Threshold"] = "NA"
    Report["Normalized"] = "NA"
    Report["Type      "] = "Binary"
    if Arguments.DataType.lower() == "categorical":
        Report["Type      "] = "Categorical"
        Report["Feature count"] = [len(UntransformedFeatures), "(" + str(len(Features)) + " total categories)"]
    
    if Arguments.DataType.lower() == "continuous":
        Report["Threshold"] = Arguments.Threshold
        Report["Normalized"] = Arguments.Normalize
        Report["Type      "] = "Continuous"
        Report["Feature count"] = [len(Features), "(" + str(len(UntransformedFeatures)) + " before transformation)"]

    return Report

def get_binary_matrix(Data, Features, Variates, Arguments):
    '''
    Transforms data to binary for MOCA. Missing values set by Arguments.NA stay the same. Zeros, which indicate
    a negative case remain zeros. All other values are assumed to indicate a positive case, and are thus 
    converted to int "1". 
    '''

    TransformedFeatures = []
    TransformedVariates = []

    for Feature in Features:
        if ":" not in Feature:
            TransformedFeatures.append(Feature + ":" + Data)
        else:
            TransformedFeatures.append(Feature)
        
        BinaryVector = [Arguments.NA if Element == Arguments.NA else (0 if Element in [0, "0"] else 1) \
                        for Element in Variates[Features.index(Feature)]]
        TransformedVariates.append(BinaryVector)

    return TransformedFeatures, TransformedVariates
    
def get_continuous_matrix(Data, Features, Variates, Arguments):
    '''
    What cutoffs should be used to binarize continuous-valued data? Provide the 'Threshold' argument with a 
    list of all desired cutoffs 
    '''
    
    TransformedFeatures = []
    TransformedVariates = []
    if not Arguments.Threshold[0]:
        print "You have data that needs to be binarized, but you haven't told MOCA how to 'Threshold' it!"
        print "Exiting..."
        exit()
    
    GreaterThanEqualTo = []
    LessThan = []
    for Threshold in Arguments.Threshold:
        if Threshold[:2] == ">=":
            GreaterThanEqualTo.append(float(Threshold[2:]))
        elif Threshold[0] == "<":
            LessThan.append(float(Threshold[1:]))
        else:
            print "MOCA error: Threshold argument called, but supplied cutoffs don't all begin with >= or <"
            print "exiting..."
            exit()
    
    for Threshold in GreaterThanEqualTo:
        for Feature in Features:
            BinaryVector = [Arguments.NA if Element == Arguments.NA else (1 if float(Element) >= Threshold else 0) \
                            for Element in Variates[Features.index(Feature)]]
            TransformedFeatures.append(Feature + ">=" + str(Threshold) + ":" + Data) 
            TransformedVariates.append(BinaryVector)

    for Threshold in LessThan:
        for Feature in Features:
            BinaryVector = [Arguments.NA if Element == Arguments.NA else (1 if float(Element) < Threshold else 0) \
                            for Element in Variates[Features.index(Feature)]]
            TransformedFeatures.append(Feature + "<" + str(Threshold) + ":" + Data) 
            TransformedVariates.append(BinaryVector)
                
    return TransformedFeatures, TransformedVariates

def get_categorical_matrix(Data, Features, Variates, Arguments):
    '''
    Converts categorical data to binary. For instance, if you had smoking history from 
    clinical data coded as 0 = "never smoker", 1 = "former smoker", and 2 = "current smoker":
    this function would convert that single feature into three binary features, where the three
    possible smoking statuses would each have there own vectors of 1s and 0s across patients. 
    You need to use the CategoricalData arguement, otherwise MOCA will just assume your categories
    represent a continuum (i.e., the data will be discretized by thresholding). 
    '''
    
    TransformedFeatures = []
    TransformedVariates = []
    
    for Feature in Features:
        Categories = list(set(Variates[Features.index(Feature)])) #For each features, figure out what the categories are
        Categories = filter(lambda Category: Category != Arguments.NA, Categories) #Missing values are not categories!
        for Category in Categories:
            CategoricalVariates = [Arguments.NA if Variate == Arguments.NA else (1 if Variate == Category else 0) \
                                       for Variate in Variates[Features.index(Feature)]]
            TransformedFeatures.append(Feature + "." + str(Category) + ":" + Data)
            TransformedVariates.append(CategoricalVariates)

    return TransformedFeatures, TransformedVariates

def process_data(Data, Arguments):
    '''
    Processes all input data and initiates the "pickle", which will be used and appended to for all MOCA
    calculations, analysis, and post-processing for a single project. This is a necessary first step to 
    any MOCA analysis. 
    '''

    if not Arguments.DataType or Arguments.DataType.lower() not in ["binary", "continuous", "categorical"]:
        print "If Mode = PreProcess, then your arguments file must specificy the 'DataType'."
        print "DataType options include: 'Binary', 'Continuous', or 'Categorical'."
        print "See the UsersManual.pdf if you don't understand what these datatypes mean."
        print "Exiting..."
        
    Labels, Features, Variates = decompose_matrix(Data)
        
    #Always store an untransformed version of the data
    FormattedFeatures = [] 
    FormattedVariates = []

    for Feature in Features:
        if ":" not in Feature:
            FormattedFeatures.append(Feature + ":" + Data)
        else:
            FormattedFeatures.append(Feature)
        FormattedVariates.append(Variates[Features.index(Feature)])
        
    #Now figure out which type of data we have (categorical, binary, continuous) and transform it
    if Arguments.DataType.lower() == "binary": #just change anything other than "0" or missing data to "1"
        TransformedFeatures, TransformedVariates = get_binary_matrix(Data, Features, Variates, Arguments)

    elif Arguments.DataType.lower() == "continuous":
        if Arguments.Normalize: Variates = z_matrix(Variates) #as of now, z-score normalization is all that's available
        TransformedFeatures, TransformedVariates = get_continuous_matrix(Data, Features, Variates, Arguments)
    
    elif Arguments.DataType.lower() == "categorical": #is the data categorical (e.g., stage of cancer)?
        TransformedFeatures, TransformedVariates = get_categorical_matrix(Data, Features, Variates, Arguments)

    else:
       print "If Mode = PreProcess, then your arguments file must specificy the 'DataType'."
       print "DataType options include: 'Binary', 'Continuous', or 'Categorical'."
       print "See the UsersManual.pdf if you don't understand what these datatypes mean."
       print "Exiting..."
       exit()
    
    return Labels, FormattedFeatures, FormattedVariates, TransformedFeatures, TransformedVariates

def pickle_matrices(Arguments):
    '''
    Saves the processed data in python pickles. Python pickles load very quickly relative to human-readable text. 
    This pickled data is required to do any real MOCA calculations. 
    '''

    for Data in Arguments.Data:
        DataDict = {}
        Labels, FormattedFeatures, FormattedVariates, TransformedFeatures, TransformedVariates \
            = process_data(Data, Arguments)
        DataDict["Labels"] = Labels
        DataDict["FormattedFeatures"] = FormattedFeatures
        DataDict["FormattedVariates"] = FormattedVariates
        DataDict["TransformedFeatures"] = TransformedFeatures
        DataDict["TransformedVariates"] = TransformedVariates
        DataDict["Report"] = make_report(Data, Labels, TransformedFeatures, FormattedFeatures, FormattedVariates, Arguments)

        if Arguments.Filename.lower() == "default":
            cPickle.dump(DataDict, open(get_path("MOCA.data") + "/" + Data, "wb"), -1)
        elif Arguments.Filename.lower() != "default" and len(Arguments.Data) > 1:
            cPickle.dump(DataDict, open(get_path("MOCA.data") + "/" + Arguments.Filename + "." + Data, "wb"), -1)
        else:
            cPickle.dump(DataDict, open(get_path("MOCA.data") + "/" + Arguments.Filename, "wb"), -1)

    return

def load_data(Arguments):
    '''
    This loads the processed and pickled input data, which is required for any MOCA calculation. 
    '''

    DataDict = dot_dict() #See dot_dict class to see why this is done
    FeatureDict = dot_dict()
    VariateDict = dot_dict()
    TransformedDict = dot_dict()
    TransformedFeatureDict = dot_dict()
    TransformedVariateDict = dot_dict()
    Report = dot_dict()

    #We only want labels common to all of the requested datatypes
    ReducedLabels = [set(cPickle.load(open(get_path("MOCA.data") + "/" + Data, "rb"))["Labels"]) \
              for Data in Arguments.Data]
    
    if Arguments.LabelList: ReducedLabels.append([Label.strip() for Label in file(get_path(Arguments.LabelList))])

    #if we're combining matrices, 'reduced' labels are all labels
    if Arguments.CombineData:
        ReducedLabels = list(chain(*ReducedLabels))
    else:
        ReducedLabels = list(set.intersection(*ReducedLabels))

    #Get the feature list if called
    if Arguments.FeatureList: FeatureList = [Feature.strip() for Feature in file(get_path(Arguments.FeatureList))]
    
    for Data in Arguments.Data:
        Report[Data] = cPickle.load(open(get_path("MOCA.data") + "/" + Data, "rb"))["Report"]
        
        Labels = cPickle.load(open(get_path("MOCA.data") + "/" + Data, "rb"))["Labels"]

        #First get the untransformed data
        Features = cPickle.load(open(get_path("MOCA.data") + "/" + Data, "rb"))["FormattedFeatures"]
        Variates = cPickle.load(open(get_path("MOCA.data") + "/" + Data, "rb"))["FormattedVariates"]
        Variates = get_ordered_matrix(ReducedLabels, Labels, Variates) #Reduce the variates too!
        
        FilteredFeatures = [] #Any feature not filtered by FeatureList...
        FilteredVariates = [] #...and the corresponding variates

        for Feature in Features:
            if Arguments.FeatureList:
                if Feature in FeatureList:
                    FilteredFeatures.append(Feature)
                    FilteredVariates.append(Variates[Features.index(Feature)])
            else:
                FilteredFeatures.append(Feature)
                FilteredVariates.append(Variates[Features.index(Feature)])  

        FeatureDict[Data] = FilteredFeatures
        VariateDict[Data] = FilteredVariates

        #Now deal with the previously transformed data
        Features = cPickle.load(open(get_path("MOCA.data") + "/" + Data, "rb"))["TransformedFeatures"]
        Variates = cPickle.load(open(get_path("MOCA.data") + "/" + Data, "rb"))["TransformedVariates"]
        Variates = get_ordered_matrix(ReducedLabels, Labels, Variates) #Reduce the variates too!
        FilteredFeatures = [] #Features that pass FeatureMin and not filtered by FeatureList
        FilteredVariates = [] #...and the corresponding variates

        for Feature in Features:
            if Arguments.FeatureList:
                if Feature in FeatureList and sum(filter(lambda Variate: Variate != Arguments.NA,
                                                         Variates[Features.index(Feature)])) >= Arguments.FeatureMin:
                    FilteredFeatures.append(Feature)
                    FilteredVariates.append(Variates[Features.index(Feature)])
            elif sum(filter(lambda Variate: Variate != Arguments.NA,
                Variates[Features.index(Feature)])) >= Arguments.FeatureMin:
                FilteredFeatures.append(Feature)
                FilteredVariates.append(Variates[Features.index(Feature)])
            
                
        TransformedFeatureDict[Data] = FilteredFeatures
        TransformedVariateDict[Data] = FilteredVariates

    DataDict["Report"] = Report

    DataDict["Labels"] = ReducedLabels
    DataDict["Features"] = FeatureDict
    DataDict["Variates"] = VariateDict

    TransformedDict["Features"] = TransformedFeatureDict
    TransformedDict["Variates"] = TransformedVariateDict
    DataDict["Transformed"] = TransformedDict

    return DataDict

def get_supervised_dataset(Arguments):
    '''
    This is the way data needs to look for most MOCA calculations (aside from unsupervised
    pairwise comparisons). Basically, we need to separate the phenotypes (and response features)
    from the markers (explanatory features). 
    '''

    Data = load_data(Arguments)

    Features = list(chain(*Data.Transformed.Features.values()))
    Variates = list(chain(*Data.Transformed.Variates.values()))

    Markers = []
    for Type in Arguments.Data:
        if Type == Arguments.Phenotype:
            Phenotypes = Data.Transformed.Features[Type]
        else:
            Markers.extend(Data.Transformed.Features[Type])

    return Data.Labels, Features, Variates, Phenotypes, Markers

def get_priors(Arguments):
    '''
    Read in a priors file. Priors are useful to focus the search, when making setworks, around features
    you think might be important. Make your own Priors file, or use MakePriors to make a priors file 
    from a previous Setworks run. 
    '''

    UnionPriors, IntersectionPriors, DifferencePriors = (), (), ()

    Priors = [Line.split() for Line in open(get_path("Priors")).readlines()]

    for Prior in Priors:
        if Prior[0].lower() == "union":
            UnionPriors = tuple(Prior[2:])
        if Prior[0].lower() == "intersection":
            IntersectionPriors = tuple(Prior[2:])
        if Prior[0].lower() == "difference":
            DifferencePriors = tuple(Prior[2:])
    
    return UnionPriors, IntersectionPriors, DifferencePriors

def combine_data(Arguments):
    '''
    '''

    AllData = {}
    for Data in Arguments.Data:
        ThisData = {}
        Labels, FormattedFeatures, FormattedVariates, TransformedFeatures, TransformedVariates \
            = process_data(Data, Arguments)

        ThisData["Labels"] = Labels
        ThisData["FormattedFeatures"] = [Feature[:Feature.index(":")] for Feature in FormattedFeatures]
        ThisData["FormattedVariates"] = FormattedVariates
        ThisData["TransformedFeatures"] = [Feature[:Feature.index(":")] for Feature in TransformedFeatures]
        ThisData["TransformedVariates"] = TransformedVariates

        AllData[Data] = ThisData

    Labels = list(chain(*[AllData[Data]["Labels"] for Data in Arguments.Data]))

    #First we get the features and their variates and join them
    TransformedFeatures, TransformedVariates = [], []   
    for Feature in set.intersection(*[set(AllData[Data]["TransformedFeatures"]) for Data in Arguments.Data]):
        TransformedFeatures.append(Feature + ":" + Arguments.Filename)
        TransformedVariates.append(list(chain(*[AllData[Data]["TransformedVariates"][AllData[Data]["TransformedFeatures"].index(Feature)] \
                                                for Data in Arguments.Data])))
    FormattedFeatures, FormattedVariates = [], [] 
    for Feature in set.intersection(*[set(AllData[Data]["FormattedFeatures"]) for Data in Arguments.Data]):
        FormattedFeatures.append(Feature + ":" + Arguments.Filename)
        FormattedVariates.append(list(chain(*[AllData[Data]["FormattedVariates"][AllData[Data]["FormattedFeatures"].index(Feature)] \
                                              for Data in Arguments.Data])))
        
    DataDict = {}
    DataDict["Labels"] = Labels
    DataDict["FormattedFeatures"] = FormattedFeatures
    DataDict["FormattedVariates"] = FormattedVariates
    DataDict["TransformedFeatures"] = TransformedFeatures
    DataDict["TransformedVariates"] = TransformedVariates
    DataDict["Report"] = make_report(Data, Labels, TransformedFeatures, FormattedFeatures, FormattedVariates, Arguments)
    
    cPickle.dump(DataDict, open(get_path("MOCA.data") + "/" + Arguments.Filename, "wb"), -1)

    if Arguments.Phenotype:
        Phenotypes, Variates = [], []
        for Data in Arguments.Data:
            if Data in Arguments.Phenotype:
                Phenotypes.append(Data + ":" + Arguments.Filename)
                Variates.append([1 if Label in AllData[Data]["Labels"] else 0 for Label in Labels])

        PhenotypeDict = {}
        PhenotypeDict["Labels"] = Labels
        PhenotypeDict["FormattedFeatures"] = Phenotypes
        PhenotypeDict["FormattedVariates"] = Variates
        PhenotypeDict["TransformedFeatures"] = Phenotypes
        PhenotypeDict["TransformedVariates"] = Variates
        PhenotypeDict["Report"] = make_report(Data, Labels, Phenotypes, Phenotypes, Variates, Arguments)
    
        cPickle.dump(PhenotypeDict, open(get_path("MOCA.data") + "/" + Arguments.Filename + ".Phenotypes", "wb"), -1)

    return 
