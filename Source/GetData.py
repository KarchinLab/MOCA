'''
Get python modules
'''
from os import environ
from collections import Counter
from ast import literal_eval

'''
Get third-party modules
'''

'''
Get MOCA modules
'''
from Statistics import binary_matrix, z_matrix

def get_labels(File, Row=0, Start=0):
    '''
    Gets desired "Row" (default=0). Initially written to get
    header (i.e., Row=0), but can get any row. "Start" (default=
    0) is the column to start at. Useful if the 0,0 coordinate
    of a matrix is not a column label, for instance. 
    '''
            
    return open(File).readlines()[Row].split()[Start:]

def get_features(File, Column=0, Start=1):
    '''
    Gets desired 'Column' (default=0). Initially written to get
    row names (e.g., gene or drug names), but can get any column.
    'Start' (default=1) is the row to start at (defaults to '1' 
    because it is assumned there is a header (e.g., sample or 
    patient IDs).
    '''

    return [Line.split()[Column] for Line in file(File)][Start:]

def get_variates(File, Row=1, Column=1, Type=str, NA="NA"):
    '''
    Returns the coordinates of an nxn matrix. Assumes there are header
    and row names (i.e., zeroth row and column are NOT part of the matrix
    and therefore they both default to '1'). 'Type' defaults to 'str',
    but can also be 'float' or 'int'. 
    '''

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

def get_file(DataType, PathFile):
    '''
    For a particular data type (e.g., drug or expression) find the 
    file containing said data using moca paths file. Absolute paths
    can start with ${HOME} to indicate "home/user"
    '''

    Path = [Line.split()[2] for Line in open(PathFile).readlines() \
                if Line.strip() and Line.split()[0] == DataType][0]
    if Path[:7] == "${HOME}":
        Path = Path.replace("${HOME}", environ["HOME"])
    
    return Path

def is_binary(Variates, NA="NA"):
    '''
    Reduces the variate matrix to a non redundant vector. If the input
    matrix is binary, then the reduced vector must either be only 2
    (i.e., [0, 1]) or 1 (i.e., [0] or [1]) elements in length. 
    '''

    ReducedMatrix = list(set([Variate for Vector in Variates for Variate in Vector]))
    ReducedMatrix = filter(lambda Variate: Variate != NA, ReducedMatrix)
    for Variate in ReducedMatrix:
        if not Variate.isdigit():
            IsBinary = False
            break
        else:
            IsBinary = True

    return IsBinary

def get_data_type(Feature):
    '''
    Given a feature with that might include the MOCA-naming convention 
    (i.e., ":Negative/Positive") returns the data type (e.g., Expression or Drug
    '''

    return Feature[Feature.index(":")+1:].replace("Positive", " Positive").replace("Negative", " Negative").split()[0]

def get_data(Arguments):
    '''
    Tell the function which types of data to fetch and what z-score
    cutoff use to convert continuous data to binary. Returns labels
    (e.g., cell lines or tumors), variates (e.g., GI50s or microarray 
    values), and features (e.g., genes or drugs) in dictionary form
    '''
    
    Labels = [] #e.g., tumors or cell lines
    DataDict = dict()

    #Get labels and only consider those present in each data type
    for Datum in Arguments.Data:
        Labels.append(set(get_labels(get_file(Datum, Arguments.PathFile))))
    Labels = list(set.intersection(*Labels))

    if Arguments.TissueSpecific:
        Labels = [Label for Label in Labels if Label[:Label.index(":")] == Arguments.TissueSpecific]
    DataDict["Labels"] = Labels
   
    for Datum in Arguments.Data:
        if is_binary(get_variates(get_file(Datum, Arguments.PathFile))):
            #Get the variates and corresponding features for binary data (e.g., mutation)
            FeatureDict = []
            VariateDict = []
            Features = get_features(get_file(Datum, Arguments.PathFile))
            Variates = get_ordered_matrix(Labels, get_labels(get_file(Datum, Arguments.PathFile)), \
                                              get_variates(get_file(Datum, Arguments.PathFile), Type=int))

            for Feature in Features:
                if sum([1 if Variate == 1 else 0 for Variate in Variates[Features.index(Feature)]]) >= Arguments.FeatureMin:
                    if ":" not in Feature:
                        FeatureDict.append(Feature + ":" + Datum)
                        VariateDict.append(Variates[Features.index(Feature)])
                    else:
                        FeatureDict.append(Feature)
                        VariateDict.append(Variates[Features.index(Feature)])
            
            DataDict[Datum + "Features"] = FeatureDict 
            DataDict[Datum + "Variates"] = VariateDict

        elif Arguments.Continuous:
            #You want to leave the data as is, without any binarization
            FeatureDict = []
            VariateDict = []
            Features = get_features(get_file(Datum, Arguments.PathFile))
            Variates = get_ordered_matrix(Labels, get_labels(get_file(Datum, Arguments.PathFile)), \
                                              get_variates(get_file(Datum, Arguments.PathFile), Type=float))

            for Feature in Features:
                if ":" not in Feature:
                    FeatureDict.append(Feature + ":" + Datum)
                    VariateDict.append(Variates[Features.index(Feature)])
                else:
                    FeatureDict.append(Feature)
                    VariateDict.append(Variates[Features.index(Feature)])
            
            DataDict[Datum + "Features"] = FeatureDict 
            DataDict[Datum + "Variates"] = VariateDict
       
        else:
            '''
            Get features and corresponding variates. Because this is continuous data
            we need to convert to a z-score matrix and then discretize using a z-score
            cutoff. We also need to return a "positive" and "negative" matrix to 
            account for both extremes of the datum
            '''
            
            #Get the variates and corresponding features for continuous data (e.g., drug)
            FeatureDict = []
            VariateDict = []
            
            Features = get_features(get_file(Datum, Arguments.PathFile))
            Variates = get_variates(get_file(Datum, Arguments.PathFile))
        
            #Z-tranform
            VariatePositive = get_ordered_matrix(Labels, get_labels(get_file(Datum, Arguments.PathFile)), \
                                                     binary_matrix(z_matrix(Variates), \
                                                                       CutOff=Arguments.ZMin, GreaterThan=True))
            VariateNegative = get_ordered_matrix(Labels, get_labels(get_file(Datum, Arguments.PathFile)), \
                                                     binary_matrix(z_matrix(Variates), \
                                                                       CutOff=-1*Arguments.ZMin, GreaterThan=False))
            
            for Feature in Features:
                if sum([1 if Variate == 1 else 0 for Variate in VariatePositive[Features.index(Feature)]]) >= Arguments.FeatureMin:
                    FeatureDict.append(Feature + ":" + Datum + "Positive")
                    VariateDict.append(VariatePositive[Features.index(Feature)])
                if sum([1 if Variate == 1 else 0 for Variate in VariateNegative[Features.index(Feature)]]) >= Arguments.FeatureMin:
                    FeatureDict.append(Feature + ":" + Datum + "Negative")
                    VariateDict.append(VariateNegative[Features.index(Feature)])
            DataDict[Datum + "Features"] = FeatureDict 
            DataDict[Datum + "Variates"] = VariateDict
    
    return DataDict

def get_master_matrix(Arguments):
    '''
    Gets data and makes dictionary matrix that holds
    ALL data!
    '''

    Data = get_data(Arguments)
    DataDict = dict()
    FeatureDict = []
    VariateDict = []

    for Datum in set(Arguments.Data):
        Features = Data[Datum + "Features"]
        Variates = Data[Datum + "Variates"]
        for Feature in Features:
            FeatureDict.append(Feature)
            VariateDict.append(Variates[Features.index(Feature)])

    DataDict["Labels"] = Data["Labels"]
    DataDict["Features"] = FeatureDict
    DataDict["Variates"] = VariateDict

    return DataDict

def make_binary_matrix(Arguments):
    '''
    Make a binary matrix out of your raw or Z-score data. This can be very useful for two reasons:
    1) If you matrix is huge, you don't want to binarize everytime, as it takes a long time.
    2) If you want to use differential Z-score cutoffs and feature mins for each data type you are
    going to have to do the conversion this way.
    '''
    
    Labels = [] #e.g., tumors or cell lines
 
    #Get labels and only consider those present in each data type
    for Datum in Arguments.Data:
        Labels.append(set(get_labels(get_file(Datum, Arguments.PathFile))))
    Labels = list(set.intersection(*Labels))

    if Arguments.TissueSpecific:
        Labels = [Label for Label in Labels if Label[:Label.index(":")] == Arguments.TissueSpecific]
    
    #MultiprocessMode support if called
    Node = int(Arguments.MultiProcessMode[0]) 
    TotalNodes = int(Arguments.MultiProcessMode[1])
    
    for Datum in Arguments.Data:
        if is_binary(get_variates(get_file(Datum, Arguments.PathFile))):

            #Get the variates and corresponding features for binary data (e.g., mutation)
            Features = get_features(get_file(Datum, Arguments.PathFile))[Node::TotalNodes]
            Variates = get_ordered_matrix(Labels, get_labels(get_file(Datum, Arguments.PathFile)), \
                                              get_variates(get_file(Datum, Arguments.PathFile), Type=int))[Node::TotalNodes]

            BinaryMatrix = open(Datum, "w")
            BinaryMatrix.write("%s \n" %(" ".join(Labels)))
            for Feature in Features:
                if ":" not in Feature:
                    BinaryMatrix.write("%s %s \n" %(Feature + ":" + Datum, \
                                                        " ".join(map(str, Variates[Features.index(Feature)]))))
                else:
                    BinaryMatrix.write("%s %s \n" %(Feature, \
                                                        " ".join(map(str, Variates[Features.index(Feature)]))))
            BinaryMatrix.close()
            
        else:
            
            '''
            Get features and corresponding variates. Because this is continuous data
            we need to convert to a z-score matrix and then discretize using a z-score
            cutoff. We also need to return a "positive" and "negative" matrix to 
            account for both extremes of the datum
            '''
        
            #Get the variates and corresponding features for continuous data (e.g., drug)
            Features = get_features(get_file(Datum, Arguments.PathFile))[Node::TotalNodes]
            Variates = get_variates(get_file(Datum, Arguments.PathFile))[Node::TotalNodes]
            
            if Arguments.ZData: #already Z-transformed; don't call "z_matrix"
                VariatePositive = get_ordered_matrix(\
                    Labels, get_labels(get_file(Datum, Arguments.PathFile)), \
                        binary_matrix(Variates, CutOff=Arguments.ZMin, GreaterThan=True))
                VariateNegative = get_ordered_matrix(\
                    Labels, get_labels(get_file(Datum, Arguments.PathFile)), \
                        binary_matrix(Variates, CutOff=-1*Arguments.ZMin, GreaterThan=False))
        
            else: #Not Z-transformed; call "z_matrix"
                VariatePositive = get_ordered_matrix(\
                    Labels, get_labels(get_file(Datum, Arguments.PathFile)), \
                        binary_matrix(z_matrix(Variates), CutOff=Arguments.ZMin, GreaterThan=True))
                VariateNegative = get_ordered_matrix(\
                    Labels, get_labels(get_file(Datum, Arguments.PathFile)), \
                        binary_matrix(z_matrix(Variates), CutOff=-1*Arguments.ZMin, GreaterThan=False))
                
            BinaryMatrix = open(Datum + "(ZMin=" + str(Arguments.ZMin) + ")", "w")
            BinaryMatrix.write("%s \n" %(" ".join(Labels)))
            for Feature in Features:
                BinaryMatrix.write("%s %s \n" %(Feature + ":" + Datum + "Positive" + str(Arguments.ZMin), \
                                                    " ".join(map(str, VariatePositive[Features.index(Feature)]))))
                BinaryMatrix.write("%s %s \n" %(Feature + ":" + Datum + "Negative" + str(Arguments.ZMin), \
                                                    " ".join(map(str, VariateNegative[Features.index(Feature)]))))
            BinaryMatrix.close()
    
    return

def get_priors(Arguments):
    '''
    Read in a priors file. Priors are useful to focus the search, when making setworks, around features
    you think might be important. Make your own Priors file, or use MakePriors to make a priors file 
    from a previous Setworks run. 
    '''

    UnionPriors, IntersectionPriors, DifferencePriors = (), (), ()

    Priors = [Line.split() for Line in open(get_file("Priors", Arguments.PathFile)).readlines()]

    for Prior in Priors:
        if Prior[0].lower() == "union":
            UnionPriors = tuple(Prior[2:])
        if Prior[0].lower() == "intersection":
            IntersectionPriors = tuple(Prior[2:])
        if Prior[0].lower() == "difference":
            DifferencePriors = tuple(Prior[2:])
    
    return UnionPriors, IntersectionPriors, DifferencePriors
    
