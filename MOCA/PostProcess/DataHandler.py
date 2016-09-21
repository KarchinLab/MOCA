'''
Get python modules
'''
from itertools import chain
import csv

'''
Get third-party modules
'''

'''
Get MOCA modules
'''
from MOCA.DataHandler import load_data, get_path
from Setworks import load_setworks

def load_results_file(ResultsFile):
    '''
    '''

    Results = [Line for Line in csv.reader(open(ResultsFile,"rb"))]
    
    return Results[0], Results[1:]

def make_priors(Arguments):
    '''
    This function can make a Priors file from any setwork output (probably from either a Setworks run, or a 
    ValidateBiomarkers run). This Priors file can then be used to bias future searches. Can be useful for 
    focusing setwork selection from "big data". For instance, do a Setworks run to see which individual features
    are being selected with high frequency. Next, make the Priors file and rerun Setworks forcing MOCA to bias 
    its search around the individual features predicted to be important in your previous run. This can be useful, 
    because some data is too big to search thru with the normal MOCA optimization process, without a little help--
    this would be that help. 
    '''

    Header, Markers = load_results_file(get_path("Priors"))
    
    UnionFeatures, IntersectionFeatures, DifferenceFeatures = [], [], []
    for Marker in Markers:
        if Marker[0]: UnionFeatures.extend(map(lambda x: x.strip(","), Marker[0].split()))
        if Marker[1]: IntersectionFeatures.extend(map(lambda x: x.strip(","), Marker[1].split()))
        if Marker[2]: DifferenceFeatures.extend(map(lambda x: x.strip(","), Marker[2].split()))

    Priors = open(Arguments.Filename, "w")
    Priors.write("%s %s \n" %("Union =", " ".join(UnionFeatures)))
    Priors.write("%s %s \n" %("Intersection =", " ".join(IntersectionFeatures)))
    Priors.write("%s %s \n" %("Difference =", " ".join(DifferenceFeatures)))
    Priors.close()

    return
    
def write_matrix(Arguments):
    '''
    Take any data that has been preprocessed (i.e., exists in your projects MOCA.data folder) and write a human 
    readable version to the current working directory. Useful if you need to load the matrices into a different utility, 
    tho why would anyone need to do that :) Also useful to load previously selected setworks back into MOCA for further
    optimization. You can also use the LabelList and FeatureList arguments to further customize what goes into these 
    matrices (that part is handled in the main DataHandler:load_data). 
    '''

    Data = load_data(Arguments)
    
    Labels = Data.Labels

    if set(Arguments.UntransformedData).difference(set(Arguments.Data)):
        print "It appears you have some datatypes passed to 'Untransformed = ' that you didn't load via 'Data ='."
        print "Only datatypes loaded via 'Data = ' can be written to matrix...transformed or not!"

    if set(Arguments.Data).difference(set(Arguments.UntransformedData)):
        Features = [Data.Transformed.Features[DataType] for DataType in Arguments.Data \
                    if DataType not in Arguments.UntransformedData]
        Variates = [Data.Transformed.Variates[DataType] for DataType in Arguments.Data \
                    if DataType not in Arguments.UntransformedData]
    
        Features = list(chain(*Features))
        Variates = list(chain(*Variates))
    
        CSVfile = open(Arguments.Filename + ".csv", "wb")
        CSVwriter = csv.writer(CSVfile, dialect='excel')
        CSVwriter.writerow(Labels)
        for Feature in Features:
             CSVwriter.writerow([Feature] + Variates[Features.index(Feature)])
        CSVfile.close()
        
    if Arguments.UntransformedData:
        Features = [Data.Features[DataType] for DataType in Arguments.Data \
                    if DataType in Arguments.UntransformedData]
        Variates = [Data.Variates[DataType] for DataType in Arguments.Data \
                    if DataType in Arguments.UntransformedData]
    
        Features = list(chain(*Features))
        Variates = list(chain(*Variates))

        CSVfile = open(Arguments.Filename + "_untransformed.csv", "wb")
        CSVwriter = csv.writer(CSVfile, dialect='excel')
        CSVwriter.writerow(Labels)
        for Feature in Features:
            CSVwriter.writerow([Feature] + Variates[Features.index(Feature)])
        CSVfile.close()
    
    return

def write_setworks(Arguments):
    '''
    Useful to load previously selected setworks back into MOCA for further
    optimization. 
    '''

    Barcodes, PValues, QValues, Performances, Interactions, \
        UnionFeatures, IntersectionFeatures, DifferenceFeatures, \
        FeatureVectors, SampleCounts, CaseCounts, EffectSizes, Report, Labels = load_setworks(Arguments)

    CSVfile = open(Arguments.Filename + ".csv", "wb")
    CSVwriter = csv.writer(CSVfile, dialect='excel')

    CSVwriter.writerow(Labels)
    for Barcode in Barcodes:
        Row = [" ".join([str(tuple([" | ".join(UnionFeatures[Barcode])])), " & ", \
                                      str(tuple([" & ".join(IntersectionFeatures[Barcode])])), " - ", \
                                      str(tuple([" | ".join(DifferenceFeatures[Barcode])]))])] + FeatureVectors[Barcode]
        CSVwriter.writerow(Row)

    CSVfile.close()
    
    return


