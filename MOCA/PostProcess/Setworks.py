'''
Get python modules
'''
import cPickle
import csv

'''
Get third-party modules
'''

'''
Get MOCA modules
'''
from MOCA.DataHandler import get_path
from MOCA.Statistics import minimum_performance, rank, confidence_interval
                        
def load_setworks(Arguments):
    '''
    Load the setworks pickled dictionary of dictionaries, and return lists that have been truncated
    to specified number of results (TopInteractions = 100, by default). The list can additionally be
    truncated using performance metric, thru the MinimumPerformance argument. 
    '''
   
    Results = cPickle.load(open(get_path("MOCA.results") + "/" + Arguments.Filename, "rb"))

    PValues = Results["PValues"] 
    QValues = Results["QValues"] 
    Performances = Results["Performances"]
    Interactions = Results["Interactions"] 
    FeatureVectors = Results["FeatureVectors"]
    UnionFeatures = Results["UnionFeatures"]
    IntersectionFeatures = Results["IntersectionFeatures"] 
    DifferenceFeatures = Results["DifferenceFeatures"]
    SampleCounts = Results["SampleCounts"]
    CaseCounts = Results["CaseCounts"]
    EffectSizes = Results["EffectSizes"]
    Barcodes = Results["Barcodes"]
    Report = Results["Report"]
    Labels = Results["Labels"]

    #Get rid of barcodes that for setworks that didn't pass a performance threshold (if provided)
    Barcodes = [Barcode for Barcode in Barcodes if minimum_performance(Performances[Barcode], Arguments)]

    #Initial sort will be done by decreasing balanced accuracy
    Barcodes = rank(Performances, Arguments.RankMethod, Arguments.TopInteractions)
    
    return Barcodes, PValues, QValues, Performances, Interactions, \
        UnionFeatures, IntersectionFeatures, DifferenceFeatures, \
        FeatureVectors, SampleCounts, CaseCounts, EffectSizes, Report, Labels


def setworks(Arguments):
    '''
    Write a csv file of your setworks, plus a ton of potentially relevant stats. Suggested for viewing in MS excel (cross OS), 
    number (os x), or the freely available OpenOffice Calc program (cross OS).
    '''

    Barcodes, PValues, QValues, Performances, Interactions, \
        UnionFeatures, IntersectionFeatures, DifferenceFeatures, \
        FeatureVectors, SampleCounts, CaseCounts, EffectSizes, Report, Labels = load_setworks(Arguments)

    CSVfile = open(Arguments.Filename + ".csv", "wb")
    CSVwriter = csv.writer(CSVfile, dialect='excel')

    #Right the excel header
    CSVwriter.writerow(["Union","Intersection","Difference","Interaction","Phenotype","P-value",
                        "FDR",  "Odds Ratio", "Effect Size", "Sensitivity", "95% CI", "Specificity",
                        "95% CI", "PPV", "95% CI", "NPV", "95% CI", "Accuracy", "MCC", "Sample Count", "Case Count"])
    
    for Barcode in Barcodes:
        p = Performances[Barcode]
        Sens, Spec, PPV, NPV, Accuracy, MCC = p.sensitivity, p.specificity, p.PPV, p.NPV, p.accuracy, p.MCC
        CSVwriter.writerow([", ".join(UnionFeatures[Barcode]), ", ".join(IntersectionFeatures[Barcode]),
                            ", ".join(DifferenceFeatures[Barcode]), Interactions[Barcode], Report["Phenotype"],
                            "%0.2e" %PValues[Barcode], ["NA" if QValues[Barcode] == "NA" else "%0.2e" %QValues[Barcode]][0],
                            "%0.2f" %EffectSizes[Barcode].odds_ratio, "%0.2f" %EffectSizes[Barcode].difference_of_proportions,
                            "%0.2f" %Sens, "-".join(["%0.2f" %CI for CI in confidence_interval(Sens, SampleCounts[Barcode])]),
                            "%0.2f" %Spec, "-".join(["%0.2f" %CI for CI in confidence_interval(Spec, SampleCounts[Barcode])]),
                            "%0.2f" %PPV, "-".join(["%0.2f" %CI for CI in confidence_interval(PPV, SampleCounts[Barcode])]),
                            "%0.2f" %NPV, "-".join(["%0.2f" %CI for CI in confidence_interval(NPV, SampleCounts[Barcode])]),
                            "%0.2f" %Accuracy, "%0.2f" %MCC, SampleCounts[Barcode], CaseCounts[Barcode]])

    CSVfile.close()

    return

