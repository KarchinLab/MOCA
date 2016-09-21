'''
Get python modules
'''
import cPickle
import os
import csv

'''
Get third-party modules
'''
from fisher import pvalue as fisher

'''
Get MOCA modules
'''
from MOCA.DataHandler import get_path, get_supervised_dataset
from MOCA.Statistics import Performance, contingency_table, EffectSize
from MOCA.Setworks import assemble_setwork
from DataHandler import load_results_file

def validate_markers(Arguments):
    '''
    '''

    Labels, Features, Variates, Phenotypes, Markers = get_supervised_dataset(Arguments)
    Header, Results = load_results_file(get_path("ValidateMarkers"))

    CSVfile = open(Arguments.Filename + ".csv", "wb")
    CSVwriter = csv.writer(CSVfile, dialect='excel')
    CSVwriter.writerow(["Union","Intersection","Difference","Interaction", "Phenotype", "P-value", "Odds Ratio", "Effect Size",
                        "Sensitivity","Specificity","PPV","NPV","Accuracy", "MCC", "Sample Count", "Case Count"])
    for Phenotype in Phenotypes:
        Response = Variates[Features.index(Phenotype)]
        for Marker in Results:
            try:
                Predictor = assemble_setwork(Features, Variates,
                                             filter(None, Marker[Header.index("Union")].split(", ")),
                                             filter(None, Marker[Header.index("Intersection")].split(", ")),
                                             filter(None, Marker[Header.index("Difference")].split(", ")), Arguments)
                
                TP,FP,FN,TN = contingency_table(Predictor, Response, NA=Arguments.NA)
                performance = Performance(Marker[Header.index("Interaction")], TP,FP,FN,TN)
                effect_size = EffectSize(Marker[Header.index("Interaction")], TP,FP,FN,TN)
                CSVwriter.writerow([Marker[Header.index("Union")], Marker[Header.index("Intersection")], Marker[Header.index("Difference")],
                                    Marker[Header.index("Interaction")], Phenotype[:Phenotype.index(":")], "%0.2e" %fisher(TP,FP,FN,TN).two_tail,
                                    "%0.2f" %effect_size.odds_ratio, "%0.2f" %effect_size.difference_of_proportions, "%0.2f" %performance.sensitivity,
                                    "%0.2f" %performance.specificity, "%0.2f" %performance.PPV, "%0.2f" %performance.NPV,
                                    "%0.2f" %performance.accuracy, "%0.2f" %performance.MCC, TP+FP+FN+TN, TP+FN])
            except ValueError:
                CSVwriter.writerow([Marker[Header.index("Union")], Marker[Header.index("Intersection")], Marker[Header.index("Difference")], "NA"])
                
    CSVfile.close()
    
    return
