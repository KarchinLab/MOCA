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

def load_pairwise_results(Results):
    '''
    Simply load the results from your MOCA pairwise calculations
    '''

    PValues = Results["PValues"]
    FDRs = Results["FDRs"]
    Interactions = Results["Interactions"]
    Performances = Results["Performances"]
    SampleCounts = Results["SampleCounts"]
    CaseCounts = Results["CaseCounts"]
    EffectSizes = Results["EffectSizes"]

    return PValues, FDRs, Interactions, Performances, SampleCounts, CaseCounts, EffectSizes

def unsupervised(Results, Arguments):
    '''
    Write a csv file of your unsupervised pairwise interactions and a few related states. Suggested for viewing in 
    MS excel (cross OS), number (os x), or the freely available OpenOffice Calc program (cross OS).
    '''

    PValues, FDRs, Interactions, Performances, SampleCounts, CaseCounts, EffectSizes = load_pairwise_results(Results)

    CSVfile = open(Arguments.Filename + ".csv", "wb")
    CSVwriter = csv.writer(CSVfile, dialect='excel')
    
    #Right the excel header
    CSVwriter.writerow(["Feature 1", "Interaction", "Feature 2", "P-value", "FDR", "Sample Count"])

    #sort by p-value, then alphabetically on the features
    for Feature, PValue in sorted(PValues.items(), key=lambda x: (x[1], x[0][0], x[0][1]))[:Arguments.TopInteractions]: 
        CSVwriter.writerow([Feature[0],  Interactions[Feature],  Feature[1],
                            "%0.2e" %PValue, "%0.2e" %FDRs[PValue], SampleCounts[Feature]])

    CSVfile.close()
    
    return

def supervised(Results, Arguments):
    '''
    Write a csv file of your unsupervised pairwise interactions and a few related states. Suggested for viewing in 
    MS excel (cross OS), number (os x), or the freely available OpenOffice Calc program (cross OS).
    '''

    PValues, FDRs, Interactions, Performances, SampleCounts, CaseCounts, EffectSizes = load_pairwise_results(Results)

    Report = Results["Report"] #Need to get the phenotype name because this is 'supervised'
    Phenotype = Report["Supervised"][1].replace("(","").replace(")", "")

    CSVfile = open(Arguments.Filename + ".csv", "wb")
    CSVwriter = csv.writer(CSVfile, dialect='excel')
    
    #Right the excel header
    CSVwriter.writerow(["Feature", "Interaction", "Phenotype", "P-value", "FDR", "Odds Ratio", "Effect Size",
                        "Sensitivity", "95% CI", "Specificity", "95% CI", "PPV", "95% CI", "NPV", "95% CI",
                        "Accuracy", "MCC", "Sample Count", "Case Count"])
    
    #just sort by p-value
    for Feature in rank(Performances, Arguments.RankMethod, Arguments.TopInteractions):
        p = Performances[Feature]
        Sens, Spec, PPV, NPV, Accuracy, MCC = p.sensitivity, p.specificity, p.PPV, p.NPV, p.accuracy, p.MCC
        if minimum_performance(Performances[Feature], Arguments):
            CSVwriter.writerow([Feature, Interactions[Feature], Phenotype, "%0.2e" %PValues[Feature], "%0.2e" %FDRs[PValues[Feature]],
                                "%0.2f" %EffectSizes[Feature].odds_ratio, "%0.2f" %EffectSizes[Feature].difference_of_proportions,
                                "%0.2f" %Sens, "-".join(["%0.2f" %CI for CI in confidence_interval(Sens, SampleCounts[Feature])]),
                                "%0.2f" %Spec, "-".join(["%0.2f" %CI for CI in confidence_interval(Spec, SampleCounts[Feature])]),
                                "%0.2f" %PPV, "-".join(["%0.2f" %CI for CI in confidence_interval(PPV, SampleCounts[Feature])]),
                                "%0.2f" %NPV, "-".join(["%0.2f" %CI for CI in confidence_interval(NPV, SampleCounts[Feature])]),
                                "%0.2f" %Accuracy, "%0.2f" %MCC, SampleCounts[Feature], CaseCounts[Feature]])

    CSVfile.close()
    
    return 

def pairwise_continuous(Results, Arguments):
    '''
    '''

    PValues = Results["PValues"]
    FDRs = Results["FDRs"]
    Correlations = Results["Correlations"]

    CSVfile = open(Arguments.Filename + ".csv", "wb")
    CSVwriter = csv.writer(CSVfile, dialect='excel')

    #Right the excel header
    CSVwriter.writerow(["Feature 1", "Feature 2", "Pearson P-value", "FDR", "Pearson Correlation"])

    #sort by p-value, then alphabetically on the features
    for Feature, PValue in sorted(PValues.items(), key=lambda x: (x[1], x[0][0], x[0][1]))[:Arguments.TopInteractions]: 
        CSVwriter.writerow([Feature[0],  Feature[1],
                            "%0.2e" %PValue, "%0.2e" %FDRs[PValue], "%0.2f" %Correlations[Feature]])

    CSVfile.close()
    
    return

def pairwise(Arguments):
    '''
    Print results from you MOCA pairwise runs. 
    '''

    Results = cPickle.load(open(get_path("MOCA.results") + "/" + Arguments.Filename, "rb"))

    if not Results["Report"]["Pairwise"]:
        print "You have 'Pairwise = True' in your arguments file, yet the file you pointed to (MOCA.results/",
        print Arguments.Filename + ") did not result from a pairwise MOCA calculation."
        print "Exiting..."
        exit()

    if Results["Report"]["Continuous-valued correlation"]:
        pairwise_continuous(Results, Arguments)
    elif type(Results["Report"]["Supervised"]) == list:
        supervised(Results, Arguments)
    else:
        unsupervised(Results, Arguments)
    
    return 
