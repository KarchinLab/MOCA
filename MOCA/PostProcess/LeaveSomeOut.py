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

def make_report(Labels, Report):
    '''
    Having selected markers from the leave-some-out calculation, we are now gonna save a standard setworks
    output file. This is function just changes one of the leave-one-out-calculation reports to a standard 
    setworks report. 
    '''

    Report["Entries"].remove("Cases")
    Report["Entries"].remove("Controls")
    Report["Entries"].remove("CrossValidation")
    Report["Entries"].insert(2, "Label count")
    Report["Label count"] = len(Labels)
    
    return Report

def load_validation_data(Arguments):
    '''
    Load the python pickles you made from LeaveSomeOut cross-validation calculations
    '''

    CrossValidations = dict([(File, {}) for File in os.listdir(get_path("MOCA.results")) if Arguments.Filename + ".Validation" in File])
    for CrossValidation in CrossValidations.keys():
        Results = cPickle.load(open(get_path("MOCA.results") + "/" + CrossValidation, "rb"))
    
        CrossValidations[CrossValidation]["Cases"] = Results["Report"]["Cases"][1]
        CrossValidations[CrossValidation]["Controls"] = Results["Report"]["Controls"][1]
        CrossValidations[CrossValidation]["Barcodes"] = Results["Barcodes"]  
        CrossValidations[CrossValidation]["PValues"] = Results["PValues"] 
        CrossValidations[CrossValidation]["QValues"] = Results["QValues"] 
        CrossValidations[CrossValidation]["Performances"] = Results["Performances"]
        CrossValidations[CrossValidation]["Interactions"] = Results["Interactions"]
        CrossValidations[CrossValidation]["FeatureVectors"] = Results["FeatureVectors"]
        CrossValidations[CrossValidation]["UnionFeatures"] = Results["UnionFeatures"]
        CrossValidations[CrossValidation]["IntersectionFeatures"] = Results["IntersectionFeatures"] 
        CrossValidations[CrossValidation]["DifferenceFeatures"] = Results["DifferenceFeatures"]
        CrossValidations[CrossValidation]["Report"] = Results["Report"]

    return CrossValidations  

def leave_some_out(Arguments):
    ''' 
    Function for selecting biomarkers from cross-validation output. Strict in that it only returns 
    biomarkers that were selected during each cross validation. Might have the benefit, relative 
    to vote-based prediction, that these were so predictive that they will translate better to future 
    predictions. Also lends itself to simple clinical use because they require little or no computational 
    support for subsequent prediction...it is simply the marker
    '''

    CrossValidations = load_validation_data(Arguments)
    Labels, Features, Variates, Phenotypes, Markers = get_supervised_dataset(Arguments)
    Phenotype = Phenotypes[0]
    Response = Variates[Features.index(Phenotype)]
    Setworks = {}
    for CrossValidation in CrossValidations:
        for Barcode in CrossValidations[CrossValidation]["Barcodes"]:
            Setwork = (CrossValidations[CrossValidation]["UnionFeatures"][Barcode], \
                           CrossValidations[CrossValidation]["IntersectionFeatures"][Barcode], \
                           CrossValidations[CrossValidation]["DifferenceFeatures"][Barcode], \
                           CrossValidations[CrossValidation]["Interactions"][Barcode])

            if Setworks.has_key(Setwork): Setworks[Setwork].append(Barcode)
            else: Setworks[Setwork] = [Barcode]

    PValues = {}
    QValues = {}
    Performances = {}
    Interactions = {}
    FeatureVectors = {}
    UnionFeatures = {}
    IntersectionFeatures = {}
    DifferenceFeatures = {}
    SampleCounts = {}
    CaseCounts = {}
    EffectSizes = {}
    
    Barcodes = []
    for Setwork in Setworks: 
        if len(Setworks[Setwork]) == len(CrossValidations): #Sework had to be selected in each cross validation!!!
            Union, Intersection, Difference, Interaction = Setwork
            Predictor = assemble_setwork(Features, Variates, Union, Intersection, Difference, Arguments)
            TP,FP,FN,TN = contingency_table(Predictor, Response, NA=Arguments.NA)

            Barcode = Setworks[Setwork][0] 
            Barcodes.append(Barcode)
            PValues[Barcode] = fisher(TP,FP,FN,TN).two_tail
            QValues[Barcode] = "NA"
            Performances[Barcode] = Performance(Interaction, TP,FP,FN,TN)
            Interactions[Barcode] = Interaction
            EffectSizes[Barcode] = EffectSize(Interactions[Barcode], TP,FP,FN,TN)
            FeatureVectors[Barcode] = Predictor
            UnionFeatures[Barcode] = Union
            IntersectionFeatures[Barcode] = Intersection
            DifferenceFeatures[Barcode] = Difference
            SampleCounts[Barcode] = TP + FP + FN + TN
            CaseCounts[Barcode] = TP + FN

    Results = {}
    Results["PValues"] = PValues
    Results["QValues"] = QValues
    Results["Performances"] = Performances 
    Results["Interactions"] = Interactions
    Results["FeatureVectors"] = FeatureVectors
    Results["UnionFeatures"] = UnionFeatures
    Results["IntersectionFeatures"] = IntersectionFeatures
    Results["DifferenceFeatures"] = DifferenceFeatures
    Results["SampleCounts"] = SampleCounts
    Results["CaseCounts"] = CaseCounts
    Results["EffectSizes"] = EffectSizes
    #Doesn't matter which index we use we just need one report. The last accessed 'CrossValidation' will do. 
    Results["Report"] = make_report(Labels, CrossValidations[CrossValidation]["Report"])
    Results["Labels"] = Labels
    Results["Barcodes"] = Barcodes
    Results["Phenotype"] = Response

    if Arguments.Filename.lower() == "default":
        DataTypes = set(Arguments.Data).difference(set([Arguments.Phenotype]))
        Pickle = "_".join(["Phenotype=" + Phenotype[:Phenotype.index(":")],
                           "_".join(sorted(DataTypes)), str(Arguments.FeatureMin), Arguments.CorrectionMethod,
                           "".join(map(str, Arguments.BooleanSets)), "".join(map(str, Arguments.Optimization))])
    else:
        Pickle = Arguments.Filename + "_" + Phenotype[:Phenotype.index(":")]
        
    cPickle.dump(Results, open(get_path("MOCA.results") + "/" + Pickle, "wb"), -1)

    
    '''
    #Get rid of barcodes that for setworks that didn't pass a performance threshold (if provided)
    Barcodes = [Barcode for Barcode in Barcodes if minimum_performance(Performances[Barcode], Arguments)]

    #Initial sort will be done by decreasing balanced accuracy
    Barcodes = sorted(Barcodes, key=lambda Barcode: \
                      (Performances[Barcode].sensitivity + Performances[Barcode].specificity)/2, reverse=True)

    CSVfile = open(Arguments.Filename + ".csv", "wb")
    CSVwriter = csv.writer(CSVfile, dialect='excel')

    #Right the excel header
    CSVwriter.writerow(["Union","Intersection","Difference","Interaction", "Phenotype",
                        "Sensitivity","Specificity","PPV","NPV","Accuracy", "Sample Count"])
    
    for Barcode in Barcodes:
        p = Performances[Barcode]
        Sens, Spec, PPV, NPV, Accuracy = p.sensitivity, p.specificity, p.PPV, p.NPV, p.accuracy 
        CSVwriter.writerow([", ".join(UnionFeatures[Barcode]), ", ".join(IntersectionFeatures[Barcode]),
                           ", ".join(DifferenceFeatures[Barcode]), Interactions[Barcode],
                            Phenotype[:Phenotype.index(":")], "%0.2f" %Sens, "%0.2f" %Spec,
                            "%0.2f" %PPV, "%0.2f" %NPV, "%0.2f" %Accuracy, SampleCount[Barcode]])

    CSVfile.close()
    '''
    
    return 

 
