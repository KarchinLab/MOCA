'''
Get python modules
'''
from math import ceil
from random import shuffle
import cPickle

'''
Get third-party modules
'''

'''
Get MOCA modules
'''
from DataHandler import get_supervised_dataset, get_ordered_matrix, get_path
from Setworks import get_setworks, assemble_setwork
from Statistics import p_adjust, MOCASet, minimum_performance

def make_report(Cases, Controls, Phenotype, CrossValidation, Arguments):
    '''
    This report is retrieved for all MOCA results files when 'Reports = Results' is set in the 
    Arguments file. Tells you from which input data the results were generated, which file the 
    results were written to, and what parameters were used during the calculation.
    '''

    Report = {}
    #Use the 'Entries' key to set the order of output when 'Reports = Results' is run.
    Report["Entries"] = ["Input data", "Phenotype", "CrossValidation", "Cases", "Controls", 
                         "Correction method", "FDR", "Optimization parameters", "Boolean sets",
                         "Eject passengers", "Bandwidth", "Seed", "Phenotype permutation test?"]
    Report["Input data"] = [Data for Data in Arguments.Data]
    Report["Phenotype"] = Phenotype[:Phenotype.index(":")]
    Report["CrossValidation"] = CrossValidation
    Report["Cases"] = ["(" + str(len(Cases[CrossValidation])) + ")", Cases[CrossValidation]]
    Report["Controls"] = ["(" + str(len(Controls[CrossValidation])) + ")", Controls[CrossValidation]]
    Report["Correction method"] = Arguments.CorrectionMethod
    Report["FDR"] = Arguments.FDR
    Report["Optimization parameters"] = Arguments.Optimization
    Report["Boolean sets"] = Arguments.BooleanSets
    if not Arguments.EjectPassengers:
        Report["Eject passengers"] = False
    else:
        Report["Eject passengers"] = Arguments.EjectPassengers
    Report["Bandwidth"] = Arguments.Bandwidth
    Report["Seed"] = Arguments.Seed
    Report["Phenotype permutation test?"] = Arguments.PermutePhenotype
    
    return Report

def cross_validation(Arguments, \
                         Labels, Features, Variates, \
                         Markers, Phenotype, CrossValidations, Cases, Controls):
    '''
    Labels, Features, and Variates have the usual meanings. Markers are any data you 
    give MOCA that you didn't tell it is the Phenotype. Phenotype is the thing your 
    trying to predict (i.e., select markers for)--only one Phenotype per LeaveSomeOut run please.
    Interger numnber of CrossValidations to do (e.g., leave ONE out or TEN-fold cross validation)
    Cases and Controls are dictionaries specifying which labels go with which CrossValidation. 
    For each cross-validation this pickles a feature matrix of the following format
    
    Label1 Label2 Label3.....
    Phenotype 0 0 1.....
    (Setwork1 InteractionType) 1 0 1.....
    (Setwork2 InteractionType) 1 1 0.....
    (Setwork3 InteractionType) 0 0 0.....
    .
    .
    .
    .
    .
    '''
    
    #Read in the setwork optimization paraments
    Trials = int(Arguments.Optimization[0])
    RepopulateFrequency = int(Arguments.Optimization[1])
    PercentToRepopulate = float(Arguments.Optimization[2])

    #Read in the Boolean set parameters
    UnionFeatures = int(Arguments.BooleanSets[0])
    IntersectionFeatures = int(Arguments.BooleanSets[1])
    DifferenceFeatures = int(Arguments.BooleanSets[2])

    #Multiprocess mode support if called. Cross-validation is so compute intensive that you 
    #can only do this for one Phenotype at a time, and the different cross-validations are 
    #distributed to different processors if MultiProcessMode is called
    Node = int(Arguments.MultiProcessMode[0]) 
    TotalNodes = int(Arguments.MultiProcessMode[1])

    for CrossValidation in range(CrossValidations)[Node::TotalNodes]:
        #Split the data is directed
        TrainLabels = list(set(Labels) - set(Cases[CrossValidation] + Controls[CrossValidation]))
        TrainVariates = get_ordered_matrix(TrainLabels, Labels, Variates)
        
        #Get setworks from the training data
        PValues, Performances, Interactions, FeatureVectors, Setworks, SampleCounts, CaseCounts, EffectSizes = \
            get_setworks(Arguments, \
                             Features, TrainVariates, Phenotype, \
                             Markers, Markers, Markers, \
                             Trials, RepopulateFrequency, PercentToRepopulate, \
                             UnionFeatures, IntersectionFeatures, DifferenceFeatures)
        
        #Make the cross-validation matrix. First row is the case-control-label header
        CrossValidationFeatureMatrix = [Cases[CrossValidation] + Controls[CrossValidation]]
        PhenotypeVector = [1 for Case in Cases[CrossValidation]] + [0 for Control in Controls[CrossValidation]]
        CrossValidationFeatureMatrix.append(PhenotypeVector) #second row is the phenotype vector

        #We only need the intersection of unique setworks passing the FDR threshold
        QValues = p_adjust(PValues, Arguments.CorrectionMethod)
        QValues = dict([(Barcode, QValues[PValue]) for Barcode, PValue in PValues.items() \
                            if QValues[PValue] < Arguments.FDR])
        Barcodes = list(set.intersection(set(Setworks.keys()), set(QValues.keys())))

        #finally, if we desire we can filter by fdr at this stage. We could do it later, but we'll get a bigger Pickle now. 
        Barcodes = [Barcode for Barcode in Barcodes if minimum_performance(Performances[Barcode], Arguments)]

        Results = {}
        Results["PValues"] = dict([(Barcode, PValues[Barcode]) for Barcode in Barcodes])
        Results["QValues"] = dict([(Barcode, QValues[Barcode]) for Barcode in Barcodes])
        Results["Performances"] = dict([(Barcode, Performances[Barcode]) for Barcode in Barcodes]) 
        Results["Interactions"] = dict([(Barcode, Interactions[Barcode]) for Barcode in Barcodes])
        Results["FeatureVectors"] = dict([(Barcode, FeatureVectors[Barcode]) for Barcode in Barcodes]) 
        Results["UnionFeatures"] = dict([(Barcode, Setworks[Barcode][0]) for Barcode in Barcodes])
        Results["IntersectionFeatures"] = dict([(Barcode, Setworks[Barcode][1]) for Barcode in Barcodes])
        Results["DifferenceFeatures"] = dict([(Barcode, Setworks[Barcode][2]) for Barcode in Barcodes])
        Results["Report"] = make_report(Cases, Controls, Phenotype, CrossValidation, Arguments)
        Results["Barcodes"] = Barcodes 

        if Arguments.Filename.lower() == "default":
            DataTypes = set(Arguments.Data).difference(set([Arguments.Phenotype]))
            Pickle = "_".join(["Phenotype=" + Phenotype[:Phenotype.index(":")],
                               "_".join(sorted(DataTypes)), str(Arguments.FeatureMin), Arguments.CorrectionMethod,
                               "".join(map(str, Arguments.BooleanSets)), "".join(map(str, Arguments.Optimization)),
                               ".Validation" + str(CrossValidation)])
        else:
            Pickle = Arguments.Filename + ".Validation" + str(CrossValidation)
        
        cPickle.dump(Results, open(get_path("MOCA.results") + "/" + Pickle, "wb"), -1)

    return

def leave_some_out(Arguments):
    '''
    Makes the data splits before sending data off for cross validation.
    Every label gets used once. To the extent possible, try and put the same 
    number of controls in every split and cases in every split. And, try and 
    balance the number of cases and controls in each individual split. Example:

    Split 1: Cancer10, Cancer28, Cancer15, Healthy2 Healthy3
    Split 2: Cancer2, Cancer12, Cancer1, Healthy10, Healthy1
    Split 3: Cancer1, Cancer3, Healthy9, Healthy0
    etc
    etc
    '''

    Labels, Features, Variates, Phenotypes, Markers = get_supervised_dataset(Arguments)

    #Only one phenotype at a time for cross validation
    Phenotype = Phenotypes[0]

    #Important: if you are in Multiprocess mode you need to set the seed or the data splits won't make sense!!!
    CrossValidations = int(ceil(len(Labels)/float(Arguments.LeaveSomeOut)))

    #Get the cases ("1"s) and controls ("0"s) from the Phenotype vector
    Cases = [Label for Label in Labels if Variates[Features.index(Phenotype)][Labels.index(Label)]]
    Controls = [Label for Label in Labels if not Variates[Features.index(Phenotype)][Labels.index(Label)]]

    shuffle(Cases) #Get rid of bias that MIGHT be inherent in the original data structure
    #Split as evenly as possible among the cross-validations
    Cases = dict([(Iteration, Cases[LabelRange:len(Cases):CrossValidations]) \
                      for Iteration, LabelRange in enumerate(range(CrossValidations))])
    
    shuffle(Controls) #Get rid of bias that MIGHT be inherent in the original data structure
    #Split as evenly as possible among the cross-validations
    Controls = dict([(Iteration, Controls[LabelRange:len(Controls):CrossValidations]) \
                         for Iteration, LabelRange in enumerate(range(CrossValidations))])
    
    #only for the leave ONE out case do we do cases and controls in series
    if Arguments.LeaveSomeOut == 1:
        if len(Cases) < len(Controls):
            Cases = dict(zip(Cases.keys(), list(reversed(Cases.values()))))
        else:
            Controls = dict(zip(Controls.keys(), list(reversed(Controls.values()))))

    cross_validation(Arguments, Labels, Features, Variates, \
                         Markers, Phenotype, CrossValidations, Cases, Controls)
    
    return 
