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
from GetData import get_master_matrix, get_ordered_matrix
from Setworks import get_setworks, assemble_setwork, sens_and_spec, balanced_accuracy
from Statistics import p_adjust, MOCASet, contingency_table, minimum_performance

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
        EmpiricalDistribution, PredictiveValues, Interaction, FeatureVectors, \
            UnionMarkers, IntersectionMarkers, DifferenceMarkers = \
            get_setworks(Arguments, \
                             Features, TrainVariates, Phenotype, \
                             Markers, Markers, Markers, \
                             Trials, RepopulateFrequency, PercentToRepopulate, \
                             UnionFeatures, IntersectionFeatures, DifferenceFeatures, \
                             sens_and_spec, balanced_accuracy)
        
        #Make the cross-validation matrix. First row is the case-control-label header
        CrossValidationFeatureMatrix = [Cases[CrossValidation] + Controls[CrossValidation]]
        PhenotypeVector = [1 for Case in Cases[CrossValidation]] + [0 for Control in Controls[CrossValidation]]
        CrossValidationFeatureMatrix.append(PhenotypeVector) #second row is the phenotype vector

        FDRs = p_adjust(EmpiricalDistribution, Arguments.CorrectionMethod) #Correct the P-values
        for Setwork, PValue in EmpiricalDistribution.items():
            FDR = FDRs[PValue]
            #See if Setwork passed either the FDR cutoff and/or minimum user-defined predictive performance requirements
            if FDR < Arguments.FDR and minimum_performance(Arguments, Interaction[Setwork], \
                                                               TrainVariates[Features.index(Phenotype)], FeatureVectors[Setwork]):
                UnionCombination, IntersectionCombination, DifferenceCombination = Setwork #Decompose setwork into components
                #We need to assemble the setwork from the full-length vectors to subsequently grab the cases and controls
                FeatureVector = assemble_setwork(Features, Variates, \
                                                     UnionCombination, IntersectionCombination, DifferenceCombination)
                #From the full-length setwork vector, get the elements for the cases and controls
                CrossValidationFeatureVector = \
                    [FeatureVector[Labels.index(Case)] for Case in Cases[CrossValidation]] + \
                    [FeatureVector[Labels.index(Control)] for Control in Controls[CrossValidation]] 
                #Setworks are bundled with the interaction type, which is required for post processing
                CrossValidationFeatureMatrix.append([(Setwork, Interaction[Setwork]), CrossValidationFeatureVector])
            
        #Pickle each validation seperately
        cPickle.dump(CrossValidationFeatureMatrix, \
                         open(Arguments.Filename + ".Validation." + str(CrossValidation), "wb"), -1)
    
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

    Data = get_master_matrix(Arguments)
    Labels = Data["Labels"]
    Features = Data["Features"]
    Variates = Data["Variates"]

    #Only one phenotype at a time for cross validation
    Phenotype = [Feature for Feature in Features if Arguments.Phenotype in Feature][0]
    Markers = [Feature for Feature in Features if Arguments.Phenotype not in Feature]

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
