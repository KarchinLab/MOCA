'''
Get python modules
'''
import cPickle
from collections import Counter
from itertools import chain

'''
Get third-party modules
'''

'''
Get MOCA modules
'''
from Source.SystemCommands import ls
from Source.GetData import get_file
from Source.Statistics import Performance, contingency_table

def load_validation_data(Arguments):
    '''
    Load the python pickles you made from LeaveSomeOut cross-validation calculations
    '''

    CrossValidations = [File for File in ls("./") if Arguments.Filename + ".Validation." in File]
    CrossValidations = [cPickle.load(open(Validation, "rb")) for Validation in CrossValidations]

    #Row one is you labels, two is the Phenotype vector, and then all the feature vectors (with interaction type)
    Labels = dict([(CrossValidation, Data[0]) for CrossValidation, Data in enumerate(CrossValidations)])
    Phenotype = dict([(CrossValidation, Data[1]) for CrossValidation, Data in enumerate(CrossValidations)])
    CrossValidations = dict([(CrossValidation, dict(Data[2:])) for CrossValidation, Data in enumerate(CrossValidations)])

    return Labels, Phenotype, CrossValidations

def biomarkers(Arguments):
    ''' 
    Function for selecting biomarkers from cross-validation output. Strict in that it only returns 
    biomarkers that were selected during each cross validation. Might have the benefit, relative 
    to vote-based prediction, that these were so predictive that they will translate better to future 
    predictions. Also lends itself to simple clinical use because they require little or no computational 
    support for subsequent prediction...it is simply the marker
    '''

    Labels, Phenotype, CrossValidations = load_validation_data(Arguments)
    
    Setworks = []
    for Validation in CrossValidations:
        for Setwork in CrossValidations[Validation]:
            Setworks.append(Setwork)
    
    Results = []
    for Setwork, Count in Counter(Setworks).items():
        if Count == len(CrossValidations): #Sework had to be selected in each cross validation!!!
            ContingencyTables = []
            #Collect the Response and Predictors for ALL cross-validations and concatenate
            for Response, Predictor in zip([Phenotype[Validation] for Validation in CrossValidations], \
                                               [CrossValidations[Validation][Setwork] for Validation in CrossValidations]):
                ContingencyTables.append(contingency_table(Response, Predictor))
            TP,FN,FP,TN = [sum(ContingencyTable) for ContingencyTable in zip(*ContingencyTables)]
            performance = Performance(TP,TN,FP,FN)
            Results.append((Setwork[0], round(performance.sensitivity, 2), round(performance.specificity, 2), Setwork[1]))
    
    Biomarkers = open(Arguments.Filename + ".Biomarkers", "w")
    for Setwork, Sensitivity, Specificity, Interaction in \
            sorted(Results, key=lambda Performance: Performance[1] + Performance[2] if Performance[3] == "Co-occurring" \
                       else (1 - Performance[1]) + (1 - Performance[2]), reverse=True):
        if Interaction == "MutuallyExclusive":
            Biomarkers.write("%s %s %0.2f %0.2f %s \n" %(Setwork, "\t", 1 - Sensitivity, 1 - Specificity, Interaction))
        else:
            Biomarkers.write("%s %s %0.2f %0.2f %s \n" %(Setwork, "\t", Sensitivity, Specificity, Interaction))

    Biomarkers.close()
    
    return 

def vote(Arguments):
    '''
    Loads your cross-validation data and lets the markers vote on their corresponding labels. The subset of 
    labels (hold outs) for in a given cross-validation have their own set of markers determined significant 
    when you ran LeaveSomeOut. Here, that set of markers votes on each label. This is carried out for each 
    of your Validation files until all labels are accounted for. Then, by majority vote, each label is predicted 
    as a "0" or "1" phenotype. Finally, those predictions are compared to the true phenotypes, and the sensitivity
    and specificity are computed. This is distinct from ValidateBiomarkers, in that that function only considers
    markers that occurred in EVERY leave-some-out cross-validation.. 
    '''
    
    Labels, Phenotype, CrossValidations = load_validation_data(Arguments)

    for Label in set(chain(*Labels.values())):
        vars()[Label] = []

    for Validation in CrossValidations:
        for Setwork in CrossValidations[Validation]:
            for Label in Labels[Validation]:
                if Setwork[1] == "MutuallyExclusive" and CrossValidations[Validation][Setwork][Labels[Validation].index(Label)]:
                    vars()[Label].append(0)
                elif Setwork[1] == "MutuallyExclusive" and not CrossValidations[Validation][Setwork][Labels[Validation].index(Label)]:
                    vars()[Label].append(1)
                else:
                    vars()[Label].append(CrossValidations[Validation][Setwork][Labels[Validation].index(Label)])
                
    #Zip the Labels with their respective Phenotypes, from each cross-validation, and make into dictionary
    Phenotypes = []
    for Validation in zip(Labels.values(), Phenotype.values()):
        Phenotypes.extend(zip(*Validation))

    Phenotypes = dict(Phenotypes)

    Response, Predictor = [], []
    for Label in set(chain(*Labels.values())):
        #This blocked-out function might be useful for determining how well we did on each sample
        #print sorted(Counter(vars()[Label]).items(), key=lambda x: x[1], reverse=True), \
            #sorted(Counter(vars()[Label]).items(), key=lambda x: x[1], reverse=True)[0][0], Phenotypes[Label]
        
        #Were the "yays" ("1"s) or "nays" ("0"s) the majority vote for this label
        Majority = sorted(Counter(vars()[Label]).items(), key=lambda x: x[1], reverse=True)[0][0]
        Response.append(Phenotypes[Label])
        Predictor.append(Majority)

    TP,FN,FP,TN = contingency_table(Response, Predictor)
    performance = Performance(TP,TN,FP,FN)
    
    print round(performance.sensitivity, 2), round(performance.specificity, 2)

    return

def make_priors(Arguments):
    '''
    This function can make a Priors file from any setwork output (probably from either a Setworks run, or a 
    ValidateBiomarkers run). This Priors file can then be used to bias future searches. Can be usefule for 
    focusing setwork selection from "big data". For instance, do a Setworks run to see which individual features
    are being selected with high frequency. Next, make the Priors file and rerun Setworks forcing MOCA to bias 
    its search around the individual features predicted to be important in your previous run. This can be useful, 
    because some data is too big to search thru with the normal MOCA optimization process, without a little help. 
    This would be that help. 
    '''

    Setworks = [eval(Setwork.split("\t")[0]) for Setwork in file(get_file("Priors", Arguments.PathFile))]
    
    UnionFeatures, IntersectionFeatures, DifferenceFeatures = [], [], []
    for Union, Intersection, Difference in Setworks:
        UnionFeatures.extend(Union)
        IntersectionFeatures.extend(Intersection)
        DifferenceFeatures.extend(Difference)

    Priors = open(Arguments.Filename + ".Priors", "w")
    Priors.write("%s %s \n" %("Union =", " ".join(UnionFeatures)))
    Priors.write("%s %s \n" %("Intersection =", " ".join(IntersectionFeatures)))
    Priors.write("%s %s \n" %("Difference =", " ".join(DifferenceFeatures)))
    Priors.close()

    return 
