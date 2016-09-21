'''
Get python modules
'''
import cPickle
from itertools import chain

'''
Get third-party modules
'''
from fisher import pvalue as fisher

'''
Get MOCA modules
'''
from Statistics import contingency_table, p_adjust, interaction, Performance, correlation, correlation_pvalue, EffectSize
from DataHandler import load_data, get_supervised_dataset, get_path

def make_report(Labels, Features, Arguments, Supervised=False):
    '''
    This report is retrieved for all MOCA results files when 'Reports = Results' is set in the 
    Arguments file. Tells you from which input data the results were generated, which file the 
    results were written to, and what parameters were used during the calculation.
    '''

    Report = {}
    #Use the 'Entries' key to set the order of output when 'Reports = Results' is run.
    Report["Entries"] = ["Pairwise", "Input data", "Label count", "Supervised", "Correction method",
                         "FDR", "Number of significant interactions", "Continuous-valued correlation"]
    Report["Input data"] = [Data for Data in Arguments.Data]
    Report["Label count"] = len(Labels)
    Report["Correction method"] = Arguments.CorrectionMethod
    Report["FDR"] = Arguments.FDR
    Report["Pairwise"] = Arguments.Pairwise
    Report["Supervised"] = Supervised #just initialized it
    if Supervised: Report["Supervised"] = [True, "(" + Supervised + ")"] #then overwrite the key if True
    Report["Number of significant interactions"] = len(Features)
    if Arguments.UntransformedData:
        Report["Continuous-valued correlation"] = True
    else:
        Report["Continuous-valued correlation"] = False

    return Report

def unsupervised(Arguments):
    '''
    Pairwise MOCA calculations that are executed if the Phenotype argument is False (the default). Similar to 
    so-called 'supervised' pairwise mode, except that no performance metrics are calculated (sens, spec, PPV, NPV, etc.).
    In unspervised mode, you can compare all inter-datatype pairs for two datatypes, or all intra-datatype pairs for 
    a single datatype. 
    '''

    if len(Arguments.Data) > 2:
        print "Unsupervised pairwise calculations can consider no more that two datatypes at a time."
        print "If you provide only one datatype, all intra-datatype pairs will be considered. If you"
        print "provide two datatypes, all inter-datatype comparisons will be made. Please change the"
        print "'Data = ' field. Exiting..."
        exit()

    Data = load_data(Arguments)

    Features = list(chain(*Data.Transformed.Features.values()))
    Variates = list(chain(*Data.Transformed.Variates.values()))
    if len(Arguments.Data) == 1:
        Features1 = Features
        Features2 = Features

    if len(Arguments.Data) == 2:
        Features1 = Data.Transformed.Features[Arguments.Data[0]]
        Features2 = Data.Transformed.Features[Arguments.Data[1]]

    PValues = {}
    Interactions = {}
    SampleCounts = {}
    CaseCounts = {} #just the positive class here
    Performances = {}
    EffectSizes  = {}
    Tested = []
    for Feature1 in Features1:
        Tested.append(Feature1)
        for Feature2 in Features2:
            if Feature2 not in Tested:
                a,b,c,d = contingency_table(Variates[Features.index(Feature1)], Variates[Features.index(Feature2)],
                                            NA=Arguments.NA)
                PValue = fisher(a,b,c,d)
                PValues[tuple([Feature1, Feature2])] = PValue.two_tail
                Interactions[tuple([Feature1, Feature2])] = interaction(PValue)
                SampleCounts[tuple([Feature1, Feature2])] = a + b + c + d
                CaseCounts[tuple([Feature1, Feature2])] = a + c
                #A placeholder solely to make pairwise post-processing generalizable
                Performances[tuple([Feature1, Feature2])] = "NA"
                EffectSizes[tuple([Feature1, Feature2])] = "NA"
                
    FDRs = p_adjust(PValues, Arguments.CorrectionMethod)
    for Pair, PValue in PValues.items():
        if FDRs[PValue] < Arguments.FDR:
            pass
        else:
            PValues.pop(Pair, None)
            Interactions.pop(Pair, None)
            SampleCounts.pop(Pair, None)
            CaseCounts.pop(Pair, None)
            Performances.pop(Pair, None)
            EffectSizes.pop(Pair, None)

    Results = {}
    Results["Report"] = make_report(Data.Labels, PValues.keys(), Arguments) 
    Results["PValues"] = PValues
    Results["Interactions"] = Interactions
    Results["FDRs"] = FDRs
    Results["SampleCounts"] = SampleCounts
    Results["CaseCounts"] = CaseCounts
    Results["Performances"] = Performances
    Results["EffectSizes"] = EffectSizes

    if Arguments.Filename.lower() == "default":
        Pickle = "_".join(["Pairwise", "_".join(sorted(Arguments.Data)), str(Arguments.FeatureMin),
                           Arguments.CorrectionMethod])
    else:
        Pickle = Arguments.Filename

    cPickle.dump(Results, open(get_path("MOCA.results") + "/" + Pickle, "wb"), -1)
        
    return

def supervised(Arguments):
    '''
    MOCA pairwise calculations executed if a 'Phenotype' is provided in the Arguments file. Not technically
    supervised 'learning', as there is no optimization (every possible pairwise comparison is tested). 
    Output includes perfomance metrics such as sensitivity, specificity, PPV, and NPV, for each features 
    ability to predict the phenotype. 
    '''
    
    Labels, Features, Variates, Phenotypes, Markers = get_supervised_dataset(Arguments)

    #Clustermode support if called
    Node = int(Arguments.MultiProcessMode[0]) 
    TotalNodes = int(Arguments.MultiProcessMode[1]) 

    for Phenotype in Phenotypes[Node::TotalNodes]:
        PValues = {}
        Interactions = {}
        Performances = {}
        SampleCounts = {}
        CaseCounts = {} #just the postive class here
        EffectSizes  = {}
        for Marker in Markers:
            TP,FP,FN,TN = contingency_table(Variates[Features.index(Marker)], Variates[Features.index(Phenotype)],
                                            NA=Arguments.NA)
            PValue = fisher(TP,FP,FN,TN)
            PValues[Marker] = PValue.two_tail
            Interaction = interaction(PValue)
            Interactions[Marker] = Interaction
            Performances[Marker] = Performance(Interaction, TP,FP,FN,TN)
            EffectSizes[Marker] = EffectSize(Interaction, TP,FP,FN,TN)
            SampleCounts[Marker] = TP + FP + FN + TN
            CaseCounts[Marker] = TP + FN

        FDRs = p_adjust(PValues, Arguments.CorrectionMethod)
        for Marker in Markers:
            FDR = FDRs[PValues[Marker]]
            if FDR < Arguments.FDR:
                pass
            else:
                PValues.pop(Marker, None)
                Interactions.pop(Marker, None)
                Performances.pop(Marker, None)
                SampleCounts.pop(Marker, None)
                CaseCounts.pop(Marker, None)
                EffectSizes.pop(Marker, None)

        if len(PValues.keys()):
            Results = {}
            Results["Report"] = make_report(Labels, PValues.keys(), Arguments, Supervised=Phenotype[:Phenotype.index(":")])
            Results["PValues"] = PValues
            Results["Interactions"] = Interactions
            Results["Performances"] = Performances
            Results["FDRs"] = FDRs
            Results["SampleCounts"] = SampleCounts
            Results["CaseCounts"] = CaseCounts
            Results["EffectSizes"] = EffectSizes
            
            if Arguments.Filename.lower() == "default":
                DataTypes = set(Arguments.Data).difference(set([Arguments.Phenotype]))
                Pickle = "_".join(["Pairwise", "Phenotype=" + Phenotype[:Phenotype.index(":")],
                                   "_".join(sorted(DataTypes)), str(Arguments.FeatureMin), Arguments.CorrectionMethod])
            else:
                Pickle = Arguments.Filename + "_" + Phenotype[:Phenotype.index(":")]
                
            cPickle.dump(Results, open(get_path("MOCA.results") + "/" + Pickle, "wb"), -1)
        
    return

def pairwise_continuous(Arguments):
    '''
    '''

    if len(Arguments.Data) > 2:
        print "Unsupervised pairwise calculations can consider no more that two datatypes at a time."
        print "If you provide only one datatype, all intra-datatype pairs will be considered. If you"
        print "provide two datatypes, all inter-datatype comparisons will be made. Please change the"
        print "'Data = ' field. Exiting..."
        exit()

    Data = load_data(Arguments)

    Features = list(chain(*Data.Features.values()))
    Variates = list(chain(*Data.Variates.values()))

    if Arguments.Phenotype:
        Features1 = [Feature for Feature in Features if Arguments.Phenotype in Feature]
        Features2 = [Feature for Feature in Features if Arguments.Phenotype not in Feature]

    else:

        if len(Arguments.Data) == 1:
            Features1 = Features
            Features2 = Features

        if len(Arguments.Data) == 2:
            Features1 = Data.Features[Arguments.Data[0]]
            Features2 = Data.Features[Arguments.Data[1]]

    PValues = {}
    Correlations = {}
    Tested = []
    for Feature1 in Features1:
        Tested.append(Feature1)
        for Feature2 in Features2:
            if Feature2 not in Tested:
                PValues[tuple([Feature1, Feature2])] = correlation_pvalue(Variates[Features.index(Feature1)],
                                                                          Variates[Features.index(Feature2)])
                Correlations[tuple([Feature1, Feature2])] = correlation(Variates[Features.index(Feature1)],
                                                                       Variates[Features.index(Feature2)])
    
    FDRs = p_adjust(PValues, Arguments.CorrectionMethod)
    for Pair, PValue in PValues.items():
        if FDRs[PValue] < Arguments.FDR:
            pass
        else:
            PValues.pop(Pair, None)
            Correlations.pop(Pair, None)

    if len(PValues.keys()):
        Results = {}
        Results["Report"] = make_report(Data.Labels, PValues.keys(), Arguments, Supervised=Arguments.Phenotype)
        Results["PValues"] = PValues
        Results["Correlations"] = Correlations
        Results["FDRs"] = FDRs

        if Arguments.Filename.lower() == "default":
            Pickle = "_".join(["_".join(sorted(Arguments.Data)), Arguments.CorrectionMethod])
        else:
            Pickle = Arguments.Filename
                
        cPickle.dump(Results, open(get_path("MOCA.results") + "/" + Pickle, "wb"), -1)
    
    return
        
def pairwise(Arguments):
    '''
    If you supply a phenotype, then pairwise MOCA is simply selecting features that
    are significantly co-occurring or mutually exclusive with the phenotype. If no
    phenotype is specified, MOCA will compare every possible pairwise interaction in the 
    provided data. 
    '''

    if Arguments.UntransformedData:
        pairwise_continuous(Arguments)
    elif Arguments.Phenotype:
        supervised(Arguments)
    else:
        unsupervised(Arguments)
    
    return
    
