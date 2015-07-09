'''
Get python modules
'''
from collections import Counter
from random import sample, shuffle
from itertools import chain

'''
Get third-party modules
'''
from fisher import pvalue as fisher

'''
Get MOCA modules
'''
from GetData import get_master_matrix, get_priors
from Statistics import combine_features, MOCASet, p_adjust, contingency_table, Performance, minimum_performance, \
    weighted_sample

def balanced_accuracy(List):
    '''
    The default function for ranking the "performance" of the setworks for predicting the 
    phenotype. Given a list of two-item tuples, where the first item is the setwork and 
    the second item is itself a two-item tuple of the sensitivity and specificity, the 
    function returns the list ranked by balance accuracy (arithmetic mean of sensitivity
    and specificity).
    '''

    return sorted(List, key=lambda Item: (Item[1][0] + Item[1][1])/2.0, reverse=True)
    
def sens_and_spec(TP, TN, FP, FN, Interaction):
    '''
    The default "scores" for setworks is the sensitivity and specificity.
    '''
    performance = Performance(TP,TN,FP,FN)

    if Interaction == "MutuallyExclusive": #logic is reversed for mutually exclusive interactions
        return (round(1 - performance.sensitivity, 2), round(1 - performance.specificity, 2))
    else:
        return (round(performance.sensitivity, 2), round(performance.specificity, 2))

def reduced_features(Setwork):
    '''
    Make sure you don't end up with two of the same features in a single setwork, 
    including features that only differ in Z-score threshold
    '''

    return [Feature.replace("Positive", " Positive").replace("Negative", " Negative").split()[0] \
                for Feature in list(chain(*Setwork))]

def priors(Arguments, Priors, UnionCombinations, IntersectionCombinations, DifferenceCombinations):
    '''
    Use priors to focus the search, while making setworks, with feature you think might be important, 
    or those selected from previous Setwork runs. Three different possible was to utilize priors:
    '''
    
    def strict(Priors, Combinations):
        '''
        If you say include Gene1 and Gene2 in the union, they are included in every single Setwork.
        Specicify as many genes as you want, and for which of the Boolean set operations you want those
        features included
        '''

        PriorCombinations = []
        for Prior, SetCombination in zip(Priors, Combinations):
            PriorCombinations.append([Prior + Combination for Combination in SetCombination])

        return PriorCombinations
    
    def weighted(Priors, Combinations):
        '''
        You'll probably only use this one if you made your Priors file using MakePriors. Because this
        occurs after an initial Setworks or ValidateBiomakers run, you could have, say 50 copies of Feature1
        and 10 copies of Feature2. Then, when randomly assigning a prior to a setwork, there is 5 times the 
        chance that Feature1 would be included, relative to Feature2. The number of features to be included
        is a randomly choosen number between 1 and 10, where selecting a single feature is 10 times 
        more likely than choosing 10...choosing two features is five times as likey as choosing 10, etc.
        '''

        #Gimme a list of numbers, used to randomly determine how many features will be selected for a single 
        #Boolean set operation. The list is biased toward choosing small numbers
        Sample = list(chain(*[range(1, Integer + 1) for Integer in range(1, 11)]))
        shuffle(Sample)

        #Randomly select the priors to use for each Boolean set operation. The number to include is random, too.
        UnionPriors = tuple(weighted_sample(Priors[0], min(sample(Sample, 1)[0], len(set(Priors[0])))))
        IntersectionPriors = tuple(weighted_sample(Priors[1], min(sample(Sample, 1)[0], len(set(Priors[1])))))
        DifferencePriors = tuple(weighted_sample(Priors[2], min(sample(Sample, 1)[0], len(set(Priors[2])))))

        PriorCombinations = []
        Priors = [UnionPriors, IntersectionPriors, DifferencePriors]
        for Prior, SetCombination in zip(Priors, Combinations):
            PriorCombinations.append([Prior + Combination for Combination in SetCombination])

        return PriorCombinations

    def stochastic(Priors, Combinations):
        '''
        Same as weighted. Except that list of priors are uniquified prior to selection. So if you priors file has
        50 copies of Feature1 and 10 copies of Feature2, there is still an equal likelyhood that either Feature1 or 
        Feature2 will be randomly selected. 
        '''

        Sample = list(chain(*[range(1, Integer + 1) for Integer in range(1, 11)]))
        shuffle(Sample)

        UnionPriors = tuple(sample(set(Priors[0]), min(sample(Sample, 1)[0], len(set(Priors[0])))))
        IntersectionPriors = tuple(sample(set(Priors[1]), min(sample(Sample, 1)[0], len(set(Priors[1])))))
        DifferencePriors = tuple(sample(set(Priors[2]), min(sample(Sample, 1)[0], len(set(Priors[2])))))

        PriorCombinations = []
        Priors = [UnionPriors, IntersectionPriors, DifferencePriors]
        for Prior, SetCombination in zip(Priors, Combinations):
            PriorCombinations.append([Prior + Combination for Combination in SetCombination])

        return PriorCombinations


    if Arguments.Priors.lower() == "strict": 
        return strict(Priors, [UnionCombinations, IntersectionCombinations, DifferenceCombinations])
    elif Arguments.Priors.lower() == "weighted": 
        return weighted(Priors, [UnionCombinations, IntersectionCombinations, DifferenceCombinations])
    elif Arguments.Priors.lower() == "stochastic": 
        return stochastic(Priors, [UnionCombinations, IntersectionCombinations, DifferenceCombinations])
    else:
        print "You called Priors but didn't specify 'strict', 'weighted', or 'stochastic'"
        print "Please correct your commands and try again"
        print "exiting...."
        exit()

def assemble_setwork(Features, Variates, UnionCombination, IntersectionCombination, DifferenceCombination):
    '''
    Given the features for each of the Boolean set operations, make the Feature vector
    '''
    
    VariateCombination = [Variates[Features.index(Feature)] for Feature in UnionCombination]
    moca_set = MOCASet(VariateCombination)
    Union = moca_set.union
                
    VariateCombination = [Variates[Features.index(Feature)] for Feature in IntersectionCombination]
    VariateCombination.append(Union)
    #None to get rid of empty lists if there were no union features
    moca_set = MOCASet(filter(None, VariateCombination))
    Intersection = moca_set.intersection
            
    if Intersection: #if not, this will ultimitely result in returning "None" (see below)
        VariateCombination = [Variates[Features.index(Feature)] for Feature in DifferenceCombination]
        VariateCombination.insert(0, Intersection)
        moca_set = MOCASet(VariateCombination)
        FeatureVector = moca_set.difference 

        return FeatureVector
        
    else:
        return None #If all we had was difference features, there is nothing to do

def get_setworks(Arguments, \
                     Features, Variates, Phenotype, \
                     UnionMarkers, IntersectionMarkers, DifferenceMarkers, \
                     Trials, RepopulateFrequency, PercentToRepopulate, \
                     UnionFeatures, IntersectionFeatures, DifferenceFeatures, \
                     score_function, evaluation_method):
    '''
    The core engine for building MOCA setworks (networks of features combined using 
    Boolean set operations). Very customizable. You can build your own implementation 
    using MyMOCA.py, or you can run the simple default mode by calling the setworks function 
    via the Arguments file or the command line (Setworks = True and --setworks True, respectively).

    Phenotype is the thing your selecting markers for
    Markers of each Boolean type represent the initial pool for that type
    Trials, RepopulateFrequency, PercentToRepopulate = Optimization parameters (see arguments or UsersManual)
    UnionFeatures, IntersectionFeatures, DifferenceFeatures = max number of each operation in a single comparison
    score_function, evaluation_method = see above-defined functions. You can also creat your own in MyMOCA.py
    '''

    EmpiricalDistribution = {}
    Scores = {}
    Interaction = {}
    FeatureVectors = {}

    #Do this outside of the for loop to prevent re-reading
    if Arguments.Priors: Priors = get_priors(Arguments)

    for Trial in range(Trials):
        if Trial and not Trial%RepopulateFrequency:
            for Marker, PredictiveValue in evaluation_method(Scores.items()) \
                    [:int(len(Scores.keys())*PercentToRepopulate)]:
                UnionMarkers.extend(Marker[0])
                IntersectionMarkers.extend(Marker[1])
                DifferenceMarkers.extend(Marker[2])
                
        #Define the min and max number of features to combine via each Boolean set operation 
        UnionSample = sample(UnionMarkers, min(UnionFeatures, len(UnionMarkers)))
        UnionCombinations = combine_features(UnionSample, 0, len(UnionSample))
        
        IntersectionSample = sample(IntersectionMarkers, min(IntersectionFeatures, len(IntersectionMarkers)))
        IntersectionCombinations = combine_features(IntersectionSample, 0, len(IntersectionSample))
        
        DifferenceSample = sample(DifferenceMarkers, min(DifferenceFeatures, len(DifferenceMarkers)))
        DifferenceCombinations = combine_features(DifferenceSample, 0, len(DifferenceSample))

        if Arguments.Priors: #Are we using priors?
            UnionCombinations, IntersectionCombinations, DifferenceCombinations = \
                priors(Arguments, Priors, UnionCombinations, IntersectionCombinations, DifferenceCombinations)

        for UnionCombination in UnionCombinations:
            for IntersectionCombination in IntersectionCombinations:
                for DifferenceCombination in DifferenceCombinations:
                    FeatureVector = assemble_setwork(Features, Variates, \
                                                         UnionCombination, IntersectionCombination, DifferenceCombination)
                    if FeatureVector:
                        Response, Predictor = zip(*[Pair for Pair in zip(Variates[Features.index(Phenotype)], FeatureVector)])
                        TP,FN,FP,TN = contingency_table(Response, Predictor)
                        PValue = fisher(TP,FN,FP,TN)
                        Setwork = [sorted(UnionCombination), sorted(IntersectionCombination), \
                                       sorted(DifferenceCombination)]
                        
                        #We don't want more than one of the same feature in a single setwork
                        if max(Counter(reduced_features(Setwork)).values()) == 1:
                            
                            Setwork = tuple(map(tuple, Setwork))
                            EmpiricalDistribution[Setwork] = PValue.two_tail
                            
                            if PValue.left_tail < PValue.right_tail:
                                Interaction[Setwork] = "MutuallyExclusive"
                            else:
                                Interaction[Setwork] = "Co-occurring"

                            Scores[Setwork] = score_function(TP, TN, FP, FN, Interaction[Setwork])
                            FeatureVectors[Setwork] = FeatureVector

    return EmpiricalDistribution, Scores, Interaction, FeatureVectors, \
        UnionMarkers, IntersectionMarkers, DifferenceMarkers
        

def setworks(Arguments):
    '''
    Default implementation for building the MOCA Boolean set networks (setworks). 
    '''

    Data = get_master_matrix(Arguments)
    Labels = Data["Labels"]
    Features = Data["Features"]
    Variates = Data["Variates"]

    Phenotypes = [Feature for Feature in Features if Arguments.Phenotype in Feature]
    Markers = [Feature for Feature in Features if Arguments.Phenotype not in Feature]

    Trials = int(Arguments.Optimization[0])
    RepopulateFrequency = int(Arguments.Optimization[1])
    PercentToRepopulate = float(Arguments.Optimization[2])

    UnionFeatures = int(Arguments.BooleanSets[0])
    IntersectionFeatures = int(Arguments.BooleanSets[1])
    DifferenceFeatures = int(Arguments.BooleanSets[2])

    #MultiProcessMode support if called. Each Phenotype gets its own node
    Node = int(Arguments.MultiProcessMode[0]) 
    TotalNodes = int(Arguments.MultiProcessMode[1]) 

    for Phenotype in Phenotypes[Node::TotalNodes]:
        EmpiricalDistribution, PredictiveValues, Interaction, FeatureVectors, \
            UnionMarkers, IntersectionMarkers, DifferenceMarkers = \
            get_setworks(Arguments, \
                             Features, Variates, Phenotype, \
                             Markers, Markers, Markers, \
                             Trials, RepopulateFrequency, PercentToRepopulate, \
                             UnionFeatures, IntersectionFeatures, DifferenceFeatures, \
                             sens_and_spec, balanced_accuracy)

        FileOut = open(Arguments.Filename + "-" + Phenotype + ".Setworks", "w")
        FDRs = p_adjust(EmpiricalDistribution, Arguments.CorrectionMethod) #Correct the P-Values
        for Setwork, PredictiveValue in sorted(PredictiveValues.items(), \
                                                   key=lambda Item: Item[1][0] + Item[1][1], reverse=True):
            PValue = EmpiricalDistribution[Setwork]
            FDR = FDRs[PValue]
            #See if Setwork passed either the FDR cutoff and/or minimum user-defined predictive performance requirements
            if FDR < Arguments.FDR and minimum_performance(Arguments, Interaction[Setwork], \
                                                               Variates[Features.index(Phenotype)], FeatureVectors[Setwork]):
                FileOut.write("%s %s %0.2e %0.2e %s %s \n" \
                                  %(Setwork, "\t", PValue, FDR, " ".join(map(str, PredictiveValue)), Interaction[Setwork]))
        FileOut.close()
    
    return
