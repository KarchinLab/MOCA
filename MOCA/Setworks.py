'''
Get python modules
'''
from collections import Counter
from random import sample, shuffle, choice, random
from itertools import chain
import string
import cPickle

'''
Get third-party modules
'''
from fisher import pvalue as fisher

'''
Get MOCA modules
'''
from DataHandler import get_supervised_dataset, get_path, get_priors
from Statistics import combine_features, MOCASet, p_adjust, contingency_table, Performance, weighted_sample, \
    minimum_performance, interaction, rank, EffectSize
    
def reduced_features(Setwork):
    '''
    Make sure you don't end up with two of the same features in a single setwork, 
    including features that only differ in Z-score threshold
    '''

    return [Feature.replace(":", " ").replace("<", " ").replace(">", " ").split()[0] for Feature in list(chain(*Setwork))]
    #return [Feature.replace(".", " ").split()[0] for Feature in list(chain(*Setwork))]

def priors(Priors, UnionCombinations, IntersectionCombinations, DifferenceCombinations, Arguments):
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

def assemble_setwork(Features, Variates, UnionCombination, IntersectionCombination, DifferenceCombination,
                     Arguments):
    '''
    Given the features for each of the Boolean set operations, make the Feature vector
    '''
    
    VariateCombination = [Variates[Features.index(Feature)] for Feature in UnionCombination]
    moca_set = MOCASet(VariateCombination, NA=Arguments.NA)
    Union = moca_set.union
                
    VariateCombination = [Variates[Features.index(Feature)] for Feature in IntersectionCombination]
    VariateCombination.append(Union)
    #None to get rid of empty lists if there were no union features
    moca_set = MOCASet(filter(None, VariateCombination), NA=Arguments.NA)
    Intersection = moca_set.intersection
            
    if Intersection: #if not, this will ultimitely result in returning "None" (see below)
        VariateCombination = [Variates[Features.index(Feature)] for Feature in DifferenceCombination]
        VariateCombination.insert(0, Intersection)
        moca_set = MOCASet(VariateCombination, NA=Arguments.NA)
        FeatureVector = moca_set.difference 

        return FeatureVector
        
    else:
        return None #If all we had was difference features, there is nothing to do

def barcode(Length=10):
    '''
    Create a unique alphanumeric code so that every interaction can be uniquely identified
    '''

    return ''.join(choice(string.ascii_uppercase + string.digits) for _ in range(Length))

def eject_passengers(LocalInteractions, Setworks, PValues, Arguments):
    '''
    Setworks yield passengers. Passengers are features not associated with the phenotype, but happen to get selected
    for because they are combined with a feature that is associated with the phenotype. For any two setworks having at
    least one feature in common, this function makes sure that setwork with more features is more significantly associated
    with the phenotype (i.e., if a feature was added to a setwork, that new setwork should be even more associated with 
    the phenotype then the previous setwork; if it's not, (r)eject the new feature as a passenger). 
    '''

    EjectPassengers = []
    for Barcode1 in LocalInteractions:
        for Barcode2 in LocalInteractions:
            Setwork1 = list(chain(*Setworks[Barcode1]))
            Setwork2 = list(chain(*Setworks[Barcode2]))

            #The two setworks underconsideration need to have at least one feature in common and not be the same length
            if set.intersection(set(Setwork1), set(Setwork2)) and len(set(Setwork1)) != len(set(Setwork2)):

                #We need the info associated with B0, S0, and P0 to correspond to the small setwork
                if len(set(Setwork1)) < len(set(Setwork2)):
                    
                    B0 = Barcode1
                    S0 = set(Setwork1)
                    P0 = PValues[Barcode1]
                    
                    B1 = Barcode2
                    S1 = set(Setwork2)
                    P1 = PValues[Barcode2]
                    
                elif len(set(Setwork2)) < len(set(Setwork1)):
                    
                    B0 = Barcode2
                    S0 = set(Setwork2)
                    P0 = PValues[Barcode2]
                    
                    B1 = Barcode1
                    S1 = set(Setwork1)
                    P1 = PValues[Barcode1]

                '''
                The default of '10' is because you are setting how many orders of magnitude more
                significant you need the P-value in order to accept the additional feature into 
                the growing setwork. 

                'Ndiff' is the how many new features are being added to the setwork. The more additional
                features, the more significant the new setwork needs to be in order to accept the new
                features. 

                Arguments.EjectPassengers is set to 0.1 by defualt. This means that the addition of a 
                single new feature requires a p-value that simply doesn't increase. This is neither 
                liberal or conservative, and could probably be increased to as much as 1.0 in most settings.

                if Ndiff == 1 and Agruments.EjectPassengers == 0.1 (the default), then the following 
                equation says if P0 <= P1, then eject the new feature as a passenger. 
                '''
                
                Ndiff = len(set.difference(S1, S0))
                if P0/(10*Ndiff*Arguments.EjectPassengers) <= P1:
                        EjectPassengers.append(B1)
                #else:
                #    print "~~~~~~~~~~~~~~~~~~~~~~"
                #    print Setwork1, PValues[Barcode1]
                #    print Setwork2, PValues[Barcode2]
                #    print Ndiff, Arguments.EjectPassengers, 10*Ndiff*Arguments.EjectPassengers, P0/(10*Ndiff*Arguments.EjectPassengers)
                #    print 
                    
    return EjectPassengers

def make_report(Labels, Phenotype, Barcodes, Arguments):
    '''
    This report is retrieved for all MOCA results files when 'Reports = Results' is set in the 
    Arguments file. Tells you from which input data the results were generated, which file the 
    results were written to, and what parameters were used during the calculation. 
    '''

    Report = {}
    #Use the 'Entries' key to set the order of output when 'Reports = Results' is run.
    Report["Entries"] = ["Input data", "Label count", "Phenotype", "Number of significant interactions",
                         "Correction method", "FDR", "Optimization parameters", "Boolean sets", "Eject passengers",
                         "Bandwidth", "Seed", "Phenotype permutation test?"]
    Report["Input data"] = [Data for Data in Arguments.Data]
    Report["Label count"] = len(Labels)
    Report["Correction method"] = Arguments.CorrectionMethod
    Report["Phenotype"] = Phenotype[:Phenotype.index(":")]
    Report["Number of significant interactions"] = len(Barcodes)
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

def get_setworks(Arguments, \
                     Features, Variates, Phenotype, \
                     UnionMarkers, IntersectionMarkers, DifferenceMarkers, \
                     Trials, RepopulateFrequency, PercentToRepopulate, \
                     UnionFeatures, IntersectionFeatures, DifferenceFeatures):
    '''
    The core engine for building MOCA setworks (networks of features combined using 
    Boolean set operations). Very customizable. You can build your own implementation 
    using MyMOCA.py, or you can run the simple default mode by calling the setworks function 
    via the Arguments file or the command line (Setworks = True and --setworks True, respectively).

    Phenotype is the thing your selecting markers for
    Markers of each Boolean type represent the initial pool for that type
    Trials, RepopulateFrequency, PercentToRepopulate = Optimization parameters (see arguments or UsersManual)
    UnionFeatures, IntersectionFeatures, DifferenceFeatures = max number of each operation in a single comparison.
    '''

    PValues = {}
    Performances = {}
    Interactions = {}
    FeatureVectors = {}
    Setworks = {}
    SampleCounts = {}
    CaseCounts = {} #just the postive class here
    EffectSizes = {}
    
    #Do this outside of the for loop to prevent re-reading
    if Arguments.Priors: Priors = get_priors(Arguments)

    #Get response outside of the loop, incase we want to permute the phenotype
    Response = Variates[Features.index(Phenotype)]
    if Arguments.PermutePhenotype: shuffle(Response)
    
    for Trial in range(Trials):
        if Trial and not Trial%RepopulateFrequency:
            for Barcode in rank(Performances, Arguments.RankMethod, int(len(Performances.keys())*PercentToRepopulate)):
                Marker = Setworks[Barcode]
                UnionMarkers.extend(Marker[0])
                IntersectionMarkers.extend(Marker[1])
                DifferenceMarkers.extend(Marker[2])
        
        #Define the min and max number of features to combine via each Boolean set operation 
        UnionSample = weighted_sample(UnionMarkers, min(UnionFeatures, len(UnionMarkers)))
        UnionCombinations = combine_features(UnionSample, 0, len(UnionSample))
        
        IntersectionSample = weighted_sample(IntersectionMarkers, min(IntersectionFeatures, len(IntersectionMarkers)))
        IntersectionCombinations = combine_features(IntersectionSample, 0, len(IntersectionSample))
        
        DifferenceSample = weighted_sample(DifferenceMarkers, min(DifferenceFeatures, len(DifferenceMarkers)))
        DifferenceCombinations = combine_features(DifferenceSample, 0, len(DifferenceSample))
            
        if Arguments.Priors: #Are we using priors?
            UnionCombinations, IntersectionCombinations, DifferenceCombinations = \
                priors(Priors, UnionCombinations, IntersectionCombinations, DifferenceCombinations, Arguments)

        LocalInteractions = [] #We'll use this to eject passengers at the end of each trial
        for UnionCombination in UnionCombinations:
            for IntersectionCombination in IntersectionCombinations:
                for DifferenceCombination in DifferenceCombinations:
                    
                    FeatureVector = assemble_setwork(Features, Variates,
                                                     UnionCombination, IntersectionCombination, DifferenceCombination,
                                                     Arguments)
                    if FeatureVector:
                        
                        Predictor = FeatureVector
                        TP,FP,FN,TN = contingency_table(Predictor, Response, NA=Arguments.NA)
                        PValue = fisher(TP,FP,FN,TN)
                        Setwork = [sorted(UnionCombination), sorted(IntersectionCombination), \
                                       sorted(DifferenceCombination)]

                        if Arguments.ForceCooccurring and interaction(PValue) == "MutuallyExclusive": break
                    
                        #We don't want more than one of the same feature in a single setwork unless we are in Bandwidth mode
                        if not Arguments.Bandwidth:
                            if max(Counter(reduced_features(Setwork)).values()) == 1:
                                pass
                            else:
                                break
                            
                        Barcode = barcode()
                        Setworks[Barcode] = tuple(map(tuple, Setwork))
                        PValues[Barcode] = PValue.two_tail
                        SampleCounts[Barcode] = TP + FP + FN + TN
                        CaseCounts[Barcode] = TP + FN
                        Interactions[Barcode] = interaction(PValue)    
                        Performances[Barcode] = Performance(Interactions[Barcode], TP,FP,FN,TN)
                        EffectSizes[Barcode] = EffectSize(Interactions[Barcode], TP,FP,FN,TN)
                        FeatureVectors[Barcode] = FeatureVector
                        LocalInteractions.append(Barcode)

        if Arguments.EjectPassengers:                    
            EjectedPassengers = eject_passengers(LocalInteractions, Setworks, PValues, Arguments)

            for Barcode in EjectedPassengers:
                Setworks.pop(Barcode, None)
                Performances.pop(Barcode, None) #Have to get rid of these to prevent repopulating with passengers
                        
    #Momentarily make the setworks keys so that the python dictionary removes redundancies
    InvertedSetworks = dict(zip(Setworks.values(), Setworks.keys())) 
    Setworks = dict(zip(InvertedSetworks.values(), InvertedSetworks.keys())) #Turn the barcodes back into keys
   
    return PValues, Performances, Interactions, FeatureVectors, Setworks, SampleCounts, CaseCounts, EffectSizes
        
def setworks(Arguments):
    '''
    Default implementation for building the MOCA Boolean set networks (setworks). 
    '''

    Labels, Features, Variates, Phenotypes, Markers = get_supervised_dataset(Arguments)
    
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
        
        PValues, Performances, Interactions, FeatureVectors, Setworks, SampleCounts, CaseCounts, EffectSizes = \
            get_setworks(Arguments, \
                             Features, Variates, Phenotype, \
                             Markers, Markers, Markers, \
                             Trials, RepopulateFrequency, PercentToRepopulate, \
                             UnionFeatures, IntersectionFeatures, DifferenceFeatures)
        
        #We only need the intersection of unique setworks passing the FDR threshold
        QValues = p_adjust(PValues, Arguments.CorrectionMethod)
        QValues = dict([(Barcode, QValues[PValue]) for Barcode, PValue in PValues.items() \
                            if QValues[PValue] < Arguments.FDR])
        Barcodes = list(set.intersection(set(Setworks.keys()), set(QValues.keys())))

        #finally, if we desire we can filter by performance at this stage. We could do it later, but we'll get a bigger Pickle now. 
        Barcodes = [Barcode for Barcode in Barcodes if minimum_performance(Performances[Barcode], Arguments)]

        if Arguments.PermutePhenotype:
            try:
                QValue = min(QValues.values())
                print "Permutation test failed: you ran with 'PermutePhenotype = True' and setworks could be generated that passed your filters!!!",
                print "This means that your current FDR cutoff is not sufficient for this data. The minimum FDR observed during",
                print "this permutation test was " + str(QValue) + ". You should do this a minimum of 10 times and set your FDR",
                print "threshold (i.e., 'FDR = threshold' in your Arguments file) AT LEAST one order of magnitude lower than the",
                print "lowest observed during permutation testing. This conservative threshold will help ensure that results",
                print "observed during your 'real' setworks run are statisically reliable. The setworks that passed your filters",
                print "for this permutation testing have been saved; if you care to see what features made it thru you can use",
                print "the standard 'Mode = PostProcess' to veiw them. Exiting..."
                
            except ValueError:
                print "You ran with 'PermutePhenotype = True' and no setworks could be generated that passed your filters --",
                print "this is a great start! You should do this a minimum of 10 times and set your FDR threshold (i.e., 'FDR = threshold'",
                print "in your Arguments file) AT LEAST one order of magnitude lower than the lowest observed during permutation testing." ,
                print "This conservative threshold will help ensure that results observed during your 'real' setworks run are statisically",
                print "reliable. Exiting..."
                exit()                

        if len(Barcodes):
        
            Results = {}
            Results["PValues"] = dict([(Barcode, PValues[Barcode]) for Barcode in Barcodes])
            Results["QValues"] = dict([(Barcode, QValues[Barcode]) for Barcode in Barcodes])
            Results["Performances"] = dict([(Barcode, Performances[Barcode]) for Barcode in Barcodes]) 
            Results["Interactions"] = dict([(Barcode, Interactions[Barcode]) for Barcode in Barcodes])
            Results["FeatureVectors"] = dict([(Barcode, FeatureVectors[Barcode]) for Barcode in Barcodes]) 
            Results["UnionFeatures"] = dict([(Barcode, Setworks[Barcode][0]) for Barcode in Barcodes])
            Results["IntersectionFeatures"] = dict([(Barcode, Setworks[Barcode][1]) for Barcode in Barcodes])
            Results["DifferenceFeatures"] = dict([(Barcode, Setworks[Barcode][2]) for Barcode in Barcodes])
            Results["SampleCounts"] = dict([(Barcode, SampleCounts[Barcode]) for Barcode in Barcodes])
            Results["CaseCounts"] = dict([(Barcode, CaseCounts[Barcode]) for Barcode in Barcodes])
            Results["EffectSizes"] = dict([(Barcode, EffectSizes[Barcode]) for Barcode in Barcodes])
            Results["Report"] = make_report(Labels, Phenotype, Barcodes, Arguments)
            Results["Labels"] = Labels
            Results["Barcodes"] = Barcodes
            Results["Phenotype"] = Variates[Features.index(Phenotype)]
            
            if Arguments.Filename.lower() == "default":
                DataTypes = set(Arguments.Data).difference(set([Arguments.Phenotype]))
                Pickle = "_".join(["Phenotype=" + Phenotype[:Phenotype.index(":")],
                                   "_".join(sorted(DataTypes)), str(Arguments.FeatureMin), Arguments.CorrectionMethod,
                                   "".join(map(str, Arguments.BooleanSets)), "".join(map(str, Arguments.Optimization))])
            else:
                Pickle = Arguments.Filename + "_" + Phenotype[:Phenotype.index(":")]
            
                cPickle.dump(Results, open(get_path("MOCA.results") + "/" + Pickle, "wb"), -1)

        else:
            print "No setworks were generated. This could mean your data set is not sufficiently powered for deriving setworks",
            print "or that you set your filters unreasonably strict. Exiting...",
            
    return
