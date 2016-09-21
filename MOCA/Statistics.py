'''
Get python modules
'''
from itertools import combinations as combos 
from itertools import permutations as permuts
from math import sqrt
from random import sample

'''
Get third-party modules
'''
from numpy import random, mean, median, std, arange, trapz, matrix, linalg
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

'''
Get MOCA modules
'''
    
def z_score(Vector, Element, NA="NA"):
    '''
    Returns z-score for specified 'Float' element with respect
    to a particular 'Vector' of floats. Sometimes data is missing
    values; you can specify the NA value (default="NA"; e.g., 0.0
    or -99999.9).
    '''

    Vector = map(float, filter(lambda Element: Element != NA, Vector))

    return (float(Element) - mean(Vector))/std(Vector)
    
def z_vector(Vector, NA="NA"):
    '''
    Converts a vector of floats in to a vector of z-scores
    Sometimes data is missing values; you can specify the 
    NA value (default="NA"; e.g., 0.0 or -99999.9).
    '''
    
    return [z_score(filter(lambda Element: Element != NA, Vector), Element) \
                if Element != NA else NA for Element in Vector]

def z_matrix(Matrix, NA="NA"):
    '''
    Converts a matrix of floats in to a matrix of z-scores
    Sometimes data is missing values; you can specify the 
    NA value (default="NA"; e.g., 0.0 or -99999.9).
    '''

    return [z_vector(Vector, NA) for Vector in Matrix]
   
def variance(Vector):
    '''
    Return the variance of vector
    '''
    
    return mean([abs(Element - mean(Vector))**2 for Element in Vector])

def mad(Vector):
    '''
    Return the mean anomolous dispersion of a vector
    '''
    
    return median([abs(Element - median(Vector))**2 for Element in Vector])

def covariance(Vector1, Vector2):
    '''
    Return the covariance of two vectors 
    '''

    N = len(Vector1)
    XBar = sum(Vector1)/N
    YBar = sum(Vector2)/N
    return sum([(X - XBar)*(Y - YBar) for X,Y in zip(Vector1, Vector2)])/(N)

def correlation(Vector1, Vector2, Method="pearson"):
    '''
    Return the two-tailed p-value associated with a pearson correlation coefficient
    '''

    Correlation = ro.r['cor.test'](ro.FloatVector(Vector1), ro.FloatVector(Vector2), alternative='two.sided', method=Method)

    return Correlation[3][0]

def correlation_pvalue(Vector1, Vector2, Method="pearson"):
    '''
    Return the two-tailed p-value associated with a pearson correlation coefficient 
    (or specify another R-compliant correlation Method
    '''

    Correlation = ro.r['cor.test'](ro.FloatVector(Vector1), ro.FloatVector(Vector2), alternative='two.sided', method=Method)

    return Correlation[2][0]

def binary_vector(Vector, CutOff=0.0, GreaterThan=True, NA="NA"):
    '''
    Given a vector, returns a binary vector(1 Element >= CutOff; 
    0 if Element <= CutOff (logic flipped if GreaterThan=False))
    '''

    if GreaterThan: return [(NA if Element == NA else (1 if float(Element) >= CutOff else 0)) for Element in Vector]
    else: return [(NA if Element == NA else (1 if float(Element) <= CutOff else 0)) for Element in Vector]

def binary_matrix(Matrix, CutOff=0.0, GreaterThan=True, NA="NA"):
    '''
    Given a matrix, returns a binary matrix (1 Element >= CutOff; 
    0 if Element <= CutOff (logic flipped if GreaterThan=False))
    '''
    
    return [binary_vector(Vector, CutOff=CutOff, GreaterThan=GreaterThan, NA=NA) \
                for Vector in Matrix]

def linear_regression(X, Y):
    '''
    Provided two vectors of equal length, fit points to the line
    y = ax + b, return a, b, and the residual (RR)
    '''

    if len(X) != len(Y): raise ValueError, "unequal length"
    N = len(X)
    Sx = Sy = Sxx = Syy = Sxy = 0.0
    for x, y in map(None, X, Y):
        Sx = Sx + x
        Sy = Sy + y
        Sxx = Sxx + x*x
        Syy = Syy + y*y
        Sxy = Sxy + x*y
    det = Sxx * N - Sx * Sx
    a, b = (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det
    meanerror = residual = 0.0
    for x, y in map(None, X, Y):
        meanerror = meanerror + (y - Sy/N)**2
        residual = residual + (y - a * x - b)**2
    RR = 1 - residual/meanerror
    ss = residual / (N-2)
    Var_a, Var_b = ss * N / det, ss * Sxx / det
    
    return [a, b, RR]

def contingency_table(Predictor, Response, NA="NA"):
    '''
    Takes two vectors and makes a contingency table (aka
    confusion matrix). Useful for fishers test and calculating 
    performance metrics (specificity, sensitivity, etc.)
    '''

    a,b,c,d = 0,0,0,0
    for Elements in zip(Predictor, Response):
        if NA in Elements: pass
        elif Elements[0] and Elements[1]: a += 1 #AKA, TP (True positive)
        elif Elements[0]: b += 1 #AKA, FP (False positive)
        elif Elements[1]: c += 1 #AKA, FN (False negative)
        else: d += 1 #AKA, TN (True negative)

    return a,b,c,d #AKA, TP,FP,FN,TN

class Performance:
    '''
    Class that computes standard performance metrics for statistical
    tests. Example:

    from Statistics import Performance

    performance = Performance(20, 3, 3, 4)
    print performance.specificity, performance.sensitivity 
    '''
    
    def __init__(self, Interaction, TP=0, FP=0, FN=0, TN=0):

        try: #sensitivity (AKA true positive rate)
            Sensitivity = float(TP) / (TP + FN)
            if Interaction == "MutuallyExclusive": self.sensitivity = 1 - Sensitivity
            elif Interaction == "Co-occurring": self.sensitivity = Sensitivity
            else: 
                print "Developer error: 'Interaction' must be defined as either 'Co-occurring' or 'MutuallyExclusive'"
                print "Code died in Stastics:Performance, exiting"
                exit()
        except ZeroDivisionError:
            self.sensitivity = 0.0

        try: #specificity (AKA 1 - false positive rate)
            Specificity = float(TN) / (FP + TN)
            if Interaction == "MutuallyExclusive": self.specificity = 1 - Specificity
            elif Interaction == "Co-occurring": self.specificity = Specificity
            else: 
                print "Developer error: 'Interaction' must be defined as either 'Co-occurring' or 'MutuallyExclusive'"
                print "Code died in Stastics:Performance, exiting"
                exit()
        except ZeroDivisionError:
            self.specificity = 0.0

        try: #positive predictive value
            if Interaction == "MutuallyExclusive":
                PPV = float(FN) / (FN + TN)
                self.PPV = PPV
            elif Interaction == "Co-occurring":
                PPV = float(TP) / (TP + FP)
                self.PPV = PPV
            else: 
                print "Developer error: 'Interaction' must be defined as either 'Co-occurring' or 'MutuallyExclusive'"
                print "Code died in Stastics:Performance, exiting"
                exit()
        except ZeroDivisionError:
            PPV = 0.0
            self.PPV = PPV

        try: #negative predictive value
            if Interaction == "MutuallyExclusive":
                NPV = float(FP) / (FP + TP)
                self.NPV = NPV
            elif Interaction == "Co-occurring":
                NPV = float(TN) / (TN + FN)
                self.NPV = NPV
            else: 
                print "Developer error: 'Interaction' must be defined as either 'Co-occurring' or 'MutuallyExclusive'"
                print "Code died in Stastics:Performance, exiting"
                exit()
        except ZeroDivisionError:
            NPV = 0.0
            self.NPV = NPV

        try: #accuracy 
            Accuracy = float(TP + TN) / (TP + FP + TN + FN)
            if Interaction == "MutuallyExclusive": self.accuracy = 1 - Accuracy
            elif Interaction == "Co-occurring": self.accuracy = Accuracy
            else: 
                print "Developer error: 'Interaction' must be defined as either 'Co-occurring' or 'MutuallyExclusive'"
                print "Code died in Stastics:Performance, exiting"
                exit()
        except ZeroDivisionError:
            self.accuracy = 0.0

        try: #mathews correlation coefficient
            MCC = ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
            if Interaction == "MutuallyExclusive": self.MCC = -1*MCC
            elif Interaction == "Co-occurring": self.MCC = MCC
            else: 
                print "Developer error: 'Interaction' must be defined as either 'Co-occurring' or 'MutuallyExclusive'"
                print "Code died in Stastics:Performance, exiting"
                exit()
        except ZeroDivisionError:
            self.MCC = 0.0

def p_adjust(PValues, Method="fdr"):
    '''
    Several multiple-testing correction procedures employed via the R Rpy2 interface. Method defs:
    BH = fdr = Benjamini and hochberg false discovery rate
    bonferroni
    holm
    hochberg
    hommel
    BY = Benjamini and Yekutieli
    qvalue = the method of Storey and Tibshirani, PNAS 2003
    none = bypass method
    '''

    Features = [Feature for Feature, PValue in sorted(PValues.items(), key=lambda Values: Values[1])]
    PValues = [PValue for Feature, PValue in sorted(PValues.items(), key=lambda Values: Values[1])]

    if Method.lower() == 'qvalue':
        qvalue = importr('qvalue')
        QValues = list(ro.r['qvalue'](ro.FloatVector(PValues))[2])
    else:
        QValues = list(ro.r['p.adjust'](ro.FloatVector(PValues), method=Method))
  
    return dict(zip(PValues, QValues))

class MOCASet:
    '''
    Boolean set operations for combining vectors of trues (1's) and falses (0's)
    '''

    def __init__(self, Vectors, NA="NA"):
        
        #The union
        if len(Vectors) > 1:
            self.union = [NA if NA in Elements and 1 not in Elements \
                              else (1 if 1 in Elements else 0) for Elements in zip(*Vectors)]
        elif len(Vectors) == 1:
            self.union = Vectors[0]
        else:
            self.union = Vectors

        #The intersection
        if len(Vectors) > 1:
            self.intersection = [NA if NA in Elements and 0 not in Elements \
                                 else (0 if 0 in Elements else 1) for Elements in zip(*Vectors)]
        elif len(Vectors) == 1:
            self.intersection = Vectors[0]
        else:
            self.intersection = Vectors
        
        #The difference 
        if len(Vectors) >= 1:
            self.difference = Vectors[0]
            if len(Vectors) > 1:
                for Vector in Vectors[1:len(Vectors)]:
                    self.difference = [0 if Element1 == 0 or Element2 == 1 \
                                           else (1 if Element1 == 1 and Element2 == 0 else NA) \
                                           for Element1, Element2 in zip(self.difference, Vector)]
        else:
            self.difference = Vectors

def weighted_sample(List, Length):
    '''
    Say you have a list that is enriched for some features (i.e., those features appear more than once
    in the list). The point of the enriched list is that there is a higher probability of sampling 
    list items that are more enriched. However, what if you wanted this enriched probability, but 
    didn't want to actually retrieve mutltiple, identical items (which is what might happen using 
    the python 'sample' function); that's the problem this function addresses. 
    '''

    if len(set(List)) <= Length: return list(set(List))

    Set = []
    while len(Set) < Length:
        Item = sample(List, 1)[0]
        if Item not in Set:
            Set.append(Item)

    return Set
        
def combine_features(Features, Min, Max):
    '''
    For a list of features return every possible combination of those features from "Min" to Max".
    So, if the list is [a,b,c,d], and you provide Min=1 and Max=3, this function would return:
    a, b, c, d, ab, ac, ad, bc, bd, cd, abc, abd, acd, bcd.
    '''

    Combinations = []
    [Combinations.extend(combos(Features, Feature)) for Feature in range(Min, Max + 1)]
   
    return sorted(Combinations)

def transpose_matrix(Labels, Features, Variates):
    '''
    Unfortunately numpy doesn't care for "NA" values, so I had to make my own. Features and labels
    just get switched.
    '''

    Matrix = [[] for Element in Variates[0]] #initialize the matrix
    for Vector in Variates:
        for Element in range(len(Vector)):
            Matrix[Element].append(Vector[Element])
    
    return Features, Labels, Matrix

def minimum_performance(performance, Arguments):
    '''
    You might wish to reject any marker with poor user-defined predictive 
    performance. Currently handles any or all of: sensitivity (can also use "sens"),
    specificity (can also use "spec"), positive predictive value (strictly "PPV"),
    and negative predictive value (strictly "NPV"). Easy to add anything from the 
    Performance function...will do later.
    '''
    
    if not Arguments.MinimumPerformance: return True

    Cutoffs = [(MinimumPerformance.split("=")[0][:3].lower(), float(MinimumPerformance.split("=")[1])) 
               for MinimumPerformance in Arguments.MinimumPerformance]

    Passed = True
    for Metric, Cutoff in Cutoffs:
        if Metric == "sen" and performance.sensitivity < Cutoff: Passed = False
        if Metric == "spe" and performance.specificity < Cutoff: Passed = False
        if Metric == "ppv" and performance.PPV < Cutoff: Passed = False
        if Metric == "npv" and performance.NPV < Cutoff: Passed = False
        if Metric == "acc" and performance.accuracy < Cutoff: Passed = False
    
    return Passed

def interaction(PValue):
    '''
    Given a p-value object determine if the corresponding interaction is mutually exclusive
    or co-occurring.
    '''

    if PValue.left_tail < PValue.right_tail:
        return "MutuallyExclusive"
    else:
        return "Co-occurring"

def rank(Scores, RankMethod, Top):
    '''
    How to rank (e.g., performance or significance) interactions. 
    'Scores' is a dictionary of feature (or some identifier) keys and values that are some sort of 
    score, such as p-value or predictive value. 'RankMethod' is likely passed by way of Arguments.RankMethod. 
    'Top' is just how many features you want returned.
    '''

    RankMethod = RankMethod.lower().replace("-", "")

    if RankMethod == "balancedaccuracy":
         SortedScores = sorted(Scores.items(), key=lambda performance: \
                                  (performance[1].sensitivity + performance[1].specificity)/2, reverse=True)[:Top]
         
    #There is no benefit in only ranking by sensitivity. Therefore we do sensitivity then specificity
    elif RankMethod == "sensitivity": 
         SortedScores = sorted(Scores.items(), key=lambda performance: \
                                  (performance[1].sensitivity, performance[1].specificity), reverse=True)[:Top]

    #There is no benefit in only ranking by specificity. Therefore we do specificity then sensitivity
    elif RankMethod == "specificity":
         SortedScores = sorted(Scores.items(), key=lambda performance: \
                                  (performance[1].specificity, performance[1].sensitivity), reverse=True)[:Top]

    elif RankMethod == "meanpredictivevalue":
         SortedScores = sorted(Scores.items(), key=lambda performance: \
                                  (performance[1].NPV + performance[1].PPV)/2, reverse=True)[:Top]

    elif RankMethod == "npv":
         SortedScores = sorted(Scores.items(), key=lambda performance: \
                                  (performance[1].NPV, performance[1].PPV), reverse=True)[:Top]

    elif RankMethod == "ppv":
         SortedScores = sorted(Scores.items(), key=lambda performance: \
                                  (performance[1].PPV, performance[1].NPV), reverse=True)[:Top]

    elif RankMethod == "accuracy":
         SortedScores = sorted(Scores.items(), key=lambda performance: \
                                  (performance[1].accuracy), reverse=True)[:Top]

    elif RankMethod == "mcc":
         SortedScores = sorted(Scores.items(), key=lambda performance: \
                                  (performance[1].MCC), reverse=True)[:Top]
         
    else:
        print "You chose a rank method that is not available. Please rerun using a valid rank method by setting"
        print "'RankMethod =' to one of the available methods (see the MOCA user's manual for all options)."
        exit()

    SortedKeys = zip(*SortedScores)[0]
    
    return SortedKeys

def confidence_interval(Probability, SampleSize, CI=0.95):
    '''
    For a given probability (probably sens, spec, PPV, or NPV) return the upper and 
    lower probabilities that define the confidence interval
    '''

    #Normal distribution table only included for CI = 90%, 95%, and 99%
    NormalDistribution = {0.90: 1.64, 0.95: 1.96, 0.99: 2.58}

    standardError = sqrt((Probability*(1 - Probability))/SampleSize)

    ConfidenceInterval = NormalDistribution[CI]*standardError
    
    return Probability - ConfidenceInterval, Probability + ConfidenceInterval
    
class EffectSize:
    '''
    Class that takes a=TP, b=FP, c=FN, d=TN from a 2x2 contingency and stores the effect size. 
    Two flavors are stored: the odds ratio and the difference of proportions. 
    '''

    def __init__(self, Interaction, a, b, c, d):

        a = float(a)
        b = float(b)
        c = float(c)
        d = float(d)

        #Because we can use the mutual exclusivity of a marker as the diagnostic, we might have flip the OR
        if Interaction == "Co-occurring":
            try:
                self.odds_ratio = (a/b)/(c/d)
            except ZeroDivisionError:
                self.odds_ratio = 0.0
                
        if Interaction == "MutuallyExclusive":
            try:
                self.odds_ratio = (c/d)/(a/b)
            except ZeroDivisionError:
                self.odds_ratio = 0.0

        try:
            n1 = a + b
            n2 = c + d
            p1 = a/n1
            p2 = c/n2
            self.difference_of_proportions = abs(p1 - p2)
        except ZeroDivisionError:
            self.difference_of_proportions = 0.0
        

