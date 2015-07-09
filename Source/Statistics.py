'''
Get python modules
'''
from itertools import combinations as combos 
from itertools import permutations as permuts
from math import sqrt

'''
Get third-party modules
'''
from numpy import random, mean, median, std, arange, trapz, matrix, linalg
import rpy2.robjects as ro
from random import sample

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

#DEPRICATED, BUT CAN BE USED TO AVOID RPY2 USE
#def correlation(Vector1, Vector2):
#    '''
#    Return the correlation coefficient of two vectors
#    '''
#
#    return covariance(Vector1, Vector2)/(std(Vector1)*std(Vector2))

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

def contingency_table(Vector1, Vector2):
    '''
    Takes two vectors and makes a contingency table (aka
    confusion matrix). Useful for fishers test and calculating 
    performance metrics (i.e., specificity, sensitivity, etc.)
    '''

    a,b,c,d = 0,0,0,0
    for Element in range(len(Vector1)):
        if "NA" in [Vector1[Element], Vector2[Element]]: pass
        elif Vector1[Element] and Vector2[Element]: a += 1
        elif Vector1[Element]: b += 1
        elif Vector2[Element]: c += 1
        else: d += 1

    return a,b,c,d

class Performance:
    '''
    Class that computes standard performance metrics for statistical
    tests. Example:

    from Statistics import Performance

    performance = Performance(20, 3, 3, 4)
    print performance.specificity, performance.sensitivity 
    '''
    
    def __init__(self, TP=0, TN=0, FP=0, FN=0, P=0, N=0):

        try: #sensitivity (AKA true positive rate)
            self.sensitivity = float(TP) / (TP + FN) 
        except ZeroDivisionError:
            self.sensitivity = 0.0

        try: #specificity (AKA 1 - false positive rate)
            self.specificity = float(TN) / (FP + TN)
        except ZeroDivisionError:
            self.specificity = 0.0

        try: #negative predictive value
            self.NPV = float(TN) / (TN + FN)
        except ZeroDivisionError:
            self.NPV = 0.0
        
        try: #positive predictive value
            self.PPV = float(TP) / (TP + FP)
        except ZeroDivisionError:
            self.PPV = 0.0

        try: #accuracy 
            self.accuracy = float(TP + TN) / (P + N)
        except ZeroDivisionError:
            self.accuracy = 0.0

        try: #mathews correlation coefficient
            self.MCC = ((TP * TN) - (FP*FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
        except ZeroDivisionError:
            self.MCC = 0.0

def ROC(Predictions, Experiments, Correlated=True):
    '''
    Given a list of predictions (continuous variable) and corresponding binary experimenal 
    results (trues and falses), returns the AUC and a lot of ROC performance metrics as 
    an object. The performance metrics object allows you to determine the cutoff the optimizes
    a particular metric. For example, the following code would return the cutoff the 
    maximizes the sum of specificity and sensitivity. Set "Correlated=False" if decreasing 
    prediction scores correlate with experimental trues

    from Statistics import ROC
    from operator import attrgetter

    AUC, Metrix = ROC(Prediction, Experiment)
    print sorted(Metrix, key=lambda Metric: Metric.Specificity + Metric.Sensitivity, reverse=True)[0].CutOff
    
    '''
    class ROCDatum:
        def __init__(self, CutOff, Sensitivity, Specificity, NPV, PPV, Accuracy, MCC):
                self.CutOff = CutOff
                self.Sensitivity = Sensitivity
                self.Specificity = Specificity
                self.NPV = NPV
                self.PPV = PPV
                self.Accuracy = Accuracy
                self.MCC = MCC
        def __repr__(self):
                return repr((self.CutOff, self.Sensitivity, self.Specificity, \
                                 self.NPV, self.PPV, self.Accuracy, self.MCC))
    
    FPR, TPR, ROCData = [],[],[]
    for CutOff in sorted(set(Predictions), reverse=Correlated):
        if Correlated: #is data correlated (see function description)?
            BinaryPredictions = [1 if Prediction >= CutOff else 0 for Prediction in Predictions] 
        else: BinaryPredictions = [1 if Prediction <= CutOff else 0 for Prediction in Predictions] 
        TP, FN, FP, TN = contingency_table(Experiments, BinaryPredictions)
        performance = Performance(TP, TN, FP, FN, sum(Experiments), \
                                      len(Experiments) - sum(Experiments))
        FPR.append(1 - performance.specificity)
        TPR.append(performance.sensitivity)
        ROCData.append(ROCDatum(CutOff, \
                                    round(performance.sensitivity, 2), \
                                    round(performance.specificity, 2), round(performance.PPV, 2), \
                                    round(performance.NPV, 2), round(performance.accuracy, 2), \
                                    round(performance.MCC, 2)))
    
    return [trapz(TPR, FPR), ROCData] #AUC, Object of ROC data

def p_adjust(EmpiricalDistribution, Method="fdr"):
    '''
    Several multiple-testing correction procedures employed via the R Rpy2 interface. Method defs:
    BH = fdr = Benjamini and hochberg false discovery rate
    bonferroni
    holm
    hochberg
    hommel
    BY = Benjamini and Yekutieli
    none = bypass method
    '''

    Features = [Feature for Feature, PValue in sorted(EmpiricalDistribution.items(), key=lambda Values: Values[1])]
    PValues = [PValue for Feature, PValue in sorted(EmpiricalDistribution.items(), key=lambda Values: Values[1])]
    QValues = list(ro.r['p.adjust'](ro.FloatVector(PValues), method=Method))
  
    return dict(zip(PValues, QValues))

class MOCASet:
    '''
    Boolean set operations for combining vectors of trues (1's) and falses (0's)
    '''

    def __init__(self, Vectors):
        
        #The union
        if len(Vectors) > 1:
            self.union = ["NA" if "NA" in Elements else (1 if sum(Elements) >= 1 else 0) for Elements in zip(*Vectors)]
        elif len(Vectors) == 1:
            self.union = Vectors[0]
        else:
            self.union = Vectors

        #The intersection
        if len(Vectors) > 1:
            self.intersection = ["NA" if "NA" in Elements else (1 if sum(Elements) == len(Vectors) else 0) \
                                     for Elements in zip(*Vectors)]
        elif len(Vectors) == 1:
            self.intersection = Vectors[0]
        else:
            self.intersection = Vectors
        
        #The difference 
        if len(Vectors) >= 1:
            self.difference = Vectors[0]
            if len(Vectors) > 1:
                for Vector in Vectors[1:len(Vectors)]:
                    self.difference = ["NA" if "NA" in [Element1, Element2] else \
                                           (1 if Element1 - Element2 == 1 else 0) \
                                           for Element1, Element2 in zip(self.difference, Vector)]
        else:
            self.difference = Vectors

def weighted_sample(List, Length):
    '''
    Say you have a list that is enriched for some features (i.e., those features appear more than once
    in the list). The point of the enriched list is that there is a higher probability of sampling 
    list items that are more enriched. However, what if you wanted this enriched probability, but 
    didn't want to actually retrienve mutltiple, identical items (which is what might happen using 
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
   
    return Combinations

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

def minimum_performance(Arguments, Interaction, Response, Predictor):
    '''
    You might wish to reject any marker with poor user-defined predictive 
    performance. Currently handles any or all of: sensitivity (can also use "sens"),
    specificity (can also use "spec"), positive predictive value (strictly "PPV"),
    and negative predictive value (strictly "NPV"). Easy to add anything from the 
    Performance function...will do later.
    '''
    
    if not Arguments.MinimumPerformance: return True

    Cutoffs = [(MinimumPerformance.split("=")[0][:4].lower(), float(MinimumPerformance.split("=")[1])) 
               for MinimumPerformance in Arguments.MinimumPerformance]

    if Interaction == "MutuallyExclusive": Predictor = [0 if Element else 1 for Element in Predictor]

    TP,FN,FP,TN = contingency_table(Response, Predictor)
    performance = Performance(TP,TN,FP,FN)
    
    Passed = True
    for Metric, Cutoff in Cutoffs:
        if Metric == "sens" and performance.sensitivity < Cutoff: Passed = False
        if Metric == "spec" and performance.specificity < Cutoff: Passed = False
        if Metric == "ppv" and performance.PPV < Cutoff: Passed = False
        if Metric == "npv" and performance.NPV < Cutoff: Passed = False
    
    return Passed

    
    
        

