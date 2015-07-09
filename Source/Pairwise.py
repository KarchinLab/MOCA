'''
Get python modules
'''
from itertools import combinations
from collections import Counter
import cPickle
import os

'''
Get third-party modules
p'''
from fisher import pvalue as fisher
from numpy import arange, random

'''
Get MOCA modules
'''
from Statistics import contingency_table, p_adjust
from GetData import get_master_matrix, get_file, get_data_type


def pairwise(Arguments):
    '''
    Standard pairwise moca. Check every pairwise interactions for the 
    provided data type(s) to find the significant interactions 
    '''

    #Clustermode support if called
    Node = int(Arguments.MultiProcessMode[0]) 
    TotalNodes = int(Arguments.MultiProcessMode[1]) 

    Data = get_master_matrix(Arguments)
    Labels = Data["Labels"]
    Features = Data["Features"]
    Variates = Data["Variates"]

    #If you specify a phenotype, all features will just be compared with that phenotype...
    if Arguments.Phenotype:
        Features1 = [Feature for Feature in Features if Arguments.Phenotype in Feature]
        Features2 = [Feature for Feature in Features if Arguments.Phenotype not in Feature]
    #...Otherwise, every feature gets compared with every feature
    else:
        Features1 = Features
        Features2 = Features

    Results = []
    for Feature1 in Features1[Node::TotalNodes]:
        EmpiricalDistribution = {}
        for Feature2 in Features2:
            if Feature2 != Feature1:
                Response, Predictor = \
                    zip(*[Pair for Pair in zip(Variates[Features.index(Feature1)], Variates[Features.index(Feature2)])])
                a,b,c,d = contingency_table(Response, Predictor)
                EmpiricalDistribution[Feature2] = fisher(a,b,c,d).two_tail
       
        FDRs = p_adjust(EmpiricalDistribution, Arguments.CorrectionMethod)
        for Feature2 in Features2:
            if Feature2 != Feature1:
                FDR = FDRs[EmpiricalDistribution[Feature2]]
                if FDR < Arguments.FDR:
                    Results.append((Feature1, Feature2, FDR))
        
    Output = open(Arguments.Filename + ".pairs", "w")
    for Feature1, Feature2, FDR in Results:
        Output.write("%s %s %0.3e \n" %(Feature1, Feature2, FDR))
    Output.close()
    
    return
    
