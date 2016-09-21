'''
Get python modules
'''
from itertools import chain

'''
Get third-party modules
'''
from numpy import arange
from rpy2 import robjects as ro

'''
Get MOCA modules
'''

def string_vector(Vector, Arguments):
    '''
    Convert a python vector to an R vector of strings 
    '''

    return ro.StrVector(["NaN" if Element == Arguments.NA else Element for Element in Vector])

def int_vector(Vector, Arguments):
    '''
    Convert a python vector to an R vector of integers
    '''

    return ro.IntVector(["NaN" if Element == Arguments.NA else Element for Element in Vector])

def float_vector(Vector, Arguments):
    '''
    Convert a python vector to an R vector of floating-point numbers
    '''

    return ro.FloatVector(["NaN" if Element == Arguments.NA else Element for Element in Vector])

def string_matrix(Labels, Features, Variates, Arguments):
    '''
    Takes the three standard MOCA matrix structures (labels, features, and variates) and converts them to 
    R string matrix, including the header and row names
    '''

    Variates = [["NaN" if Element == Arguments.NA else Element for Element in Variate] for Variate in Variates]
    Dimnames = ro.r['list'](string_vector(Labels, Arguments), string_vector(Features, Arguments))
    
    return ro.r['t'](ro.r['matrix'](string_vector(list(chain(*Variates)), Arguments), nrow=len(Labels), dimnames=Dimnames))

def int_matrix(Labels, Features, Variates, Arguments):
    '''
    Takes the three standard MOCA matrix structures (labels, features, and variates) and converts them to 
    R matrix of integers, including the header and row names
    '''

    Variates = [["NaN" if Element == Arguments.NA else Element for Element in Variate] for Variate in Variates]
    Dimnames = ro.r['list'](string_vector(Labels, Arguments), string_vector(Features, Arguments))
    
    return ro.r['t'](ro.r['matrix'](int_vector(list(chain(*Variates)), Arguments), nrow=len(Labels), dimnames=Dimnames))

def float_matrix(Labels, Features, Variates, Arguments):
    '''
    Takes the three standard MOCA matrix structures (labels, features, and variates) and converts them to 
    R matrix of float, including the header and row names
    '''

    Variates = [["NaN" if Element == Arguments.NA else Element for Element in Variate] for Variate in Variates]
    Dimnames = ro.r['list'](string_vector(Labels, Arguments), string_vector(Features, Arguments))

    return ro.r['t'](ro.r['matrix'](float_vector(list(chain(*Variates)), Arguments), nrow=len(Labels), dimnames=Dimnames))

def transpose_matrix(Matrix):
    '''
    Transposes a matrix that is alread in R-matrix format
    '''

    return ro.r['t'](Matrix)

