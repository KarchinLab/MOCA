'''
Get python modules
'''

'''
Get third-party modules
'''

'''
Get MOCA modules
'''
from Setworks import biomarkers, vote, make_priors
from Pairwise import make_filtered_network

def post_process(Arguments):
    ''' 
    This is simply a distibutor for all functions involved in 
    postprocessing MOCA results from previous calculations, which 
    is required to make sense of most MOCA output
    '''

    if Arguments.PostProcess == "MakeFilteredNetwork":
        make_filtered_network(Arguments)

    if Arguments.PostProcess == "ValidateBiomarkers":
        biomarkers(Arguments)

    if Arguments.PostProcess == "ValidateVote":
        vote(Arguments)

    if Arguments.PostProcess == "MakePriors":
        make_priors(Arguments)

    return 

