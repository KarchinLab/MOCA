'''
Get python modules
'''

'''
Get third-party modules
'''

'''
Get MOCA modules
'''
#from LeaveSomeOut import biomarkers, vote, make_priors
from Pairwise import pairwise
from Setworks import setworks
from LeaveSomeOut import leave_some_out
from R.Barplot import barplot
from R.Heatmap import heatmap
from DataHandler import write_matrix, make_priors
from Validate import validate_markers

def post_process(Arguments):
    ''' 
    This is simply a distibutor for all functions involved in 
    postprocessing MOCA results from previous calculations, which 
    is required to make sense of most MOCA output
    '''

    if Arguments.PostProcess:
        if Arguments.PostProcess.lower() == "barplot":
            barplot(Arguments)
        if Arguments.PostProcess.lower() == "heatmap":
            heatmap(Arguments)
        if Arguments.PostProcess.lower() == "writematrix":
            write_matrix(Arguments)
        if Arguments.PostProcess.lower() == "writesetworks":
            write_setworks(Arguments)
        if Arguments.PostProcess.lower() == "validatemarkers":
            validate_markers(Arguments)
        if Arguments.PostProcess.lower() == "makepriors":
            make_priors(Arguments)
    
    elif Arguments.Pairwise:
        pairwise(Arguments)
        
    elif Arguments.LeaveSomeOut:
        leave_some_out(Arguments)
    
    else:
        setworks(Arguments)

    return 

