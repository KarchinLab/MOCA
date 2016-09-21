'''
Get python modules
'''
import pdb
import random
import cProfile
import os

'''
Get third-party modules
'''
from numpy import random as numdom

'''
Get MOCA modules
'''
from MOCA.Arguments import get_arguments
from MOCA.DataHandler import pickle_matrices, get_path, combine_data
from MOCA.MyMOCA import my_moca
from MOCA.Pairwise import pairwise
from MOCA.Reports import reports
from MOCA.PostProcess.PostProcess import post_process
from MOCA.Setworks import setworks
from MOCA.LeaveSomeOut import leave_some_out

def main(Arguments):
    '''
    This routine exists just to call the appropriate 
    functions to run the desired MOCA mode.
    '''
    
    #Step thru the code line by line with python debug mode!
    if Arguments.Debug:
        pdb.set_trace()

    #If a seed isn't provided, set one and report to user.
    if Arguments.Seed:
        pass
    else:
        Arguments.Seed = int(1e5*random.random()) 
        
    random.seed( Arguments.Seed )
    numdom.seed( Arguments.Seed ) #we need to make control numpy's randomness as well!

    #For arguments that can take list, see if they were supplied one or fewer arguments and correct the format.
    if type(Arguments.Data) != list:
        Arguments.Data = [Arguments.Data]
    if type(Arguments.UntransformedData) != list:
        Arguments.UntransformedData = filter(lambda DataType: DataType != None, [Arguments.UntransformedData])
    if type(Arguments.OrderBy) != list:
        Arguments.OrderBy = filter(lambda Feature: Feature != None, [Arguments.OrderBy])
    if type(Arguments.Heatmap) != list:
        Arguments.Heatmap = filter(lambda Option: Option != None, [Arguments.Heatmap])
    if type(Arguments.Threshold) != list:
        Arguments.Threshold = filter(lambda Threshold: Threshold != None, [Arguments.Threshold])
        
    #Begin determining the mode we are in
    if Arguments.MyMOCA:
        my_moca(Arguments)

    elif Arguments.Reports:
        reports(Arguments)

    elif Arguments.Mode:
        if Arguments.Mode == "PreProcess":
            if Arguments.CombineData:
                combine_data(Arguments)
            else:
                pickle_matrices(Arguments)
            
        elif Arguments.Mode == "MOCA":
            if Arguments.Pairwise:
                pairwise(Arguments)
            elif Arguments.LeaveSomeOut:
                leave_some_out(Arguments)
            else:
                setworks(Arguments)
            
        elif Arguments.Mode == "PostProcess":
            post_process(Arguments)
        
        else:
            print "'PreProcess', 'MOCA', or 'PostProcess' are the only valid options for 'Mode'."
            print "Please fix your Arguments and run again."
            print "Exiting..."
            exit()
    else:
        print "To run MOCA you need to choose either a mode, request a report, or call your own function from MyMOCA."
        print "Please see the 'Arguments file' and/or execute moca.py --help. Exiting..."
        exit()

#Check to see that you have a paths file called 'Paths' in your current working directory
if not os.path.isfile("Paths"):
    print "You must have a paths file in your current-working directory, and it must be called 'Paths'"
    print "Exiting..."
    exit()

#MOCA data and results directories must exist for any MOCA calculation
try:
    os.makedirs(get_path("MOCA.data"))
except OSError:
    pass

try:
    os.makedirs(get_path("MOCA.results"))
except OSError:
    pass

# vvv Doesn't work well in MultiprocessMode vvv
#if not os.path.exists(get_path("MOCA.results")): os.makedirs(get_path("MOCA.results"))

#Test the speed of MOCA...
if get_arguments().Profile:
    cProfile.run( "main(get_arguments())", sort="cumulative" )
else: #...or just do a MOCA run
    main(get_arguments())
