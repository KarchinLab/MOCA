'''
Get python modules
'''
import random

'''
Get third-party modules
'''
from numpy import random as numdom

'''
Get MOCA modules
'''
from Source.Pairwise import pairwise
from Source.GetArguments import get_arguments
from Source.MyMOCA import my_moca
from Source.GetData import make_binary_matrix
from Source.PostProcess.PostProcess import post_process
from Source.Setworks import setworks
from Source.LeaveSomeOut import leave_some_out

def main(Arguments):
    '''
    This routine exists just to call the appropriate 
    functions to run the desired MOCA mode.
    '''

    if Arguments.Data == None:
        print "You didn't specify any data!!!"
        print "Exiting..."
        exit()
    if type(Arguments.Data) != list:
        Arguments.Data = [Arguments.Data]
    if Arguments.Filename == "Default":
        Arguments.Filename = "".join(Arguments.Data)

    if Arguments.Seed:
        random.seed( Arguments.Seed )
        numdom.seed( Arguments.Seed )

    if Arguments.MyMOCA:
        my_moca(Arguments)

    elif Arguments.Setworks:
        setworks(Arguments)

    elif Arguments.LeaveSomeOut:
        leave_some_out(Arguments)

    elif Arguments.Pairwise:
        pairwise(Arguments)

    elif Arguments.MakeBinaryMatrix:
        make_binary_matrix(Arguments)

    elif Arguments.PostProcess:
        post_process(Arguments)

    else:
        print "You ran MOCA without turning on a mode! Please see the 'Arguments file' and/or execute moca.py --help"
        print "exiting..."

main(get_arguments())
