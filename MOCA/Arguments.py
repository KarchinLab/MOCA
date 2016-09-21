'''
Get python modules
'''
from argparse import ArgumentParser as argument_parser
from sys import argv
from ast import literal_eval
import re
import os.path

'''
Get third-party modules
'''

'''
Get MOCA modules
'''

def get_argument_file(Arguments):
    '''
    This function is used to read in an agrument file, instead of reading arguments
    from the command line. First, all defaults are set from the get_arguments function, 
    then overide defaults with any command that appears in the Arguments file. 
    '''

    class ArgumentObject:

        def __init__(self, Arguments):

            for Argument in Arguments:
                try: #first see if the argument is number
                    Value = float(Argument[1])
                except ValueError:
                    Value = Argument[1]
                if Argument[1].isdigit(): #specifically, is it an int?
                    vars(self)[Argument[0]] = int(Argument[1])
                elif Argument[1] == "True" or Argument[1] == "False": #should we convert to bool?
                    vars(self)[Argument[0]] = literal_eval(Value)
                elif type(Value) == float:
                    vars(self)[Argument[0]] = Value
                else: #do we need to convert to a python list?
                    #provided it isn't filename, which can have a space, but never be a list
                    if len(re.sub(r"[?]", " ",  Value).split()) > 1 and Argument[0] != "Filename":
                        vars(self)[Argument[0]] = re.sub(r"[?]", " ",  Value).split()
                    else:
                        vars(self)[Argument[0]] = Value

    #Make the argument-value pairs
    UserArguments = [(Argument.split()[0], " ".join(Argument.split()[2:]).replace(",", "")) \
                     for Argument in file("Arguments") \
                     if Argument.strip() and Argument.split()[1] == "="]
    
    #Now, overide any default arguments with those specified in the arguments file
    for Key, Value in ArgumentObject(UserArguments).__dict__.items():
        Arguments.__dict__[Key] = Value
    
    return Arguments

def get_arguments():
    '''
    Read in the MOCA defaults. Override defaults with from the command line or arguments file
    See UsersManual.pdf for more detail on commands
    '''

    parser = argument_parser()
    parser.add_argument("-a", "--pairwise", dest="Pairwise", action="store_true", default=False,
                        help="Run regular pairwise MOCA")
    parser.add_argument("-b", "--bandwidth", dest="Bandwidth", action="store_true", default=False,
                        help="Allow more than one of the same continuous-valued features, but with different cutoffs, into the same setwork. See MOCA user's manual for all options (default=False)")
    parser.add_argument("-c", "--datatype", dest="DataType", nargs=1, default=False,
                        help="For the purpose of preprocessing, tell MOCA if data if 'Binary', 'Continuous', or 'Categorical'")
    parser.add_argument("-d", "--data", dest="Data", nargs="*", 
                        help="Tell me which data types to run (e.g., Mutation, Drug, etc.).")
    parser.add_argument("-e", "--ejectpassengers", dest="EjectPassengers", default=0.1, type=float,
                        help="Remove features from setworks that aren't actually contributing to the significance of the setwork. See the MOCA user's manual for more details (default=0.1).")
    parser.add_argument("-f", "--filename", dest="Filename", default="Default", 
                        help="Input/output filname prefix for reading/writing files. If you don't specify, MOCA will use other parameters to make the filename")
    parser.add_argument("-g", "--featuremin", dest="FeatureMin", nargs="*", default=3, type=int, 
                        help="Minimum number of postives (i.e., values of '1') required in a vector for MOCA to consider it in any comparison. See MOCA user's manual for all options (default=3)")
    #h is reserved for help and cannot be used!
    parser.add_argument("-i", "--booleansets", dest="BooleanSets", nargs=3, default=[0, 0, 0], 
                        help="For a given setwork, this is the max number of features to be combined using the union, intersection, and difference Boolean set operations. See MOCA user's manual for all option")
    parser.add_argument("-j", "--minimumperformance", dest="MinimumPerformance", nargs="*", default=False, 
                        help="When writing output files, only include features that meet a certain performance criteria. For instance, spec=0.8 sens=0.8. See MOCA user's manual for all option (default=False).")
    parser.add_argument("-k", "--correctionmethod", dest="CorrectionMethod", default="BH",
                        help="Which method do you want to use to correct p-values? See MOCA user's manual for all options (default=Benjamini & Hochberg).")
    parser.add_argument("-l", "--leavesomeout", dest="LeaveSomeOut", default=0, type=int, \
                        help="Specificy the number to leave out for a leave-some-out cross-validation")
    parser.add_argument("-m", "--mode", dest="Mode", default="PreProcess", 
                        help="'PreProcess' formats data for subsequent calculations. 'PostProcess' for analyzing results. See MOCA user's manual for all options (default=PreProcess)")
    parser.add_argument("-n", "--normalize", dest="Normalize", action="store_true", default=False, 
                        help="Normalize the data before applying a 'Threshold'? As of now, the only type of normalization possible is z-score normalization")
    parser.add_argument("-o", "--optimization", dest="Optimization", nargs=3, default=[1000, 100, 0.01], 
                        help="Total cycles, frequency of feature-pool repopulation, and what fraction of top performers to repopulate with. See MOCA user's manual for more details (default=[1000, 100, 0.01].")
    parser.add_argument("-p", "--phenotype", dest="Phenotype",  action="store_true", default=False, 
                        help="Name of the phenotype (e.g., Drug, Cyst, Survival). Required for sets, optional for pairs (this would be the feature(s) you'd like to select markers for).")
    parser.add_argument("-q", "--fdr", dest="FDR", default=0.05, type=float, 
                        help="Minimum, adjusted p-value threshold -- anything larger is filtered. Default = 0.05.")
    parser.add_argument("-r", "--reports", dest="Reports", action="store_true", default=False, 
                        help="Tell me about the data in MOCA.data or MOCA.results. Options are 'Data' or 'Results'. Default=False")
    parser.add_argument("-s", "--seed", dest="Seed", default=0, type=int, 
                        help="Supply an interger seed (controls stochastic components of algorithm) -- if you don't, MOCA will. All seeds can be accessed using the Reports = Results argument")
    parser.add_argument("-t", "--threshold", dest="Threshold", nargs="*", default=False,
                        help="During 'PreProcessing' mode, tell MOCA how to threshold and dichotomize continuous-valued input. See MOCA user's manual for details (default=False)")
    parser.add_argument("-u", "--multiprocessmode", dest="MultiProcessMode", nargs=2, default=[0, 1], 
                        help="If distributing jobs accross multiple processors, the first option determines which node is running the job, and the second argument determines how many nodes are running")
    parser.add_argument("-v", "--na", dest="NA", default="NA", \
                        help="Tells MOCA what convention is used to delineate any missing values, such as 'nan' or 'NA' (default)")
    parser.add_argument("-w","--permutephenotype",dest="PermutePhenotype", action="store_true", default=False,
                        help="Permutation testing via shuffling the phenotype vector(s). Useful to see if, given your feature vectors, you could derive biomarkers for a permutation if your phenotype.")
    parser.add_argument("-x","--priors",dest="Priors", default=False,
                        help="Give MOCA a list of features to bias the search. Typically employed to focus the search if you know some features are important.  See users manual for details")
    parser.add_argument("-y","--labellist",dest="LabelList", default=False,
                        help="Provide MOCA with a list (one label per line), and the calculation will be restricted to that subset of labels")
    parser.add_argument("-z","--featurelist",dest="FeatureList", default=False,
                        help="Provide MOCA with a list (one feature per line), and the calculation will be restricted to that subset of features")
    parser.add_argument("--forcecooccurring",dest="ForceCooccurring", default=False,
                        help="If set to true, this will reject all MutuallyExclusive interactions -- some might find this easier to interpret")
    parser.add_argument("--combinedata",dest="CombineData", default=False,
                        help="")
    
    #Post-processing options get capital letters
    parser.add_argument("-H", "--heatmap", dest="Heatmap", nargs="*",
                        help="Supply arguments to alter the default behavior for making heatmaps in MOCA (see the user's manual for options).")
    parser.add_argument("-O", "--orderby", dest="OrderBy", nargs="*",
                        help="For various plotting functions, specify the ordering of features (eg., columns in a heatmap). Options are complex, see the user's manual for details.")
    parser.add_argument("-P","--postprocess",dest="PostProcess", default=False,
                        help="")
    parser.add_argument("-R", "--rankmethod", dest="RankMethod", default="BalancedAccuracy",
                        help="Tell MOCA how to sort and optimize features (default=BalancedAccuracy (arithmetic mean of sensitivity and specificity))")
    parser.add_argument("-T", "--topinteractions", dest="TopInteractions", nargs="*", default=100, type=int, 
                        help="Number of 'top' interactions (associations) to return during various post-processing tasks. 'Top' is defined by the score function you provide. (default=100)")
    parser.add_argument("-U", "--untransformeddata", dest="UntransformedData", nargs="*",
                        help="Useful for postprocessing and making plots, particular with data that was initially continuous valued. Tells MOCA which data to use in it's untransformed state.")
    
    #Developer options
    parser.add_argument("--mymoca", dest="MyMOCA", action="store_true", default=False, 
                        help="Build your own MOCA module! See MyMOCA.py")
    parser.add_argument("--profile", dest="Profile", action="store_true", default=False,
                        help="Profile to find inefficient code blocks")
    parser.add_argument("--debug", dest="Debug", action="store_true", default=False,
                        help="run in python debug mode 'pdb'")
    
    Arguments = parser.parse_args()
    
    #Overwrite defaults and command line arguments if the Arguments file exists
    if os.path.isfile("Arguments"): Arguments = get_argument_file(Arguments)
    
    return Arguments  
    
