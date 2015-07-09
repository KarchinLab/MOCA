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
    This mode will run if there are either no commands, in which case looks for ./Arguments--
    if it finds this file, it will read it in. Alternatively, you can tell MOCA the path/name
    of the arguments file from the command line with the --argumentfile (-a) command.
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
                    if len(re.sub(r"[?]", " ",  Value).split()) > 1:
                        vars(self)[Argument[0]] = re.sub(r"[?]", " ",  Value).split()
                    else:
                        vars(self)[Argument[0]] = Value

    if len(argv) == 1:
        ArgumentFile = "Arguments.txt"
    else:
        ArgumentFile = argv[2]

    #Make the argument-value pairs
    UserArguments = [(Argument.split()[0], " ".join(Argument.split()[2:]).replace(",", "")) \
                     for Argument in file(Arguments.ArgumentFile) \
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
    parser.add_argument("-a", "--argumentfile", dest="ArgumentFile", default="Arguments", \
                            help="Use an argument file rather than the command line. Supply the name of the argument file. Note: MOCA automatically looks for './Arguments'.")
    parser.add_argument("-b", "--makebinarymatrix", dest="MakeBinaryMatrix", action="store_true", default=False, \
                            help="MOCA will look to binarize your continuous-valued data everytime. So, if you are doing many different runs with, say and expression matrix, binarize once with this command, and save yourself a lot of time on future runs")
    parser.add_argument("-c", "--seed", dest="Seed", default=0, type=int, \
                            help="Supply an interger seed when output reproducility is required (controls stochastic components of algorithm)")
    parser.add_argument("-d", "--data", dest="Data", nargs="*", \
                            help="Tell me which data types to run (e.g., Mutation, Drug, etc.); you must use the same name in your path file so that MOCA can find it.")
    parser.add_argument("-e", "--booleansets", dest="BooleanSets", nargs=3, default=[0, 0, 0], \
                            help="For a given setwork, this is the max number of features to be combined using the union, intersection, and difference Boolean set operations in a single comparison")
    parser.add_argument("-f", "--featuremin", dest="FeatureMin", nargs="*", default=3, type=int, \
                            help="Minimum number of postives (i.e., values of '1') required in a vector for MOCA to consider it in any comparison. Will speed up time if your data matrices contain many vectors too sparse to be meaningful")
    parser.add_argument("-g", "--continuous", dest="Continuous", action="store_true", default=False, \
                            help="Leave data as continuous, rather than binarizing? Has no meaning with setworks!")
    #h is reserved for help and cannot be used!
    parser.add_argument("-i", "--zdata", dest="ZData", nargs="*", \
                            help="MOCA assumes continuous-valued input data has to be normalized. If you're providing MOCA with already normalized data, use this option. If any data is already normalized, please list a bool for each data type (e.g., --data Drug Expression --ZData True False).")
    parser.add_argument("-j", "--filename", dest="Filename", default="Default", \
                            help="Input/output filname prefix for reading/writing files. By default, filename prefix's are just all requested datatypes joined (e.g., DrugTissue.Pair or MutationExpression.sif")
    parser.add_argument("-k", "--phenotype", dest="Phenotype",  action="store_true", default=False, \
                            help="Name of the phenotype (e.g., Drug, Cyst, Survival) Required for sets, optional for pairs (this would be the feature(s) you'd like to select markers for)")
    parser.add_argument("-l", "--leavesomeout", dest="LeaveSomeOut", default=0, type=int, \
                            help="Specificy the number to leave out for a leave-some-out cross-validation")
    parser.add_argument("-m", "--mymoca", dest="MyMOCA", action="store_true", default=False, \
                            help="Build your own MOCA module! See MyMOCA.py")
    parser.add_argument("-n", "--postprocess", dest="PostProcess", action="store_true", default=False, \
                            help="Choose a method for postproccessing your MOCA results")   
    parser.add_argument("-o", "--optimization", dest="Optimization", nargs=3, default=[1000, 100, 0.01], \
                            help="The first param specifies the number of optimization cycles for making MOCA setworks, the second is the frequency with which to repopulate the feature pool with top-predicting features, and the third is the percent of top features repopulate with. So the default says to do 1000 total cycles, and every 100 cycles put the top 1 percent of features into the feature pool.")
    parser.add_argument("-p", "--pathfile", dest="PathFile", default="Paths", \
                            help="Show me where the path file is. Default is './Paths'")
    parser.add_argument("-q", "--fdr", dest="FDR", default=0.05, type=float, \
                            help="Minimum, adjusted p-value threshold")
    parser.add_argument("-r", "--minimumperformance", dest="MinimumPerformance", nargs="*", default=False, \
                            help="When writing output files, only include features that meet a certain performance criteria. For instance, spec=0.8 sens=0.8. Paricularly useful for cross-validation voting")
    parser.add_argument("-s", "--setworks", dest="Setworks", action="store_true", default=False, \
                            help="Make the MOCA multi-feature, Boolean-set-operation networks (i.e., setworks)")
    parser.add_argument("-t","--tissuespecific",dest="TissueSpecific", default=False, \
                            help="Use this option to tell moca to consider only one of the tissue types you have represented in your header (labels). Requires colon notation (e.g., Liver:Sample1)")
    parser.add_argument("-u", "--multiprocessmode", dest="MultiProcessMode", nargs=2, default=[0, 1], \
                            help="If distributing jobs accross multiple processors, the first option determines which node is running the job, and the second argument determines how many nodes are running")
    parser.add_argument("-v", "--correctionmethod", dest="CorrectionMethod", default="BH", \
                            help="Which method do you want to use to correct p-values? Options are: BH=Benjamini & Hochberg; bonferroni, holm, hochberg, hommel, BY = Benjamini and Yekutieli; none = bypass method")
    parser.add_argument("-w","--priors",dest="Priors", default=False, \
                            help="Use this option to tell moca to use priors for making setworks. You will have to point MOCA to a Priors file, using the Paths file, that will tell MOCA which features to use as priors for any or all of the Boolean set operations. If not False, Priors must be set to one of three options: Strict, Weighted, or Stochastic. See users manual for details")
    parser.add_argument("-x", "--pairwise", dest="Pairwise", action="store_true", default=False, \
                            help="Run regular pairwise MOCA")
    
    parser.add_argument("-z", "--zmin", dest="ZMin", nargs="*", default=0.8, type=float, \
                            help="Z-score cutoff to convert continuous variates to discrete")
    
    Arguments = parser.parse_args()
    
    #Overwrite  any defaults from the Arguments file
    if os.path.isfile(Arguments.ArgumentFile): Arguments = get_argument_file(Arguments)
    
    return Arguments  
    
