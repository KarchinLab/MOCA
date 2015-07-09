'''
Get python modules
'''

'''
Get third-party modules
'''
from fisher import pvalue as fisher

'''
Get MOCA modules
'''
from Source.Statistics import p_adjust, contingency_table, Performance, minimum_performance
from Source.GetData import get_master_matrix

def make_filtered_network(Arguments):
    '''
    This will make you a .sif file for all significant pairwise interactions determined 
    when you ran pairwise moce and generated the .pairs file. The resulting .sif file
    will be much more informaticve and can be directly fed into CytoScape
    '''

    Data = get_master_matrix(Arguments)
    Features = Data["Features"]
    Variates = Data["Variates"]

    EmpiricalDistribution = {}
    
    try: #find the .pairs file you should've created in pairwise moca
        Interactions = [Interaction.split() for Interaction in file(Arguments.Filename + ".pairs")]
        for Feature1, Feature2, FDR in Interactions: EmpiricalDistribution[tuple(sorted([Feature1, Feature2]))] = float(FDR)
    except IOError:
        print "Looking for", Arguments.Filename + ".pairs in current directory. If you have not yet created the pairs",
        print "please do so before making filtered networks. If you have made it, please place it in the current working",
        print "directory, and give it the correct name or use the 'Filename' argument to tell me which file to open."
        print "Exiting..."
        exit()
        
    Output = open(Arguments.Filename + ".sif", "w")
    for Feature1, Feature2 in EmpiricalDistribution:
        FDR = EmpiricalDistribution[tuple(sorted([Feature1, Feature2]))]
        if Arguments.Data[0] in Feature1:
            X = Variates[Features.index(Feature1)]
            Y = Variates[Features.index(Feature2)]
        else:
            X = Variates[Features.index(Feature2)]
            Y = Variates[Features.index(Feature1)]
        Response, Predictor = zip(*[Pair for Pair in zip(X, Y)])
        TP,FN,FP,TN = contingency_table(Response, Predictor)
        PValue = fisher(TP,FN,FP,TN)
        if PValue.left_tail < PValue.right_tail:
            Interaction = "MutuallyExclusive"
        else:
            Interaction = "Co-occurring"

        #You might wish to shrink the network size to only those interactions that met you predictive 
        #performance criteria. Probably only relevant when you ran pairwise moca with a Phenotype
        #specified (i.e., you were doing biomarker selection). 
        if minimum_performance(Arguments, Interaction, Response, Predictor):
            performance = Performance(TP,TN,FP,FN) 
            if Interaction == "MutuallyExclusive":
                Output.write("%s %s %s %0.2e %0.2e %0.2f %0.2f \n" %(Feature1, Interaction, Feature2, PValue.two_tail, FDR, 1.0 - performance.sensitivity, 1.0 - performance.specificity))
            if Interaction == "Co-occurring":
                Output.write("%s %s %s %0.2e %0.2e %0.2f %0.2f \n" %(Feature1, Interaction, Feature2, PValue.two_tail, FDR, performance.sensitivity, performance.specificity))
    Output.close()
    
    return
