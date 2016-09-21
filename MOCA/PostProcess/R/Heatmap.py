'''
Get python modules
'''
from itertools import chain
import os
from random import shuffle
from math import log10

'''
Get third-party modules
'''
from rpy2 import robjects as ro
from numpy import arange
import rpy2.rinterface as ri
ri.initr()

'''
Get MOCA modules
'''
from MOCA.DataHandler import load_data, get_ordered_matrix, get_path
from DataHandler import float_matrix, string_matrix, transpose_matrix
from Colors import Colors

'''
Load R libraries
'''
ro.r['source'](os.path.dirname(os.path.realpath(__file__)) + "/Functions.R")
#ro.r['library']('gplots')
#ro.r['library']('devtools')

def reorder_data(Data, Arguments):
    '''
    '''

    Features = list(chain(*Data.Features.values()))
    Variates = list(chain(*Data.Variates.values()))

    TransformedFeatures = list(chain(*Data.Transformed.Features.values()))
    TransformedVariates = list(chain(*Data.Transformed.Variates.values()))

    VariatesToSort = []
    for Feature in Arguments.OrderBy:
        #If the users requests it, categorical data can be sorted by the original strings, rather than the processed features
        if Feature in list(chain(*[Data.Features[DataType] for DataType in Arguments.CategoricalData])):
            VariatesToSort.append(Variates[Features.index(Feature)])
        else:
            try:
                VariatesToSort.append(TransformedVariates[TransformedFeatures.index(Feature)])
            except ValueError:
                VariatesToSort.append(map(float, Variates[Features.index(Feature)]))
            
    #'zip' our sort-by variates together, then zip them with the Labels and make a dictionary
    LabelsToSort = dict(zip(Data.Labels, zip(*VariatesToSort)))
    
    #Reorder the labels by sorting by the variates of every requested feature!
    Labels = []
    for Label, Variate in sorted(LabelsToSort.items(), key=lambda Variate: Variate[1][:]):
        Labels.append(Label)

    #Now we simply have to reorder ALL variates to match our new label order
    for DataType in Arguments.Data:
        Data.Variates[DataType] = get_ordered_matrix(Labels, Data.Labels, Data.Variates[DataType])
        Data.Transformed.Variates[DataType] = get_ordered_matrix(Labels, Data.Labels, Data.Transformed.Variates[DataType])

    Data.Labels = Labels

    return Data
    
def get_heatmap_data(Data, Arguments):
    '''
    This function gets the data to-be plotted on the heatmap itself, not the color sidebars.
    The user must specificy which data that is, using the 'UntransformedData = ' argument, and 
    that data has to be continuous. 
    '''
    
    if not Arguments.UntransformedData:
        print "You haven't specified which datatype(s) are to be plotted in the heatmap itself."
        print "Because heatmap data (i.e., everything aside from the colorbars) needs to be continuous,"
        print "you must tell MOCA which data that is using the 'UntransformedData =' argument. If any of"
        print "of this doesn't make sense, please see the MOCA user's manual."
        print "Exiting..."
        exit()

    for DataType in Arguments.UntransformedData:
        if Data.Report[DataType]["Type      "] != "Continuous":
            print DataType, "Is not continuous-valued. Please rerun, supplying only continuous-valued datatypes"
            print "to the 'UntransformedData =' argument. This is the data that gets plotted on the actual heatmap." 
            print "Exiting..."
            exit()

    Features = [Data.Features[DataType] for DataType in Arguments.Data \
                if DataType in Arguments.UntransformedData]
    Variates = [Data.Variates[DataType] for DataType in Arguments.Data \
                if DataType in Arguments.UntransformedData]
    
    return list(chain(*Features)), list(chain(*Variates))

def get_colsidematrix(Data, Arguments, colors=Colors().palette_one):
    '''
    '''

    if "randomcolors" in [Option.lower() for Option in Arguments.Heatmap]: shuffle(colors)

    #Binary features are always black ("1") and white ("0")
    ColorDict = {0:"white", 1:"black", Arguments.NA:"white"}

    Features = list(chain(*Data.Features.values()))
    Variates = list(chain(*Data.Variates.values()))

    TransformedFeatures = list(chain(*Data.Transformed.Features.values()))
    TransformedVariates = list(chain(*Data.Transformed.Variates.values()))

    ColorFeatures, ColorVariates, Legend = [], [], []
    for DataType in set(Arguments.Data) - set(Arguments.UntransformedData):
        if DataType in Arguments.CategoricalData:
            for Feature in Data.Features[DataType]:
                ColorVariates.append(Variates[Features.index(Feature)])
                Legend.append(Variates[Features.index(Feature)])
                ColorFeatures.append(Feature)
        else:
            for Feature in Data.Transformed.Features[DataType]:
                ColorVariates.append(TransformedVariates[TransformedFeatures.index(Feature)])
                ColorFeatures.append(Feature)

    #Uniquify the list of variates that get colors, zip it with some colors, and make a dictionary legend
    Legend = dict(zip(set(chain(*Legend)), colors[:len(set(chain(*Legend)))]))
    #Add the legend to the color dict so we can make the color matrix
    ColorDict = dict(ColorDict.items() + Legend.items())

    #Make the color matrix
    ColorMatrix = [[ColorDict[Variate] for Variate in ColorVariate] for ColorVariate in ColorVariates]
    ColorFeatures = [Feature[:Feature.index(":")] for Feature in ColorFeatures]
    
    return ColorFeatures, ColorMatrix, Legend

def get_rowsidematrix(Features, Arguments):
    '''
    '''
    
    GrayScale = {-9:"gray10", -8:"gray20", -7:"gray30", -6:"gray40", -5:"gray50", -4:"gray60", -3:"gray70", -2:"gray80"}

    Matrix = dict([(Row.split()[0], Row.split()[1]) for Row in file(get_path("RowSideColors"))])
    
    PValues = []
    for Feature in Features:
        PValues.append(Matrix[Feature])

    Label = ["Mesenchymal"]
    Variates = [int(("%0.2e" %float(PValue)).split("e")[1]) for PValue in PValues]
    Variates = ["black" if Variate <= -10 else ("white" if Variate >= -1 else GrayScale[Variate]) for Variate in Variates]

    return string_matrix(Features, Label, [Variates], Arguments)
   
def get_heatmap(
        Heatmap,
        ColSideMatrix=False,
        RowSideMatrix=False,
        ColorScheme="bluered",
        BreakBegin=-1.0,
        BreakEnd=1.0,
        scale="row",
        key=True,
        keysize=1.0,
        symbreaks=False,
        density="none",
        symkey=False, 
        trace="none", 
        cexRow=0.75, 
        cexCol=0.01,
        Rowv=True,
        Colv=True,
        BottomMargin=5,
        RightMargin=10,
        Legend=False):
    '''
    '''

    ro.r['source'](os.path.dirname(os.path.realpath(__file__)) + "/Heatmap3.R")
    ro.r['library']('gplots')
    ro.r['library']('devtools')

    Breaks = ro.FloatVector(list(arange(BreakBegin, BreakEnd+0.1, 0.001)))
    Cluster=ro.r('function(c) {hclust(c,method="average")}')
    Distance=ro.r('function(c) {dist(c,method="euclidean")}')

    if ColSideMatrix and not RowSideMatrix:
        ro.r['heatmap.3'](
            Heatmap,
            ColSideColors=ColSideMatrix,
            col=ro.r[ColorScheme](len(Breaks) -1),
            breaks=Breaks,
            hclustfun=Cluster,
            distfun=Distance,
            scale=scale,
            key=key,
            keysize=keysize,
            symbreaks=symbreaks,
            density=density,
            symkey=symkey, 
            trace=trace, 
            cexRow=cexRow, 
            cexCol=cexCol,
            Rowv=Rowv,
            Colv=False,
            margins=ro.IntVector([BottomMargin,RightMargin])
        )

    elif RowSideMatrix:
        ro.r['heatmap.3'](
            Heatmap,
            ColSideColors=ColSideMatrix,
            RowSideColors=RowSideMatrix,
            col=ro.r[ColorScheme](len(Breaks) -1),
            breaks=Breaks,
            hclustfun=Cluster,
            distfun=Distance,
            scale=scale,
            key=key,
            keysize=keysize,
            symbreaks=symbreaks,
            density=density,
            symkey=symkey, 
            trace=trace, 
            cexRow=cexRow, 
            cexCol=cexCol,
            Rowv=Rowv,
            Colv=False,
            margins=ro.IntVector([BottomMargin,RightMargin])
        )
        
    else:
        ro.r['heatmap.3'](
            Heatmap,
            col=ro.r[ColorScheme](len(Breaks) -1),
            breaks=Breaks,
            hclustfun=Cluster,
            distfun=Distance,
            scale=scale,
            key=key,
            keysize=keysize,
            symbreaks=symbreaks,
            density=density,
            symkey=symkey, 
            trace=trace, 
            cexRow=cexRow, 
            cexCol=cexCol,
            Rowv=Rowv,
            Colv=Colv,
            margins=ro.IntVector([BottomMargin,RightMargin])
        )

    if Legend:
        ncol = 1
        if len(Legend.values()) >= 6: ncol = 2
        Fill = ro.StrVector(Legend.values())
        Legend = ro.StrVector(Legend.keys())
        ro.r['legend'](
            "topright",
            legend=Legend,
            fill=Fill,
            border=False,
            bty="n",
            cex=0.6,
            ncol=ncol,
            **{'y.intersp':0.7}
        )
        
    return
    
def heatmap(Arguments):
    '''
    '''

    Data = load_data(Arguments)
    if Arguments.OrderBy: Data = reorder_data(Data, Arguments)

    #This is specifically the part you want plotted on the heatmap, NOT the sidebars
    Features, Variates = get_heatmap_data(Data, Arguments)
    Features = [Feature[:Feature.index(":")] for Feature in Features]

    ClusteredMatrix = ro.r["hierarchical_cluster_distance_matrix"](float_matrix(Data.Labels, Features, Variates, Arguments))
    Columns = list(ClusteredMatrix.colnames)
    Rows = list(ClusteredMatrix.rownames)
    Matrix = np.array(ClusteredMatrix)
    
    ZMatrix = z_matrix(Matrix)
    ColorPalette = list(ro.r["colorRampPalette"](ro.StrVector(["blue", "white", "red"]))(100))
    ColorMatrix = [float2color(Vector, ColorPalette) for Vector in ZMatrix]
    exit()
    
    Heatmap = float_matrix(Data.Labels, Features, Variates, Arguments)
    
    Legend = False #overwrite if ColSideColors is called

    if "colsidecolors" in [Option.lower() for Option in Arguments.Heatmap]:
        ColorFeatures, ColSideMatrix, Legend = get_colsidematrix(Data, Arguments)
        ColSideMatrix = transpose_matrix(string_matrix(Data.Labels, ColorFeatures, ColSideMatrix, Arguments))
    else:
        ColSideMatrix = False

    if "rowsidecolors" in [Option.lower() for Option in Arguments.Heatmap]:
        RowSideMatrix = get_rowsidematrix(Features, Arguments)
    else:
        RowSideMatrix = False

    ro.r['pdf'](Arguments.Filename + '.pdf')
    if Legend: ro.r['par'](mar=ro.IntVector([1,1,1,1]))
    ro.r['par'](**{'cex.axis':0.8})
    get_heatmap(Heatmap, ColSideMatrix=ColSideMatrix, RowSideMatrix=RowSideMatrix, Legend=Legend)
    ro.r['dev.off']()
    
    return
