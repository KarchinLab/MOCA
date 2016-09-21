'''
Get python modules
'''
from collections import Counter

'''
Get third-party modules
'''
from rpy2 import robjects as ro

'''
Get MOCA modules
'''
from MOCA.PostProcess.Setworks import load_setworks
from DataHandler import int_matrix

def group_continuous_features(Features):
    '''
    '''

    GreaterThan = {}
    for Feature in Features:
        if ">=" in Feature:
            Type, Threshold = Feature.replace(".>=", " ").split()
            GreaterThan.setdefault(Type, []).append(Threshold)

    LessThan = {}
    for Feature in Features:
        if "<" in Feature:
            Type, Threshold = Feature.replace(".<", " ").split()
            LessThan.setdefault(Type, []).append(Threshold)

    Features = [Feature for Feature in Features if ">=" not in Feature and "<" not in Feature]
    
    for k, v in GreaterThan.items():
        Features.extend([k + ".>=" + sorted(v)[0] for i in v])

    for k, v in LessThan.items():
        Features.extend([k + ".<" + sorted(v, reverse=True)[0] for i in v])
    
    return Features

def get_counts(Features):
    '''
    '''

    Features = group_continuous_features(Features)

    #
    Features = [Feature[:Feature.index(":")] for Feature in Features]
    
    SortedFeatures, Counts = [], []
    for Feature, Count in sorted(Counter(Features).items(), key=lambda x: x[1], reverse=True):
        SortedFeatures.append(Feature)
        Counts.append(Count)

    return SortedFeatures, Counts
        
def get_barplot(FeatureCounts, Main, Arguments):
    '''
    '''

    Features, Counts = FeatureCounts

    if not len(Features):
        ro.r["plot"](1, type='n', xlab='', ylab='', main=Main, axes=False)

    else:
        ro.r['barplot'](
            int_matrix(Features, ["row1"], [Counts], Arguments),
            las=3,
            main=Main,
            **{'cex.names':0.75}
        )
    
    return

def barplot(Arguments):
    '''
    '''

    Barcodes, PValues, QValues, Performances, Interactions, \
        UnionFeatures, IntersectionFeatures, DifferenceFeatures, \
        FeatureVectors, SampleCount, Report = load_setworks(Arguments)

    ro.r['pdf'](Arguments.Filename + '.pdf')
    ro.r['par'](mar=ro.IntVector([8,3,1,1]), mfrow=ro.IntVector([3,3]))
    
    CoUnion, CoInt, CoDiff, MutUnion, MutInt, MutDiff = [],[],[],[],[],[]
    for Barcode in Barcodes:
        if Interactions[Barcode] == "Co-occurring":
            CoUnion.extend(UnionFeatures[Barcode])
            CoInt.extend(IntersectionFeatures[Barcode])
            CoDiff.extend(DifferenceFeatures[Barcode])
        elif Interactions[Barcode] == "MutuallyExclusive":
            MutUnion.extend(UnionFeatures[Barcode])
            MutInt.extend(IntersectionFeatures[Barcode])
            MutDiff.extend(DifferenceFeatures[Barcode])

    get_barplot(get_counts(CoUnion), "Union:Co-occurring", Arguments)
    get_barplot(get_counts(CoInt), "Intersection:Co-occurring", Arguments)
    get_barplot(get_counts(CoDiff), "Difference:Co-occurring", Arguments)
    get_barplot(get_counts(MutUnion), "Union:MutuallyExclusive", Arguments)
    get_barplot(get_counts(MutInt), "Intersection:MutuallyExclusive", Arguments)
    get_barplot(get_counts(MutDiff), "Difference:MutuallyExclusive", Arguments)
    
    ro.r['dev.off']()
    
    return
