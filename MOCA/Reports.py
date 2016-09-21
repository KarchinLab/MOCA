'''
Get python modules
'''
import os
import cPickle

'''
Get third-party modules
'''

'''
Get MOCA modules
'''
from DataHandler import get_path

def data(Arguments):
    '''
    Get info for the files in your MOCA.data directory.
    '''
    
    Path = get_path("MOCA.data") + "/"
    Files =  os.listdir(Path)
    
    if not len(Files):
        print "You set 'Reports = Data', yet your MOCA.data directory is empty"
        print "Try processing some data first. Exiting..."
        
    for File in Files:
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print " "
        print "Processed data file =", Path + File
        print " "
        Report = cPickle.load(open(Path + File, "rb"))["Report"]
        Entries = Report["Entries"]
        for Entry in Entries:
            if type(Report[Entry]) == list:
                print "\t", Entry, "\t", " ".join(map(str, Report[Entry]))
            else:
                print "\t", Entry, "\t", Report[Entry]
                
        print " "
        print " "

    return

def results(Arguments):
    '''
    Get info for the files in your MOCA.results directory
    '''

    Path = get_path("MOCA.results") + "/"
    Files =  os.listdir(Path)
    
    if not len(Files):
        print "You set 'Reports = Results', yet your MOCA.results directory is empty"
        print "Try runnning some MOCA calculations first. Exiting..."
        exit()
        
    for File in Files:
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print " "
        print "MOCA results file =", Path + File
        print " "
        Report = cPickle.load(open(Path + File, "rb"))["Report"]
        Entries = Report["Entries"]
        for Entry in Entries:
            if type(Report[Entry]) == list:
                print "\t", Entry, "\t", " ".join(map(str, Report[Entry]))
            else:
                print "\t", Entry, "\t", Report[Entry]
                
        print " "
        print " "

    return

def reports(Arguments):
    '''
    The function prints out the contents of your MOCA.data and MOCA.results directories. 
    The print-out is accompanied by useful, file-specific info, which can be used to determine
    which processed data you want to do some MOCA calculations on, or which previously completed
    MOCA calculations you want to do some post-processing on. 
    '''

    try:
        if Arguments.Reports.lower() == "data":
            data(Arguments)

        elif Arguments.Reports.lower() == "results":
            results(Arguments)
        
    except AttributeError: 
        print "Only valid arguments for 'Reports' is either 'False', 'Data', or 'Results'"
        print "Please fix your arguments file. Exiting..."
        exit()
        
    return 
