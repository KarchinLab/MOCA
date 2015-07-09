'''
Get python modules
'''
from os import popen, listdir, getcwd

'''
Get third-party modules
'''

'''
Get MOCA modules
'''

def mkdir(Directory):
    '''
    Function to simulate Linux mkdir command. "-p" to avoid error if directory exists
    '''
    Command = "mkdir -p " + Directory
    popen(Command)
    
    return 
        
def mv(File, Destination):
    '''
    Function to simulate Linux mv command. 
    '''
    Command = "mv " + File + " " + Destination  
    popen(Command)
    
    return 
        
def rm(Thing):
    '''
    Function to simulate Linux rm command. "-rf" in case directory
    '''
    Command = "rm -rf " + Thing
    popen(Command)
    
    return 

def cp(Thing, Destination):
    '''
    Function to simulate Linux cp command. "-r" in case directory
    '''
    Command = "cp -r " + Thing + " " + Destination
    popen(Command)
    
    return 

def ls(Path):
    '''
    Function to simulate Linux ls command. 
    '''

    return listdir(Path)
    
def NoParens(Command):
    """
    Allows calling functions without parens (if no argument required).
    See "pwd", for instance. 
    """
    class Parens:
        def __init__(self, Command):
            self.Command = Command
        def __repr__(self):
            PythonCommand =apply(self.Command)
            if PythonCommand:
                return repr(PythonCommand)
            else:
                return ""
    return Parens(Command)


def pwd():
    """
    Get the current working directory
    """

    return getcwd()
pwd = NoParens(pwd)
