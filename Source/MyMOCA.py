'''
Get python modules
'''

'''
Get third-party modulesz
'''

'''
Get MOCA modules
'''

##########################################################
#
# Just an example to illustrate the MyMOCA function.
# 
# This exists so that users can combine moca functions in 
# different ways, and extend their usage here
#
# Under the "my_moca" function, put in as many new functions 
# as you want to call, and add them to this file. 
#
# See the functions below as examples
#
# To call, set --mymoca fnc (from the command line) or 
# MyMOCA = fnc (from the Arguments file), where fnc is either
# "HelloWord" or "Goodbye" in the examples below
#
###########################################################

def hello_world(Arguments):
    '''
    Explain function here
    '''
   
    print "hello world"
   
    return 

def good_bye(Arguments):
    '''
    Explain function here
    '''
   
    print "good bye"
   
    return 

def my_moca(Arguments):
    
    if Arguments.MyMOCA == "HelloWorld":
        hello_world(Arguments)

    if Arguments.MyMOCA == "GoodBye":
        good_bye(Arguments)

    return
