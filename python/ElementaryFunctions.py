#       __DUBINSAIRPLANEFUNCTIONS__
#       A small set of helping functions
#
#       Authors: 
#       Kostas Alexis (konstantinos.alexis@mavt.ethz.ch)

import __builtin__
    
def max(a, b):
    # just a simple max function between 2 values
    if a >=b:
        return a
    else:
        return b

def min(a, b=0, nargout=0):
    # just a simple min function between 2 values
    if a <= b:
        return a
    else:
        return b
