#       __PLOTTINGTOOLS__
#       Plotting functions
#
#       Authors: 
#       Kostas Alexis (konstantinos.alexis@mavt.ethz.ch)

import matplotlib.pyplot as mpplot
from matplotlib import pyplot
import pylab
from mpl_toolkits.mplot3d import Axes3D
    
def plot3(a,b,c,mark="o",col="r"):
    # mimic matlab plot3
    from matplotlib import pyplot
    import pylab
    pylab.ion()
    fig = pylab.figure()
    ax = Axes3D(fig)
    ax.invert_zaxis()
    ax.invert_xaxis()
    ax.set_aspect('equal', 'datalim')
    ax.plot(a, b,c,color=col,marker=mark)
    fig.show()