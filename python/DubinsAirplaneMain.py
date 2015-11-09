#       __DUBINSAIRPLANEMAIN__
#       This is the main file to execute examples of the Dubins Airplane mode
#       that supports 16 cases of possible trajectories
#
#       Authors: 
#       Kostas Alexis (konstantinos.alexis@mavt.ethz.ch)

from DubinsAirplaneFunctions import * 
from PlottingTools import plot3
import numpy as np
import time
import sys

pi = np.pi
dubins_case = 0
verbose_flag = 0
plot_flag = 1

class ExecutionFlags(object):
    """
        Execution flags
    """  
    def __init__(self, verbose_flag, plot_flag):
        self.verbose = verbose_flag
        self.plot = plot_flag

class VehicleParameters(object):
    """
        Vehicle Parameters
    """ 
    def __init__(self, Vairspeed_0, Bank_max, Gamma_max):
        self.Vairspeed_0 = Vairspeed_0
        self.Bank_max = Bank_max
        self.Gamma_max = Gamma_max
    
def main():
    # Example main for the 16 cases Dubins Airplane paths
    t0 = time.clock()
    VehiclePars = VehicleParameters( 15, pi/4, pi/6 ) 
    ExFlags = ExecutionFlags( verbose_flag,plot_flag )
    
    flag_nc = 0
    fname = 'path_dubins_solution.txt'
    
    if dubins_case == 1: # short climb RSR
        print '### Path Type: short climb RSR' 
        start_node  = np.array( [0,   0,   -100,   0*pi/180,    VehiclePars.Vairspeed_0] )
        end_node    = np.array( [0, 200,   -125, 270*pi/180,    VehiclePars.Vairspeed_0] )
    elif dubins_case == 2: # short climb RSL
        print '### Path Type: short climb RSL'
        start_node  = np.array([ 0,   0,    -100, -70*pi/180,    VehiclePars.Vairspeed_0] )
        end_node    = np.array( [100, 100,  -125, -70*pi/180,    VehiclePars.Vairspeed_0] )
    elif dubins_case == 3: # short climb LSR
        print '### Path Type: short climb LSR'
        start_node  = np.array( [0,   0,    -100, 70*pi/180,    VehiclePars.Vairspeed_0] )
        end_node    = np.array( [100, -100, -125, 70*pi/180,    VehiclePars.Vairspeed_0] )
    elif dubins_case == 4: # short climb LSL
        print '### Path Type: short climb LSL'
        start_node  = np.array( [0,   0,    -100,  70*pi/180,   VehiclePars.Vairspeed_0] )
        end_node    = np.array( [100, -100, -125, -135*pi/180,  VehiclePars.Vairspeed_0] )
    elif dubins_case == 5: # long climb RSR
        print '### Path Type: long climb RSR'
        start_node  = np.array( [0,   0,   -100,   0*pi/180,    VehiclePars.Vairspeed_0] )
        end_node    = np.array( [0, 200,   -250, 270*pi/180,    VehiclePars.Vairspeed_0] )
    elif dubins_case == 6: # long climb RSL
        print '### Path Type: long climb RSL'
        start_node  = np.array( [0,   0,    -100, -70*pi/180,    VehiclePars.Vairspeed_0] )
        end_node    = np.array( [100, 100,  -350, -70*pi/180,    VehiclePars.Vairspeed_0] )
    elif dubins_case == 7: # long climb LSR
        print '### Path Type: long climb LSR'
        start_node  = np.array( [0,   0,    -350, 70*pi/180,    VehiclePars.Vairspeed_0] )
        end_node    = np.array( [100, -100, -100, 70*pi/180,    VehiclePars.Vairspeed_0] )
    elif dubins_case == 8: # long climb LSL
        print '### Path Type: long climb LSL'
        start_node  = np.array( [0,   0,    -350,  70*pi/180,   VehiclePars.Vairspeed_0] )
        end_node    = np.array( [100, -100, -100, -135*pi/180,  VehiclePars.Vairspeed_0] )
    elif dubins_case == 9: # intermediate climb RLSR (climb at beginning)
        print '### Path Type: intermediate climb RLSR (climb at beginning)'
        start_node  = np.array( [0,   0,   -100,   0*pi/180,     VehiclePars.Vairspeed_0] )
        end_node    = np.array( [0, 200,   -200, 270*pi/180,    VehiclePars.Vairspeed_0] )
    elif dubins_case == 10: # intermediate climb RLSL (climb at beginning)
        print '### Path Type: intermediate climb RLSL (climb at beginning)'
        start_node  = np.array( [0,   0,   -100, 0*pi/180,     VehiclePars.Vairspeed_0] )
        end_node    = np.array( [100, 100, -200, -90*pi/180,   VehiclePars.Vairspeed_0] )
    elif dubins_case == 11: # intermediate climb LRSR (climb at beginning)
        print '### Path Type: intermediate climb LRSR (climb at beginning)'
        start_node  = np.array( [0,   0,   -100, 0*pi/180,     VehiclePars.Vairspeed_0] )
        end_node    = np.array( [100, -100, -200, 90*pi/180,    VehiclePars.Vairspeed_0] )
    elif dubins_case == 12: # intermediate climb LRSL (climb at beginning)
        print '### Path Type: intermediate climb LRSL (climb at beginning)'
        start_node  = np.array( [0,   0,   -100, 0*pi/180,     VehiclePars.Vairspeed_0] )
        end_node    = np.array( [100, -100, -200, -90*pi/180,   VehiclePars.Vairspeed_0] )
    elif dubins_case == 13: # intermediate climb RSLR (descend at end)
        print '### Path Type: intermediate climb RSLR (descend at end)'
        start_node  = np.array( [0,   0,   -200, 0*pi/180,     VehiclePars.Vairspeed_0] )
        end_node    = np.array( [100, 100, -100, 90*pi/180,    VehiclePars.Vairspeed_0] )
    elif dubins_case == 14: # intermediate climb RSRL (descend at end)
        print '### Path Type: intermediate climb RSRL (descend at end)'
        start_node  = np.array( [0,   0,   -200, 0*pi/180,     VehiclePars.Vairspeed_0] )
        end_node    = np.array( [100, 100, -100, -90*pi/180,   VehiclePars.Vairspeed_0] )
    elif dubins_case == 15: # intermediate climb LSLR (descend at end)
        print '### Path Type: intermediate climb LSLR (descend at end)'
        start_node  = np.array( [0,   0,   -200, 70*pi/180,     VehiclePars.Vairspeed_0] )
        end_node    = np.array( [100, -100, -100, 90*pi/180,    VehiclePars.Vairspeed_0] )
    elif dubins_case == 16: # intermediate climb LSRL (descend at end)
        print '### Path Type: intermediate climb LSRL (descend at end)'
        start_node  = np.array( [0,   0,   -150, 0*pi/180,     VehiclePars.Vairspeed_0] )
        end_node    = np.array( [100, -100,-100, -90*pi/180,   VehiclePars.Vairspeed_0] )
    elif dubins_case == 0: # for fixing errors
        print '### Path Type: for fixing errors'
        start_node  = np.array( [0,   0,   0, 0, VehiclePars.Vairspeed_0] )
        end_node    = np.array( [40, -140,  100, 11*pi/9, VehiclePars.Vairspeed_0] ) # LSRL
        #end_node    = np.array( [40, -140,  140, 2*pi/9, VehiclePars.Vairspeed_0] ) # LSLR

        #end_node    = np.array( [40, 140,  140, 11*pi/9, VehiclePars.Vairspeed_0] ) # RSLR
        #end_node    = np.array( [40, 140,  140, 1*pi/9, VehiclePars.Vairspeed_0] ) # RSRL

        #end_node    = np.array( [40, 140,  -140, 11*pi/9, VehiclePars.Vairspeed_0] ) # RLSR
        end_node    = np.array( [60, 140,  -140, 0*pi/14, VehiclePars.Vairspeed_0] ) # RLSL
        #end_node    = np.array( [40, -140,  -100, 11*pi/9, VehiclePars.Vairspeed_0] ) # LRSL        
        end_node    = np.array( [40, -140,  -100, 10*pi/180, VehiclePars.Vairspeed_0] ) # LRSR        
        
    if dubins_case > 16:
        flag_nc = 1
        print 'Not a case'
    
    if flag_nc == 0:
        print '### Case loaded.'
        print '### computing path...'
        R_min = MinTurnRadius_DubinsAirplane( VehiclePars.Vairspeed_0, VehiclePars.Bank_max )
        
        # Check if start and end node are too close. Since spiral-spiral-spiral (or curve-curve-curve) paths are not considered, the optimal path may not be found... (see Shkel, Lumelsky, 2001, Classification of the Dubins set, Prop. 5/6). Below is a conservative bound, which seems (by experiments) to assure a unproblematical computation of the dubins path.
        if ( np.linalg.norm(end_node[0:2] - start_node[0:2],ord=2) < 6*R_min ):
            print "!!!!!!!!!!!!!!!!"
            print "Conservative condition (end_node[0:2] - start_node[0:2],ord=2) < 6*R_min) not fulfilled!"
            print "Start and end pose are close together. Path of type RLR, LRL may be optimal"
            print "May fail to compute optimal path! Aborting"
            print "!!!!!!!!!!!!!!!!"
            sys.exit()
        
        
        DubinsAirplaneSolution = DubinsAirplanePath( start_node, end_node, R_min, VehiclePars.Gamma_max )
        if ExFlags.verbose :
            PrintSolutionAgainstMATLAB( DubinsAirplaneSolution )
        
        path_dubins_airplane = ExtractDubinsAirplanePath( DubinsAirplaneSolution )
        path_dubins_airplane = path_dubins_airplane.T
        print '### Execution time = ' + str( time.clock() - t0 )
        np.savetxt( fname, path_dubins_airplane.T, delimiter = ',' ) 
        print '### Dubins airplane solution saved in ' + fname
        if ExFlags.plot :
            print '### Dubins airplane solution plot'
            plot3( path_dubins_airplane[:,0], path_dubins_airplane[:,1], path_dubins_airplane[:,2], 'o', 'g' )
        
        # print 'Press any button to continue'
        # raw_input()


def testAllCases():
    #for dubins_case in xrange(1, 16):
    main()
    print 'Press any button to continue'
    raw_input()
        

if __name__ == "__main__":
    testAllCases()
