# -*- coding: utf-8 -*-
"""

@author: Rachel Schlossman

"""

from math import pi
from pyomo.environ import Constraint, Objective, Var, ConstraintList
from pyomo.environ import *
from pyomo.dae import *
import numpy as np
from pyomo.util.model_size import build_model_size_report
import matplotlib.pyplot as plt

from PlanarVehicle import PlanarVehicle
from CVX_obstacles import CVX_obstacles
from utilities.saveOptimizationVariables import saveOptimizationVariables
from utilities.loadOptimizationVariables import loadOptimizationVariables
from utilities.VarContainer import VarContainer
from utilities.plotter import plotResults

RAD_TO_DEG = 180 / pi
DEG_TO_RAD = pi / 180

def main():
    
    # Initial and final conditions
    x0=0
    y0=0
    theta0=0 
    v0=0
    omega0=0 
    xf=4
    yf=8
    thetaf=None
    vf=None
    omegaf=None
    
    for iter_num in range(2):
        
        # Create planar vehicle
        veh = PlanarVehicle(x0=x0, y0=y0, theta0=theta0, v0=v0, omega0=omega0, \
                            xf=xf, yf=yf, thetaf=thetaf,  vf=vf, omegaf=omegaf)
    
        # No obstacle in first iteration
        if iter_num == 0:
            CVX_obstacles(veh.m, obstacles_present=False)
        else:
            CVX_obstacles(veh.m, obstacles_present=True)
        
        # Boundary Conditions
        veh.m.bcCon         = ConstraintList(rule= veh.BCs)
        
        # Phase 1
        n1 = [1] # flag for Phase 1
        veh.m.ForceControlCon1     = Constraint(n1, veh.m.tau1, rule=veh.ForceControl)
        veh.m.TorqueControlCon1    = Constraint(n1, veh.m.tau1, rule=veh.TorqueControl)
        veh.m.ForceControlDot1     = Constraint(n1, veh.m.tau1, rule=veh.ForceControlDot)
        veh.m.TorqueControlDot1    = Constraint(n1, veh.m.tau1, rule=veh.TorqueControlDot)
        veh.m.xdotCon1             = Constraint(n1, veh.m.tau1, rule=veh.xdot)
        veh.m.ydotCon1             = Constraint(n1, veh.m.tau1, rule=veh.ydot)
        veh.m.thetadotCon1         = Constraint(n1, veh.m.tau1, rule=veh.thetadot)
        veh.m.vdotCon1             = Constraint(n1, veh.m.tau1, rule=veh.vdot)
        veh.m.omegadotCon1         = Constraint(n1, veh.m.tau1, rule=veh.omegadot)   
        veh.m.JdotCon1             = Constraint(n1, veh.m.tau1, rule=veh.JdotTotalTime)
        
        # Phase 2
        n2 = [2] # flag for Phase 2
        veh.m.ForceControlCon2     = Constraint(n2, veh.m.tau2, rule=veh.ForceControl)
        veh.m.TorqueControlCon2    = Constraint(n2, veh.m.tau2, rule=veh.TorqueControl)
        veh.m.ForceControlDot2     = Constraint(n2, veh.m.tau2, rule=veh.ForceControlDot)
        veh.m.TorqueControlDot2    = Constraint(n2, veh.m.tau2, rule=veh.TorqueControlDot)
        veh.m.xdotCon2             = Constraint(n2, veh.m.tau2, rule=veh.xdot)
        veh.m.ydotCon2             = Constraint(n2, veh.m.tau2, rule=veh.ydot)
        veh.m.thetadotCon2         = Constraint(n2, veh.m.tau2, rule=veh.thetadot)
        veh.m.vdotCon2             = Constraint(n2, veh.m.tau2, rule=veh.vdot)
        veh.m.omegadotCon2         = Constraint(n2, veh.m.tau2, rule=veh.omegadot)   
        veh.m.JdotCon2             = Constraint(n2, veh.m.tau2, rule=veh.JdotTotalTime)
        
        veh.m.obj  = Objective(expr = veh.m.J1[1] + veh.m.J2[1], sense=minimize)

        # Warm start initial guess from previous solution
        if iter_num > 0:
            loadOptimizationVariables(veh, myPyomoVars)
        
        # Solve optimization problem
        solver = SolverFactory('ipopt')
        solver.options["halt_on_ampl_error"] = "yes"
        solver.options["tol"] = 1e-6
        solver.options["max_iter"] = 2000
        solver.options["linear_solver"] = "ma27"
        solver.options["ma27_pivtol"] = 1e-6
        solver.options["linear_scaling_on_demand"] = "yes"
        
        results = solver.solve(veh.m, tee=True)

        # Save all data for next warm start
        _, myPyomoVars = saveOptimizationVariables(veh)

        
        trajectory_vars = VarContainer(veh.m)
        if iter_num == 1:
            plotResults(trajectory_vars, veh.m)
        
    plt.show()
    return veh.m
    
if __name__ == '__main__':
    m = main()