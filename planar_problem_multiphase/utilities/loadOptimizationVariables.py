# -*- coding: utf-8 -*-
from pyomo.environ import *
from pyomo.dae import *
import numpy as np

# For warm start
def loadOptimizationVariables(vehicle, myPyomoVars):

    for myVar in vehicle.m.component_objects(Var, active=True):
        print("Initializing Variable: ", myVar)
        if len(myVar) == 1:
            myVar.set_value(myPyomoVars[myVar.name]())
        else:
            myVar.set_values(dict(zip(list(myVar.index_set()), np.interp(myVar.index_set(), list(myPyomoVars[myVar.index_set().getname()]), myPyomoVars[myVar.name][:]()))))

    return
