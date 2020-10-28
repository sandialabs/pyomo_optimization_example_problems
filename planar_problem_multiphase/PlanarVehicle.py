# -*- coding: utf-8 -*-
"""

@author: Rachel Schlossman

"""

from math import pi
from pyomo.environ import *
from pyomo.dae import *

from utilities.getPhaseVariables import getPhaseVariables
from AeroBehavior import AeroBehavior

class PlanarVehicle(AeroBehavior):
    
    def __init__(self, 
                 x0=0, y0=0, theta0=None, v0=None, omega0=None, 
                 xf=None, yf=None, thetaf=None, vf=None, omegaf=None):
        
        super().__init__()
        
        # model
        self.m      = ConcreteModel()
        
        # scaling
        self.m.x_scale = Param(initialize = 10)
        self.m.y_scale = Param(initialize = 20)
        self.m.theta_scale = Param(initialize = 20)
        self.m.v_scale = Param(initialize = 5)
        self.m.omega_scale = Param(initialize = 20)
        self.m.F_scale = Param(initialize = 5)
        self.m.T_scale = Param(initialize = 5)
        
        # mass
        self.m.mass = Param(initialize = 10)
        # moment of inertial
        self.m.I = Param(initialize = 40)
        
        #=============Phase 1===========================================
        
        # time
        self.m.tau1  = ContinuousSet(bounds=(0,1))
        self.m.tf1 = Var(initialize = 30, bounds=(0.01, 60)) 

        
        # state variables
        self.m.x1     = Var(self.m.tau1)   
        self.m.y1     = Var(self.m.tau1)
        self.m.theta1 = Var(self.m.tau1)
        self.m.v1     = Var(self.m.tau1, bounds=(0,None))
        self.m.omega1 = Var(self.m.tau1)
        self.m.J1     = Var(self.m.tau1)
        
        # derivatives of state variables with respect to tau1
        self.m.dx_dtau1     = DerivativeVar(self.m.x1, wrt=self.m.tau1)
        self.m.dy_dtau1     = DerivativeVar(self.m.y1, wrt=self.m.tau1)
        self.m.dtheta_dtau1 = DerivativeVar(self.m.theta1, wrt=self.m.tau1)
        self.m.dv_dtau1     = DerivativeVar(self.m.v1, wrt=self.m.tau1)
        self.m.domega_dtau1 = DerivativeVar(self.m.omega1, wrt=self.m.tau1)
        self.m.dJ_dtau1     = DerivativeVar(self.m.J1, wrt=self.m.tau1)
        
        # Controllable inputs 
        self.m.F_control1 = Var(self.m.tau1, bounds=(-5/self.m.F_scale,5/self.m.F_scale)) 
        self.m.T_control1 = Var(self.m.tau1, bounds=(-5/self.m.T_scale,5/self.m.T_scale))
        
        self.m.F_controldot1 = Var(self.m.tau1, bounds=(-50/self.m.F_scale,50/self.m.F_scale)) 
        self.m.T_controldot1 = Var(self.m.tau1, bounds=(-50/self.m.F_scale,50/self.m.F_scale)) 
        
        self.m.F_virtual1 = Var(self.m.tau1)
        self.m.T_virtual1 = Var(self.m.tau1)
        
        # derivatives of control variables with respect to tau
        self.m.dF_control_dtau1     = DerivativeVar(self.m.F_control1, wrt=self.m.tau1)
        self.m.dT_control_dtau1     = DerivativeVar(self.m.T_control1, wrt=self.m.tau1)
        
        #================Phase 2==============================================
        
        # time
        self.m.tau2  = ContinuousSet(bounds=(0,1))
        self.m.tf2 = Var(initialize = 30, bounds=(0.01, 60)) 

        
        # state variables
        self.m.x2     = Var(self.m.tau2)   
        self.m.y2     = Var(self.m.tau2)
        self.m.theta2 = Var(self.m.tau2)
        self.m.v2     = Var(self.m.tau2, bounds=(0,None))
        self.m.omega2 = Var(self.m.tau2)
        self.m.J2     = Var(self.m.tau2)
        
        # derivatives of state variables with respect to tau2
        self.m.dx_dtau2     = DerivativeVar(self.m.x2, wrt=self.m.tau2)
        self.m.dy_dtau2     = DerivativeVar(self.m.y2, wrt=self.m.tau2)
        self.m.dtheta_dtau2 = DerivativeVar(self.m.theta2, wrt=self.m.tau2)
        self.m.dv_dtau2     = DerivativeVar(self.m.v2, wrt=self.m.tau2)
        self.m.domega_dtau2 = DerivativeVar(self.m.omega2, wrt=self.m.tau2)
        self.m.dJ_dtau2     = DerivativeVar(self.m.J2, wrt=self.m.tau2)
        
        # Controllable inputs 
        self.m.F_control2 = Var(self.m.tau2, bounds=(-5/self.m.F_scale,5/self.m.F_scale)) 
        self.m.T_control2 = Var(self.m.tau2, bounds=(-5/self.m.T_scale,5/self.m.T_scale))
        
        self.m.F_controldot2 = Var(self.m.tau2, bounds=(-50/self.m.F_scale,50/self.m.F_scale)) 
        self.m.T_controldot2 = Var(self.m.tau2, bounds=(-50/self.m.F_scale,50/self.m.F_scale)) 
        
        self.m.F_virtual2 = Var(self.m.tau2)
        self.m.T_virtual2 = Var(self.m.tau2)
        
        # derivatives of control variables with respect to tau
        self.m.dF_control_dtau2     = DerivativeVar(self.m.F_control2, wrt=self.m.tau2)
        self.m.dT_control_dtau2     = DerivativeVar(self.m.T_control2, wrt=self.m.tau2)    
    
        # discretize problem
        discretizer = TransformationFactory('dae.finite_difference')
        discretizer.apply_to(self.m, nfe=200, wrt=self.m.tau1, scheme='BACKWARD')
        discretizer = TransformationFactory('dae.finite_difference')
        discretizer.apply_to(self.m, nfe=200, wrt=self.m.tau2, scheme='BACKWARD')
        
        # initial and final conditions
        self.x0=x0
        self.y0=y0
        self.theta0=theta0
        self.v0=v0
        self.omega0=omega0
        
        self.xf=xf
        self.yf=yf
        self.thetaf=thetaf
        self.vf=vf
        self.omegaf=omegaf
        
        return
    
    # Dynamics    
    def xdot(self, m, n,t):

        (tf, x, y, theta, v, omega, J, dx_dtau, dy_dtau, \
         dtheta_dtau, dv_dtau, domega_dtau, dJ_dtau, \
         F_control, T_control, F_controldot, T_controldot, \
         F_virtual, T_virtual, dF_control_dtau, dT_control_dtau) = getPhaseVariables(m, n, t)
        
        return (dx_dtau) == tf * (v) * cos(theta)
    
    def ydot(self, m, n, t):

        (tf, x, y, theta, v, omega, J, dx_dtau, dy_dtau, \
         dtheta_dtau, dv_dtau, domega_dtau, dJ_dtau, \
         F_control, T_control, F_controldot, T_controldot, \
         F_virtual, T_virtual, dF_control_dtau, dT_control_dtau) = getPhaseVariables(m, n, t)
        
        return (dy_dtau) == tf * (v) * sin(theta)
    
    def thetadot(self, m, n, t):

        (tf, x, y, theta, v, omega, J, dx_dtau, dy_dtau, \
         dtheta_dtau, dv_dtau, domega_dtau, dJ_dtau, \
         F_control, T_control, F_controldot, T_controldot, \
         F_virtual, T_virtual, dF_control_dtau, dT_control_dtau) = getPhaseVariables(m, n, t)
        
        return (dtheta_dtau) == tf * (omega)
    
    def ForceControl(self, m, n, t):
        
        (tf, x, y, theta, v, omega, J, dx_dtau, dy_dtau, \
         dtheta_dtau, dv_dtau, domega_dtau, dJ_dtau, \
         F_control, T_control, F_controldot, T_controldot, \
         F_virtual, T_virtual, dF_control_dtau, dT_control_dtau) = getPhaseVariables(m, n, t)
        
        (cv, cvw) = self.getCvCvw(v)    
        
        F_aero = 0.1 * (omega) + cv * (v)
        
        return F_control == -F_aero + F_virtual
        
    def TorqueControl(self, m, n, t):
        
        (tf, x, y, theta, v, omega, J, dx_dtau, dy_dtau, \
         dtheta_dtau, dv_dtau, domega_dtau, dJ_dtau, \
         F_control, T_control, F_controldot, T_controldot, \
         F_virtual, T_virtual, dF_control_dtau, dT_control_dtau) = getPhaseVariables(m, n, t)
        
        (cv, cvw) = self.getCvCvw(v)    
        
        T_aero = 0.2 * (theta)  + cvw * (v)
        
        return T_control == -T_aero + T_virtual
    
    def ForceControlDot(self, m, n, t):
        
        (tf, x, y, theta, v, omega, J, dx_dtau, dy_dtau, \
         dtheta_dtau, dv_dtau, domega_dtau, dJ_dtau, \
         F_control, T_control, F_controldot, T_controldot, \
         F_virtual, T_virtual, dF_control_dtau, dT_control_dtau) = getPhaseVariables(m, n, t)
        
        return dF_control_dtau == tf * F_controldot
    
    def TorqueControlDot(self, m, n, t):
        
        (tf, x, y, theta, v, omega, J, dx_dtau, dy_dtau, \
         dtheta_dtau, dv_dtau, domega_dtau, dJ_dtau, \
         F_control, T_control, F_controldot, T_controldot, \
         F_virtual, T_virtual, dF_control_dtau, dT_control_dtau) = getPhaseVariables(m, n, t)        
        
        return dT_control_dtau == tf * T_controldot
    
    def vdot(self, m, n, t):   
        
        (tf, x, y, theta, v, omega, J, dx_dtau, dy_dtau, \
         dtheta_dtau, dv_dtau, domega_dtau, dJ_dtau, \
         F_control, T_control, F_controldot, T_controldot, \
         F_virtual, T_virtual, dF_control_dtau, dT_control_dtau) = getPhaseVariables(m, n, t)      
        
        return m.mass * (dv_dtau) == tf * (F_virtual)
                
    def omegadot(self, m, n, t):        

        (tf, x, y, theta, v, omega, J, dx_dtau, dy_dtau, \
         dtheta_dtau, dv_dtau, domega_dtau, dJ_dtau, \
         F_control, T_control, F_controldot, T_controldot, \
         F_virtual, T_virtual, dF_control_dtau, dT_control_dtau) = getPhaseVariables(m, n, t)     
        
        return m.I * (domega_dtau) == tf * (T_virtual) 
    
    def JdotTravelTime(self, m, n, t):
        
        (tf, x, y, theta, v, omega, J, dx_dtau, dy_dtau, \
         dtheta_dtau, dv_dtau, domega_dtau, dJ_dtau, \
         F_control, T_control, F_controldot, T_controldot, \
         F_virtual, T_virtual, dF_control_dtau, dT_control_dtau) = getPhaseVariables(m, n, t)    
        
        ## Penalize total travel in (x,y):
        return dJ_dtau == tf*((x-m.x1[0])**2 + (y-m.y1[0])**2)

    def JdotTotalTime(self, m, n, t):
        
        (tf, x, y, theta, v, omega, J, dx_dtau, dy_dtau, \
         dtheta_dtau, dv_dtau, domega_dtau, dJ_dtau, \
         F_control, T_control, F_controldot, T_controldot, \
         F_virtual, T_virtual, dF_control_dtau, dT_control_dtau) = getPhaseVariables(m, n, t)    
        
        ## Penalize total time:
        return dJ_dtau == tf*1
    
    # Boundary conditions
    def BCs(self, m):
        
        # initial x, y position
        yield m.x1[0]     == self.x0 / m.x_scale
        yield m.y1[0]     == self.y0  / m.y_scale
        if self.theta0 != None:   
            yield m.theta1[0] == self.theta0  / m.theta_scale
        if self.v0 != None:
            yield m.v1[0]     == self.v0 / m.v_scale
        if self.omega0 != None:
            yield m.omega1[0] == self.omega0 / m.omega_scale
        
        # Boundary conditions 1->2
        yield m.x1[1]     == m.x2[0]
        yield m.y1[1]     == m.y2[0] 
        yield m.theta1[1] == m.theta2[0]
        yield m.v1[1]     == m.v2[0]
        yield m.omega1[1] == m.omega2[0]  
        
        # mid-trajectory constraint
        yield m.x1[1]     == 3.5 / m.x_scale
        yield m.y1[1]     == 2.5 / m.y_scale
        
        # Final conditions
        if self.xf != None:
            yield m.x2[1] == self.xf / m.x_scale
        if self.yf != None:
            yield m.y2[1] == self.yf / m.y_scale
        if self.thetaf != None:            
            yield m.theta2[1] == self.thetaf  / m.theta_scale 
        if self.vf != None:
            yield m.v2[1] == self.vf / m.v_scale
        if self.omegaf != None:
            yield m.omega2[1] == self.omegaf / m.omega_scale
                
        # Objective function initialize
        yield m.J1[0] == 0
        yield m.J2[0] == 0
        


    
    