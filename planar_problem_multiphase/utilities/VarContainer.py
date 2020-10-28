# -*- coding: utf-8 -*-
import numpy as np 

# Store the values of the optimal solution
class VarContainer():
    
    def __init__(self, m):

        # time
        self.tau1 = m.tau1
        self.tau2 = m.tau2

        self.time1 = np.dot(m.tau1,m.tf1())
        self.time2 = np.dot(m.tau2,m.tf2())+m.tf1()

        # states
        
        self.x1 = [m.x1[t]() * m.x_scale for t in m.tau1]
        self.y1 = [m.y1[t]() * m.y_scale for t in m.tau1]
        self.theta1 = [m.theta1[t]() * m.theta_scale for t in m.tau1]
        self.v1 = [m.v1[t]() * m.v_scale for t in m.tau1]
        self.omega1 = [m.omega1[t]() * m.omega_scale for t in m.tau1]
        self.F_control1 = [m.F_control1[t]() * m.F_scale for t in m.tau1]
        self.T_control1 = [m.T_control1[t]() * m.T_scale for t in m.tau1]        
        
        self.x2 = [m.x2[t]() * m.x_scale for t in m.tau2]
        self.y2 = [m.y2[t]() * m.y_scale for t in m.tau2]
        self.theta2 = [m.theta2[t]() * m.theta_scale for t in m.tau2]
        self.v2 = [m.v2[t]() * m.v_scale for t in m.tau2]
        self.omega2 = [m.omega2[t]() * m.omega_scale for t in m.tau2]
        self.F_control2 = [m.F_control2[t]() * m.F_scale for t in m.tau2]
        self.T_control2 = [m.T_control2[t]() * m.T_scale for t in m.tau2]              
        return
    
    
        