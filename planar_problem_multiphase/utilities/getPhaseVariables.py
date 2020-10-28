# -*- coding: utf-8 -*-
"""

@author: Rachel Schlossman

"""

def getPhaseVariables(m, n, t):
    
    if n == 1:
        # time
        tf = m.tf1

        # state variables
        x = m.x1[t] * m.x_scale
        y = m.y1[t] * m.y_scale
        theta = m.theta1[t] * m.theta_scale
        v = m.v1[t] * m.v_scale
        omega = m.omega1[t] * m.omega_scale
        J = m.J1[t]
        
        # derivatives of state variables with respect to tau1
        dx_dtau = m.dx_dtau1[t] * m.x_scale 
        dy_dtau = m.dy_dtau1[t] * m.y_scale
        dtheta_dtau = m.dtheta_dtau1[t] * m.theta_scale
        dv_dtau = m.dv_dtau1[t] * m.v_scale
        domega_dtau = m.domega_dtau1[t] * m.omega_scale
        dJ_dtau = m.dJ_dtau1[t]
        
        # Controllable inputs 
        F_control = m.F_control1[t] * m.F_scale
        T_control = m.T_control1[t] * m.T_scale
        
        F_controldot = m.F_controldot1[t] * m.F_scale
        T_controldot = m.T_controldot1[t] * m.T_scale
        
        F_virtual = m.F_virtual1[t] * m.F_scale
        T_virtual = m.T_virtual1[t] * m.T_scale
        
        # derivatives of control variables with respect to tau
        dF_control_dtau = m.dF_control_dtau1[t] * m.F_scale
        dT_control_dtau = m.dT_control_dtau1[t] * m.T_scale
        
    elif n == 2:
        # time
        tf = m.tf2

        # state variables
        x = m.x2[t] * m.x_scale
        y = m.y2[t] * m.y_scale
        theta = m.theta2[t] * m.theta_scale
        v = m.v2[t] * m.v_scale
        omega = m.omega2[t] * m.omega_scale
        J = m.J2[t]
        
        # derivatives of state variables with respect to tau2
        dx_dtau = m.dx_dtau2[t] * m.x_scale 
        dy_dtau = m.dy_dtau2[t] * m.y_scale
        dtheta_dtau = m.dtheta_dtau2[t] * m.theta_scale
        dv_dtau = m.dv_dtau2[t] * m.v_scale
        domega_dtau = m.domega_dtau2[t] * m.omega_scale
        dJ_dtau = m.dJ_dtau2[t]
        
        # Controllable inputs 
        F_control = m.F_control2[t] * m.F_scale
        T_control = m.T_control2[t] * m.T_scale
        
        F_controldot = m.F_controldot2[t] * m.F_scale
        T_controldot = m.T_controldot2[t] * m.T_scale
        
        F_virtual = m.F_virtual2[t] * m.F_scale
        T_virtual = m.T_virtual2[t] * m.T_scale
        
        # derivatives of control variables with respect to tau
        dF_control_dtau = m.dF_control_dtau2[t] * m.F_scale
        dT_control_dtau = m.dT_control_dtau2[t] * m.T_scale
    
    return (tf, x, y, theta, v, omega, J, \
            dx_dtau, dy_dtau, dtheta_dtau, dv_dtau, domega_dtau, dJ_dtau, \
            F_control, T_control, F_controldot, T_controldot, \
            F_virtual, T_virtual, dF_control_dtau, dT_control_dtau)
        
    