# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

import numpy as np

def plotResults(VarContainer, m):
    
    blue = '#1f77b4ff'
    orange = '#ff7f0eff'
    green = '#2ca02cff'
    
        
    plt.figure(1)
    plt.plot(VarContainer.x1, VarContainer.y1, '.-', color=blue)
    plt.plot(VarContainer.x2, VarContainer.y2, '.-', color='r')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid()
    
    ells = [Ellipse((m.ell1CenterX.value, m.ell1CenterY.value), width=m.ell1Width.value, height=m.ell1Height.value,
                    angle=0)]


    ax = plt.subplot(1, 1, 1)
    for e in ells:
        ax.add_artist(e)
        e.set_facecolor(np.dot(1 / 255, [255, 204, 204]))  # 255 204 204
        e.set_alpha(0.5)
    
    plt.figure(2)
    plt.plot(VarContainer.time1, VarContainer.theta1, '.-', color=blue)
    plt.plot(VarContainer.time2, VarContainer.theta2, '.-', color='r')
    plt.xlabel('time (s)')
    plt.ylabel('theta (rad)')
    plt.grid()
    
    plt.figure(3)
    plt.plot(VarContainer.time1, VarContainer.v1, '.-', color=blue)
    plt.plot(VarContainer.time2, VarContainer.v2, '.-', color='r')
    plt.xlabel('time (s)')
    plt.ylabel('v (m/s)')
    plt.grid()
 
    plt.figure(4)
    plt.plot(VarContainer.time1, VarContainer.omega1, '.-', color=blue)
    plt.plot(VarContainer.time2, VarContainer.omega2, '.-', color='r')
    plt.xlabel('time (s)')
    plt.ylabel('omega (rad/s)')
    plt.grid()    

    plt.figure(5)
    plt.plot(VarContainer.time1, VarContainer.F_control1, '.-', color=blue)
    plt.plot(VarContainer.time2, VarContainer.F_control2, '.-', color='r')
    plt.xlabel('time (s)')
    plt.ylabel('F control (N)')
    plt.grid()    
    
    plt.figure(6)
    plt.plot(VarContainer.time1, VarContainer.T_control1, '.-', color=blue)
    plt.plot(VarContainer.time2, VarContainer.T_control2, '.-', color='r')
    plt.xlabel('time (s)')
    plt.ylabel('T control (N-m)')
    plt.grid()    
    
    
    return
