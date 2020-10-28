# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 10:05:20 2019

@author: kwilli2

Taken from page 12-24 in "Feedback Systems, An Introduction for Scientists and Engineers" 
by Karl Johan Astrom and Richard M. Murray

This is a modified version of the problem above, in which the controller has separate gains
for the proportional and integral control: u = Kp + Ki/s

"""
from math import pi
from pyomo.environ import *
from pyomo.dae import *
import control
import numpy as np
import sympy
import numpy as np
import sympy
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rc
import matplotlib as mpl

class GainScheduler():
    
    def __init__(self):

        self.m      = ConcreteModel()
        self.m.tau  = ContinuousSet(bounds=(0,1))

        return
    
    #Boundary conditions
    def BCs(self, m):
        yield m.y1[0] == 0
        yield m.y2[0] == 0 
        yield m.y3[0] == 0
        yield m.u[0] == 0
        yield m.r[0] == 8
        yield m.J[0] == 0
        
        
    def _y1dot(self, m, t):    
        return m.dy1dt[t] == m.tf * m.y2[t]
    
    def _y2dot(self, m, t):    
        return m.dy2dt[t] == m.tf * m.y3[t]
    
    def _y3dot(self, m, t):    
        return m.dy3dt[t] == m.tf * (-8*m.y1[t] -8*m.y2[t] -4*m.y3[t] + m.u[t])
    
    def _udot(self, m, t):    
        return m.dudt[t] == m.tf * ( m.Kp*(m.drdt[t]/m.tf-m.y2[t]) + m.Ki*(m.r[t]-m.y1[t]) )
    
    def _rdot(self, m, t):    
        # return m.drdt[t] == m.tf * (0)  # r is a constant 
        # return m.drdt[t] == m.tf * (1)  # r is a slope
        return m.drdt[t] == m.tf * (cos(10*t))  # r is a sine wave
        
        
    def _ForwardPathTF(self):
        w = sympy.Symbol('w', real=True)
        s = sympy.I*w
        Kp = sympy.Symbol('Kp', real=True)
        Ki = sympy.Symbol('Ki', real=True)
        C = Kp + Ki/s
        G = 1/((s+2)*(s**2+2*s+4))
        L = C*G
        
        num,den = L.as_numer_denom()
        temp = (den*den.conjugate()).as_real_imag()[0]
        real,imag = (num*den.conjugate()/temp).as_real_imag()    
        _real = sympy.lambdify((Kp,Ki,w), real)
        _imag = sympy.lambdify((Kp,Ki,w), imag)
        return (_real, _imag)
        
    def _freqPhase180Cross0dBCross(self, m):    
        w = sympy.Symbol('w', real=True)
        s = sympy.I*w
        C = m.Kp() + m.Ki()/s
        G = 1/((s+2)*(s**2+2*s+4))
        L = C*G    
        
        num,den = L.as_numer_denom()
        temp = (den*den.conjugate()).as_real_imag()[0]
        a,b = (num*den.conjugate()/temp).as_real_imag()    
        _a = sympy.lambdify((w), a)
        _b = sympy.lambdify((w), (b.as_numer_denom()[0])**2)    
        
        real, imag = L.as_real_imag()
        phase = -180/pi* sympy.atan2(imag,real)
        magL = sympy.sqrt(real**2+imag**2)
        f1 = (-180 - phase)**2
        f2 = (1 - magL)**2
        _f1 = sympy.lambdify((w),f1)
        _f2 = sympy.lambdify((w),f2)
        L,U = 0.001*2*pi, 10  #lower and upper bounds on w180cross frequency, rad/s
        I, w1, aa, bb = lineSearch(L,U,0.01,_b)
        I, w2, aa, bb = lineSearch(L,U,0.01,_f2)
        
        return w1, w2
        
    def _w180cross(self, m):    
        _real, _imag = self._ForwardPathTF()    
        return 0 == _imag(m.Kp,m.Ki,m.w180cross)
    
    def _w0dBcross(self, m):    
        _real, _imag = self._ForwardPathTF()    
        return 1 == _real(m.Kp,m.Ki,m.w0dBcross)**2 + _imag(m.Kp,m.Ki,m.w0dBcross)**2
    
    def _gainMargindB(self, m):    
        _real, _imag = self._ForwardPathTF()    
        magL = sqrt(_real(m.Kp,m.Ki,m.w180cross)**2+_imag(m.Kp,m.Ki,m.w180cross)**2)
        return 20*log10(magL) <= -10
    
    def getActivationFunction(self, lambda_low, tau_high, epsilon):
            n = 0.01        
            a = 0.5 * ((epsilon - lambda_low) / ( (epsilon-lambda_low)**2 + n**2 )**0.5 + \
               (tau_high - epsilon) / ( (tau_high-epsilon)**2 + n**2 )**0.5 )        
            return a
    
    def _atan2(self, y, x):    
        q1 = self.getActivationFunction(-1e6, 0, x)
        q2 = pi*(self.getActivationFunction(0, 1e6, y)*2-1)
        angle = atan(y/x) + q1*q2
        return angle
    
    def _phaseMargin(self, m):    
        _real, _imag = self._ForwardPathTF()    
        phaseL = 180/pi*self._atan2(_imag(m.Kp,m.Ki,m.w0dBcross),_real(m.Kp,m.Ki,m.w0dBcross))
        return phaseL >= -100
        
    def _Jdot(self, m, t):    
        tracking_error = (m.r[t] - m.y1[t])**2
        return m.dJdt[t] == m.tf * tracking_error

def main():
    # construct gain scheduler object
    gs = GainScheduler()
    
    Kp_guess = 3
    Ki_guess = 8
    
    # Vars
    gs.m.tf = Param(initialize = 30)
    gs.m.Kp = Var(initialize = Kp_guess, bounds=(1,10))
    gs.m.Ki = Var(initialize = Ki_guess, bounds=(0,10))
    gs.m.y1 = Var(gs.m.tau)      
    gs.m.y2 = Var(gs.m.tau) 
    gs.m.y3 = Var(gs.m.tau) 
    gs.m.u  = Var(gs.m.tau) 
    gs.m.r  = Var(gs.m.tau)
    gs.m.J  = Var(gs.m.tau)
    
    gs.m.w180cross = Var(initialize = 1, bounds=(0,10))
    gs.m.w0dBcross = Var(initialize = 1.2, bounds=(0,10))
    
    # Derivative Vars 
    gs.m.dy1dt     = DerivativeVar(gs.m.y1, wrt=gs.m.tau)
    gs.m.dy2dt     = DerivativeVar(gs.m.y2, wrt=gs.m.tau)
    gs.m.dy3dt     = DerivativeVar(gs.m.y3, wrt=gs.m.tau)
    gs.m.dudt      = DerivativeVar(gs.m.u, wrt=gs.m.tau)
    gs.m.drdt      = DerivativeVar(gs.m.r, wrt=gs.m.tau)
    gs.m.dJdt      = DerivativeVar(gs.m.J, wrt=gs.m.tau)
    
    ### Discretize problem (after m.tau and all variables are defined, before constraint defs defined)       
    discretizer = TransformationFactory('dae.collocation')        
    discretizer.apply_to(gs.m, wrt=gs.m.tau, nfe=75, ncp=10, scheme='LAGRANGE-RADAU') 
    
    #% Build model
        
    
    # Phase and Gain crossing constraints
    gs.m.w180crossConstraint            = Constraint(rule=gs._w180cross)
    gs.m.w0dBcrossConstraint            = Constraint(rule=gs._w0dBcross)
    
    # Standard constraints
    gs.m.bcConstraint         = ConstraintList(rule=gs.BCs)
    gs.m.y1dotConstraint      = Constraint(gs.m.tau, rule=gs._y1dot)
    gs.m.y2dotConstraint      = Constraint(gs.m.tau, rule=gs._y2dot)
    gs.m.y3dotConstraint      = Constraint(gs.m.tau, rule=gs._y3dot)
    gs.m.udotConstraint      = Constraint(gs.m.tau, rule=gs._udot)
    gs.m.rdotConstraint      = Constraint(gs.m.tau, rule=gs._rdot)
    
    # Margin constraints
    gs.m.gainMargindBConstraint   = Constraint(rule=gs._gainMargindB)
    gs.m.phaseMarginConstraint   = Constraint(rule=gs._phaseMargin)
    
    gs.m.JdotConstraint      = Constraint(gs.m.tau, rule=gs._Jdot)
    
    gs.m.obj  = Objective(expr = gs.m.J[1], sense=minimize)
    
    
    solver = SolverFactory('ipopt')
    solver.options["halt_on_ampl_error"] = "yes"
    solver.options["tol"] = 1e-6
    solver.options["max_iter"] = 2000
    solver.options["linear_solver"] = "ma27"
    solver.options["ma27_pivtol"] = 1e-6
    #solver.options["ma27_liw_init_factor"] = 5
    solver.options["linear_scaling_on_demand"] = "yes"
    solver.options["mu_init"] = 1
    solver.solve(gs.m, tee=True)
    
    results_dict = {}
    var_list = [gs.m.Kp, gs.m.Ki]
    for var in var_list:
        print("Saving Variable: ",var) # doctest: +SKIP
        results_dict[var.name] = [value(var[index]) for index in var]
    
    #%% Plotting
    mpl.rcParams['font.size'] = 14
    mpl.rcParams['legend.fontsize'] = 14
    mpl.rcParams["font.family"] = "Times New Roman"
    csfont = {'fontname':'Times New Roman'}
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    #plt.close('all')
    
    blue = '#1f77b4ff'
    orange = '#ff7f0eff'
    green = '#2ca02cff'
     
    time = np.dot(gs.m.tau,gs.m.tf())
    y1 = [gs.m.y1[t]() for t in gs.m.tau]
    r = [gs.m.r[t]() for t in gs.m.tau]
    plt.figure(1)
    plt.plot(time, y1,'.-')
    plt.plot(time, r,'k')
    plt.xlabel('t [s]')
    plt.ylabel('y')
    plt.grid()
    plt.title(['Kp:%6.2f' %gs.m.Kp(), 'Ki:%6.2f' %gs.m.Ki(), 'J = %6.2f' %gs.m.J[1]()])
    
    print('Kp = ', gs.m.Kp())
    print('Ki = ', gs.m.Ki())
    
    sysC = control.tf([gs.m.Kp(), gs.m.Ki()],[1,0])
    sysG = control.tf([1],[1,4,8,8])
    sys = sysC * sysG
    rlist, klist = control.rlocus(sys, kvect=np.logspace(-10,1,500)) 
    
    sysC = control.tf([gs.m.Kp(), gs.m.Ki()],[1,0])
    sysG = control.tf([1],[1,4,8,8])
    sys = sysC * sysG
    r,k = control.rlocus(sys,[1])
    plt.close()
    #fig,ax = plt.subplots()
    plt.scatter(r.real,r.imag, color='k')
    
    sysC = control.tf([gs.m.Kp(), gs.m.Ki()],[1,0])
    sysG = control.tf([1],[1,4,8,8])
    sys = sysC * sysG
    plt.figure()
    mag, phase, omega = control.bode(sys,dB=True, margins=True)
    
    
    
    #%% Check gain/phase
    
    #K  = m.K()
    #s  = complex(real=0, imag=m.w180cross())
    #CG = (s+1)/(s*(s+2)*(s**2+2*s+4))
    #L  = K*CG
    #magL = sqrt(L.real**2+L.imag**2)
    #angleL = -180/pi*np.arctan2(L.imag,L.real)
    
    
    w = sympy.Symbol('w', real=True)
    s = sympy.I*w
    C = gs.m.Kp() + gs.m.Ki()/s
    G = 1/((s+2)*(s**2+2*s+4))
    L = C*G
    real, imag = L.as_real_imag()
    phase = -180/pi* sympy.atan2(imag,real)
    mag = sympy.sqrt(real**2+imag**2)
    
    _phase = sympy.lambdify((w),phase)
    _mag = sympy.lambdify((w),mag)
    
    print('w180cross = ', gs.m.w180cross(), 'rad/s')
    print('w0dBcross = ', gs.m.w0dBcross(), 'rad/s')
    print('magLdB at w0dBcross = ', 20*log10(_mag(gs.m.w0dBcross())))
    print('angleL at w180cross = ', _phase(gs.m.w180cross()))
    
    return

if __name__ == '__main__':
    main()


