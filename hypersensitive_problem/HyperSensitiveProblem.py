from pyomo.environ import *
from pyomo.dae import *
from math import pi
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from adaptiveMesh import adaptiveMeshFunc

for iter in range(4):
    warmStart = 1
    cubicSplineExtrapU0 = 1
    
    # myMethod = 'ORTHO'; nfe = 20; ncp = 11  
    myMethod = 'ORTHO'; nfe = 50; ncp = 9  
    # myMethod = 'ORTHO'; nfe = 100; ncp = 5
    
    if 'run_count' not in locals():
        run_count = 0
    run_count = run_count + 1
    
    reMesh = np.linspace(0,1,nfe+1)
    if 'myProfile' in locals():
        my_deriv_vars = [m.dxdt]
        reMesh = adaptiveMeshFunc(nfe+1,my_deriv_vars,m)            
    
    plt.figure(99)
    plt.plot(reMesh,run_count*np.ones(len(reMesh)),'-o',MarkerSize=10,alpha=0.5)
    plt.ylim([0,run_count+1])
    plt.ylabel('Iteration')
    plt.xlabel('End location of Finite Element $n$')
    
    if 'soln_time' not in locals():
        soln_time = 0
    
    model = m = ConcreteModel()
    m.tau = ContinuousSet(initialize=reMesh)    
       
    ### Params    
    m.tf = Param(initialize=10000)    
    # scaling
    m.x_scale = Param(initialize = 1e-1)
    m.u_scale = Param(initialize = 1e-1)        
    m.y_scale = Param(initialize = 1e0)    
    # IC's and TC's
    m.x_0 = Param(initialize = 1/m.x_scale)
    m.y_0 = Param(initialize = 0/m.y_scale)
    m.x_1 = Param(initialize = 1.5/m.x_scale)    
    # Control limits
    m.uBound = Param(initialize = 10)        
    
    ### Vars
    # controlled variables (states)
    m.x = Var(m.tau, within=Reals)
    m.y = Var(m.tau, within=NonNegativeReals)    
    # manipulated variables (controls)
    m.u = Var(m.tau, bounds = (-m.uBound/m.u_scale,m.uBound/m.u_scale))
    
    ### DerivativeVars
    m.dxdt = DerivativeVar(m.x, wrt=m.tau)
    m.dydt = DerivativeVar(m.y, wrt=m.tau)
    
    
    ### Initial and Terminal conditions
    @m.ConstraintList()
    def BCs(m):
        yield m.x[0] == m.x_0
        yield m.x[1] == m.x_1
        yield m.y[0] == m.y_0                   
    
    
    ### Differential equations
    @m.Constraint(m.tau)
    def xdot(m, t):
        x = m.x_scale * m.x[t]
        u = m.u_scale * m.u[t]
        return m.x_scale * m.dxdt[t] == m.tf * (-x**3 + u)
    
    @m.Constraint(m.tau)
    def ydot(m, t):
        x = m.x_scale * m.x[t]
        u = m.u_scale * m.u[t]
        return m.y_scale * m.dydt[t] == m.tf * 0.5 * ( x**2 + u**2 )
    
    
    m.obj = Objective(expr=m.y[1], sense=minimize)
    
    ### Discretize problem        
    if myMethod == 'ORTHO':        
        discretizer = TransformationFactory('dae.collocation')        
        discretizer.apply_to(m,wrt=m.tau, nfe=nfe, ncp=ncp, scheme='LAGRANGE-RADAU')   
    else:
        discretizer = TransformationFactory('dae.finite_difference')
        discretizer.apply_to(m, wrt=m.tau, nfe=nfe, scheme='backward')  
    
    #discretizer.reduce_collocation_points(m, var=m.alpha, ncp=1, contset=m.tau)
    
    vals = np.linspace(0/m.u_scale.value, 0/m.u_scale.value, len(m.tau))
    m.u.set_values(dict(zip(m.tau.value, vals)))
    
    vals = np.linspace(m.x_0.value, m.x_1.value, len(m.tau))
    m.u.set_values(dict(zip(m.tau.value, vals)))
    
    # warm-start    
    if 'oldtau' in locals() and warmStart == 1:            
        for myVar in m.component_objects(Var, active= True):
            print("Initializing Variable: ",myVar)
            if len(myVar) == 1:
                myVar.set_value(myProfile[myVar.name][0])
            else:
                myVar.set_values(dict(zip(list(m.tau), np.interp(m.tau, oldtau, myProfile[myVar.name]))))
        
   
    ### Declare all suffixes
    solver = SolverFactory('ipopt')
    solver.options["halt_on_ampl_error"] = "yes"
    solver.options["tol"] = 1e-6
    solver.options["max_iter"] = 2000
    solver.options["linear_solver"] = "ma27"
    solver.options["ma27_pivtol"] = 1e-6
    #solver.options["ma27_liw_init_factor"] = 5
    solver.options["linear_scaling_on_demand"] = "yes"
    if 'oldtau' in locals() and warmStart == 1:
        solver.options["mu_init"] = 1 #The larger mu_init is, the less ipopt will "use" the warm start guess
    
    ### FIRST SOLVE
    results = solver.solve(m, tee=True, keepfiles=True, logfile=(str("ipopt_%s.log" % nfe )))
    soln_time = results['Solver'][0]['Time'] + soln_time
    
    print("\n")
    print('Mayer-Type Objective Function Value = ',m.y[1].value*m.y_scale.value)
    print('Total Solution Time = ',soln_time)
    print('Solution Time = ',results['Solver'][0]['Time'])
    print("\n")
    
    # Recover initial control through cubic spline extrapolation
    if cubicSplineExtrapU0 == 1:
        t = [list(m.tau.value)[i] for i in np.arange(1,5)]
        u = [m.u[list(m.tau.value)[i]].value for i in np.arange(1,5)]
        z = np.polyfit(t, u, 3)
        f = np.poly1d(z)
        m.u[0].value = f(list(m.tau.value)[0])
    
    # For Warm Start
    oldtau      = m.tau            
    myProfile = {}
    for myVar in m.component_objects(Var, active= True):
        print("Saving Variable: ",myVar) 
        for index in myVar:            
            myProfile[myVar.name] = [value(myVar[index]) for index in myVar] 
    
    #%% Plotting
    myLegend = ['Uniform mesh', 'Adapted mesh']    
    plt.rcParams['font.size'] = 14
    
    plt.figure(1)
    
    ax1 = plt.subplot(2,2,1)
    ax1.plot(np.dot(m.tf.value,list(m.tau.value)), [m.x_scale.value*value(m.x[tau]) for tau in m.tau],'o-')
    plt.xlim([0, 25])
    plt.ylabel('$x(t)$')
    plt.xlabel('t')
    plt.legend(myLegend)
    plt.title('Near $t=0$')
    plt.grid(True,linestyle='--')
    
    ax2 = plt.subplot(2,2,2)        
    ax2.plot(np.dot(m.tf.value,list(m.tau.value)), [m.x_scale.value*value(m.x[tau]) for tau in m.tau],'o-')
    plt.xlim([9975, 10000])
    plt.grid(True,linestyle='--')
    plt.ylabel('$x(t)$')
    plt.xlabel('t')
    plt.legend(myLegend)
    plt.title('Near $t=t_f$')
    plt.grid(True,linestyle='--')
    
    ax3 = plt.subplot(2,2,3,sharex = ax1)
    ax3.plot(np.dot(m.tf.value,list(m.tau.value)), [m.u_scale.value*value(m.u[tau]) for tau in m.tau],'o-')
    plt.ylim([-1, 1])
    plt.ylabel('$u(t)$')
    plt.xlabel('t')
    plt.legend(myLegend)
    plt.grid(True,linestyle='--')
    
    ax4 = plt.subplot(2,2,4,sharex = ax2)
    ax4.plot(np.dot(m.tf.value,list(m.tau.value)), [m.u_scale.value*value(m.u[tau]) for tau in m.tau],'o-')
    plt.ylabel('$u(t)$')
    plt.xlabel('t')
    plt.legend(myLegend)
    plt.grid(True,linestyle='--')
    
    plt.figure(2)
    plt.plot(np.dot(m.tf.value,list(m.tau.value)), [m.x_scale.value*value(m.x[tau]) for tau in m.tau],'o-')        
    plt.grid(True,linestyle='--')
    plt.ylabel('$x(t)$')
    plt.xlabel('t')
    plt.legend(myLegend)
    plt.grid(True,linestyle='--')
        
        
