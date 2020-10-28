# -*- coding: utf-8 -*-

"""
Created on Mon Aug 26 13:24:03 2019

@author: kwilli2
"""
def adaptiveMeshFunc(nodes,my_deriv_vars,m):
    import numpy as np
    import matplotlib.pyplot as plt
    
    tau = m.tau
    M = 0.01
    U = 10
    flag = 0
    while flag == 0:
        print('M=',M) #Comment this out for speed
        print('U=',U)
        
        ## Create cumulative distribution function Fx
        myStack = np.zeros([len(tau)])
        for dxdt in my_deriv_vars:
            
            x = dxdt.get_state_var()
            #print(x.name)
            
            temp = abs(np.array(x[:]()))
            myMean = np.std(temp)+0.001  #absolute value mean, plus small  value to prevent zero
            
            # Compute regularized density increment
            myStackIncr = np.array([max(0,abs(dxdt[t]())/myMean*abs(tau.prev(max(tau[2],t))-t)) for t in tau]) 
            myStack = myStack + myStackIncr
            
            
        fx = myStack/sum(myStack)  
        Fx = np.cumsum(fx)
        Fx[-1] = 1; #Final value of Fx may be slightly different than 1 due to numerical precision, force to 1.0
        
        ## Regularization, works by adjusting slope of Fx
        e = M # Minimum allowed Fx slope 
        E = U # Maximum allowed Fx slope
        Fx = [0*i for i in range(len(tau))]    
        Fx[0] = fx[0]
        for i in range(len(tau)-1):
            Fx[i+1] = Fx[i]+fx[i+1]
            if (Fx[i+1]-Fx[i])/(list(tau.value)[i+1]-list(tau.value)[i]) < e:
                Fx[i+1] = e*(list(tau.value)[i+1]-list(tau.value)[i])+Fx[i]
            if (Fx[i+1]-Fx[i])/(list(tau.value)[i+1]-list(tau.value)[i]) > E:
                Fx[i+1] = E*(list(tau.value)[i+1]-list(tau.value)[i])+Fx[i]
        Fx = Fx/Fx[-1]            
        Fx[-1] = 1; 
        
        
        ## Distribute nodes according to inverted Fx 
        y = np.zeros(nodes)
        u = np.linspace(0,1,nodes)
        for i in range(nodes):
            y[i] = np.interp(u[i], Fx, np.array(list(tau.value)))
        y[-1] = 1; y[0] = 0;  #Numerical precision may produce values different from 0 and 1 at endpoints.  Force 0 and 1 at endpoints.
        
        ## Check sparsity/density of nodes, set flag to 1 if sparsity/density level OK, or if maximum value of regularization term met
        flag=1 # by default, prepare to terminate the while loop
        maxdiff=0
        mindiff=np.inf
        for i in range(nodes-1):
            thediff = y[i+1]-y[i]
            if thediff > maxdiff:
                maxdiff = thediff
            if thediff < mindiff:
                mindiff = thediff                    
        if maxdiff > 1.5/(nodes-1): #Heuristic: Do not allow finite elements to be further apart than 1.5 times uniform distribution
            if M > 4: #if maximum value of regularization term met, terminate while loop
                flag = 1
            else:
                M = M*2 #Increment regularization term
                flag=0
        if mindiff < 1/(nodes-1)/200:
            if U < M: 
                flag = 1
            else:
                U = U*0.9 #Increment regularization term
                flag=0        
    
    # Plot density function and cumulative distribution function
    plt.figure(101)      
    plt.subplot(2, 1, 1)
    plt.title('Node distribution along indices')  
    # plt.plot(list(tau.value),fx,'o'); plt.xlabel('tau')
    plt.plot(list(range(len(tau))),fx,'o'); plt.xlabel('node i')
    plt.ylabel('fx')
    plt.grid(True)
    
    plt.subplot(2, 1, 2)
    # plt.plot(list(tau.value),Fx,'o'); plt.xlabel('tau')
    plt.plot(list(range(len(tau))),Fx,'o'); plt.xlabel('node i')
    plt.plot(list(range(len(tau))),np.linspace(0,1,len(tau)),'k--'); plt.xlabel('node i')
    plt.ylabel('Fx')
    plt.grid(True)
    
    return y

