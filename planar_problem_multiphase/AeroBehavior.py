# -*- coding: utf-8 -*-

"""

@author: Rachel Schlossman

"""
class AeroBehavior():
    def __init__(self):
        
        return
    
    def getActivationFunction(self, lambda_low, tau_high, epsilon):
        n = 0.01        
        a = 0.5 * ((epsilon - lambda_low) / ( (epsilon-lambda_low)**2 + n**2 )**0.5 + \
           (tau_high - epsilon) / ( (tau_high-epsilon)**2 + n**2 )**0.5 )        
        return a
   
    # coefficients 
    def getCvCvw(self, x):
    
        
        bound0 = 0
        bound1 = 10
        bound2 = 20
        bound3 = 10000
        
        epsilon = abs(x)
        
        a1 = self.getActivationFunction(bound0, bound1, epsilon)
        a2 = self.getActivationFunction(bound1, bound2, epsilon)
        a3 = self.getActivationFunction(bound2, bound3, epsilon)
        
        cv = -0.1 * a1 - 0.2 * a2 - 0.3 * a3
        cvw = -0.3 * a1 - 0.2 * a2 - 0.1 * a3
        
        return (cv, cvw)