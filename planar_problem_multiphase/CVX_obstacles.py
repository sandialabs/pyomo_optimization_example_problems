from math import pi
from pyomo.environ import *
from pyomo.dae import *

RAD_TO_DEG = 180 / pi

from utilities.getPhaseVariables import getPhaseVariables


# Defines ellipses and the constraints to avoid the ellipses
class CVX_obstacles():

    def __init__(self, m, obstacles_present=False):

        # Convex Obstacle Geometries
        m.ell1Width = Param(initialize=0.5)
        m.ell1Height = Param(initialize=0.5)
        m.ell1CenterX = Param(initialize=3)
        m.ell1CenterY = Param(initialize=1.5)


        if obstacles_present:

            # No-fly-zone 1
            self.a = m.ell1Width
            self.b = m.ell1Height
            self.c = (m.ell1CenterX, m.ell1CenterY)
            
            # Add constraints to model
            phase = 1
            m.CVX1a = Constraint(m.tau1, [phase], rule=self.ellipseX)

            phase = 2
            m.CVX1b = Constraint(m.tau2, [phase], rule=self.ellipseX)

        return

    ### CVX obstacle
    def ellipseX(self, m, t, n):

        (tf, x, y, theta, v, omega, J, dx_dtau, dy_dtau, \
         dtheta_dtau, dv_dtau, domega_dtau, dJ_dtau, \
         F_control, T_control, F_controldot, T_controldot, \
         F_virtual, T_virtual, dF_control_dtau, dT_control_dtau) = getPhaseVariables(m, n, t)  

        return 1 / (self.b / 2) ** 2 * (y - self.c[1]) ** 2 + 1 / (self.a / 2) ** 2 * (x - self.c[0]) ** 2 >= 1

