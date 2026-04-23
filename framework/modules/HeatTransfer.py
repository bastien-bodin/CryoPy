import numpy as np
from math import exp, erf, sqrt, pi

from pysph.sph.equation import Equation
from compyle.api import declare

class MassCluster(Equation):
    def initialize(self, d_idx, d_MassCluster, d_m):
        d_MassCluster[d_idx] = d_m[d_idx]
    
    def loop(self, d_idx, s_m, RIJ, s_idx, d_MassCluster):
        if RIJ > 1.0e-12:
            d_MassCluster[d_idx] += s_m[s_idx]

class CenterOfMass(Equation):
    def __init__(self, dest, sources, dim=2, x0=0.0):
        self.x0 = x0
        if dim < 1:
            dim=1
        elif dim > 3:
            dim=3
        if dim == 1:
            self.xi = 4
        elif dim == 2:
            self.xi = 13
        else:
            self.xi = 44
        super(CenterOfMass, self).__init__(dest, sources)
    
    def initialize(self, d_idx, d_Surface, d_x_com, d_y_com):
        d_x_com[d_idx] = 0.0
        d_y_com[d_idx] = 0.0
        d_Surface[d_idx] = 0.0
    
    def loop_all(self, d_idx, d_Surface, d_MassCluster,
                 s_m, d_x_com, d_y_com, d_x, d_y, s_x, s_y,
                 d_h, N_NBRS, NBRS):
        i = declare('int')
        s_idx = declare('long')
        xij = declare('matrix(2)')
        r_com = 0.0
        a = 0.0
        # Compute the CoM vector
        for i in range(N_NBRS):
            s_idx = NBRS[i]
            xij[0] = d_x[d_idx] - s_x[s_idx]
            xij[1] = d_y[d_idx] - s_y[s_idx]
            rij = (xij[0]*xij[0] + xij[1]*xij[1])**0.5
            if rij > 1e-12:
                d_x_com[d_idx] += (
                    s_m[s_idx] * xij[0] / d_MassCluster[d_idx]
                )
                d_y_com[d_idx] += (
                    s_m[s_idx] * xij[1] / d_MassCluster[d_idx]
                )
        r_com = (
            d_x_com[d_idx]*d_x_com[d_idx] +
            d_y_com[d_idx]*d_y_com[d_idx]
        )**0.5
        # Criterion to select if the particle belongs to the
        # surface or not
        condition_com = (r_com >= 0.45*d_h[d_idx])
        condition_nei = (N_NBRS <= self.xi)
        condition = (condition_com and condition_nei)
        if condition and (d_x[d_idx] > self.x0 - 5*d_h[d_idx]):
            d_Surface[d_idx] = 1.0


class UpdateThermalProps(Equation):
    def __init__(self, dest, sources, ParamsLiq:tuple, ParamsSol:tuple):
        """We compute averaged thermal properties based on their
        massic ice fraction (as the density stays the same, it is
        the same as their volumic ice fraction).
        
        Every property for the liquid or solid should be a tuple:
        `(kappa, cp)`.

        Args:
            ParamsLiq (tuple): Thermal properties of the liquid lava.
            ParamsSol (tuple): Thermal properties of the solid lava.
        """
        self.kliq, self.cliq = ParamsLiq
        self.ksol, self.csol = ParamsSol
        super(UpdateThermalProps, self).__init__(dest, sources)
    
    def initialize(self, d_idx, d_qstate, d_cp, d_rho, d_kappa):
        kliq, cliq = self.kliq, self.cliq
        ksol, csol = self.ksol, self.csol
        q = d_qstate[d_idx]
        q1 = 1.0 - q
        d_kappa[d_idx] = kliq*q1 + ksol*q
        d_cp[d_idx] = cliq*q1 + csol*q
        # if d_qstate[d_idx] == 1.0:
        #     d_kappa[d_idx] = ksol
        #     d_cp[d_idx] = csol
        #     d_rho[d_idx] = rsol
        # else:
        #     d_kappa[d_idx] = kliq
        #     d_cp[d_idx] = cliq
        #     d_rho[d_idx] = rliq


class ConductionMonaghan(Equation):
    def initialize(self, d_idx, d_aTemp):
        d_aTemp[d_idx] = 0.0
    
    def loop(self, d_idx, d_Cond, d_kappa, d_rho, d_cp, d_Temp,
             s_idx, s_kappa, s_rho, s_Temp, s_m, d_aTemp,
             DWIJ, R2IJ, XIJ, EPS):
        dot = DWIJ[0]*XIJ[0] + DWIJ[1]*XIJ[1]
        Fij = dot/(R2IJ + EPS)
        Vj = s_m[s_idx]/s_rho[s_idx]
        Kij = (4*d_kappa[d_idx]*s_kappa[s_idx]/
               (d_kappa[d_idx]+s_kappa[s_idx]))
        Tij = d_Temp[d_idx] - s_Temp[s_idx]
        d_aTemp[d_idx] += Vj*Kij*Tij*Fij/d_cp[d_idx]/d_rho[d_idx]

class Conduction(Equation):
    def initialize(self, d_idx, d_Cond):
        d_Cond[d_idx] = 0.0
    
    def loop(self, d_idx, d_Cond, d_kappa, d_rho, d_cp, d_Temp,
             s_idx, s_kappa, s_rho, s_Temp, s_m, d_aTemp,
             DWIJ, R2IJ, XIJ, EPS):
        dot = DWIJ[0]*XIJ[0] + DWIJ[1]*XIJ[1]
        Fij = dot/(R2IJ + EPS)
        Vj = s_m[s_idx]/s_rho[s_idx]
        Kij = (4*d_kappa[d_idx]*s_kappa[s_idx]/
               (d_kappa[d_idx]+s_kappa[s_idx]))
        Tij = d_Temp[d_idx] - s_Temp[s_idx]
        d_Cond[d_idx] += Vj*Kij*Tij*Fij/d_cp[d_idx]/d_rho[d_idx]

class Radiation(Equation):
    def __init__(self, dest, sources, Tenv=93.0, eps = 1.0):
        self.Tenv = Tenv
        self.eps = eps
        super(Radiation, self).__init__(dest, sources)
    
    def initialize(self, d_idx, d_Rad, d_Temp, d_Surface):
        if d_Surface[d_idx] == 1.0:
            eps = self.eps
            sigma = 5.67e-8     # Stefan-Boltzmann constant
            d_Rad[d_idx] = eps*sigma*(
                d_Temp[d_idx]**4.0 - self.Tenv**4.0
            )
        else:
            d_Rad[d_idx] = 0.0

class Convection(Equation):
    def __init__(self, dest, sources, Tenv=93.0,
                 Uwind = 0.0, g = 1.352):
        self.Tenv = Tenv
        self.Uwind = Uwind
        self.g = abs(g)
        super(Convection, self).__init__(dest, sources)
    
    def initialize(self, d_idx, d_Conv, d_Temp, d_Surface):
        if d_Surface[d_idx] == 1.0:
            dT = d_Temp[d_idx] - self.Tenv
            Teff = 0.5*(d_Temp[d_idx] + self.Tenv)
            eta_atm = (0.0749*Teff - 0.00005*Teff**2.0) # in muPa
            cp_atm  = (
                1.9327e-10 * Teff**4.0 -
                7.9999e-7  * Teff**3.0 +
                0.0011407  * Teff**2.0 -
                0.4489     * Teff      +
                1087.5
            )
            rho_atm = 4202.022/(8.3145*Teff)
            k_atm   = (
                1.5207e-11 * Teff**3.0 -
                4.8574e-8  * Teff**2.0 +
                1.0184e-4  * Teff      -
                0.00039333
            )
            beta    = 1.0/Teff
            f       = 0.0036
            A       = 0.14
            f1 = 100.0*(cp_atm*self.g*beta/eta_atm)**(1.0/3.0)
            f2 = (rho_atm*k_atm)**(2.0/3.0)
            f3 = dT**(4.0/3.0)
            F_force = self.Uwind * f * dT * rho_atm * cp_atm
            F_nat   = A*f1*f2*f3
            d_Conv[d_idx] = max(F_force, F_nat)
        else:
            d_Conv[d_idx] = 0.0

class HeatTransfer(Equation):
    def __init__(self, dest, sources, dp):
        self.dp = dp
        super(HeatTransfer, self).__init__(dest, sources)

    def initialize(self, d_idx, d_aTemp,
                   d_Cond, d_Rad, d_Conv, d_m, d_cp):
        dS = self.dp
        d_aTemp[d_idx] = d_Cond[d_idx] - (
            d_Conv[d_idx] + d_Rad[d_idx]
        )*dS/(d_m[d_idx]*d_cp[d_idx])
