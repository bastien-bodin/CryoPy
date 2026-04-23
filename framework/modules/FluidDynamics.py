import numpy as np
from math import exp, erf, sqrt, pi

from pysph.sph.equation import Equation
from compyle.api import declare

class GradientCorrection(Equation):
    r"""**Kernel Gradient Correction**

    From [BonetLok1999], equations (42) and (45)

    .. math::
            \nabla \tilde{W}_{ab} = L_{a}\nabla W_{ab}

    .. math::
            L_{a} = \left(\sum \frac{m_{b}}{\rho_{b}} \nabla W_{ab}
            \mathbf{\otimes}x_{ba} \right)^{-1}
    """

    def __init__(self, dest, sources, dim=2, tol=0.1):
        self.dim = dim
        self.tol = tol
        super(GradientCorrection, self).__init__(dest, sources)
        
    def loop(self, d_idx, d_m_mat, DWIJ, HIJ):
        i, j, n, nt = declare('int', 4)
        n = self.dim
        nt = n + 1
        # Note that we allocate enough for a 3D case but may only use a
        # part of the matrix.
        temp = declare('matrix(12)')
        res = declare('matrix(3)')
        eps = 1.0e-04 * HIJ
        for i in range(n):
            for j in range(n):
                temp[nt * i + j] = d_m_mat[9 * d_idx + 3 * i + j]
            # Augmented part of matrix
            temp[nt*i + n] = DWIJ[i]

        res_mag = 0.0
        dwij_mag = 0.0
        for i in range(n):
            res_mag += abs(res[i])
            dwij_mag += abs(DWIJ[i])
        change = abs(res_mag - dwij_mag)/(dwij_mag + eps)
        if change < self.tol:
            for i in range(n):
                DWIJ[i] = res[i]


class TaitEoS(Equation):
    """Tait Equation of State for Weakly-Compressible SPH (WCSPH):
    Computes pressure `p` from density `rho`.
    
    Args:
        rho0 (float): Base density of SPH particles.
        c0 (float): Sound speed for SPH particles. By convention it is put at 10*Vmax.
        gamma (float, optional): gamma factor in EoS. Default to 7.0.
    """
    def __init__(self, dest, sources, rho0, c0, gamma=7.0):
        self.rho0=rho0
        self.c0 = c0
        self.gamma = gamma
        self.B = rho0*c0*c0/gamma
        super(TaitEoS, self).__init__(dest, sources)
    
    def initialize(self, d_idx, d_p, d_rho, d_cs):
        tmp1 = pow(d_rho[d_idx]/self.rho0, self.gamma)
        tmp2 = pow(d_rho[d_idx]/self.rho0, self.gamma - 1.0)
        d_p[d_idx] = self.B * (tmp1 - 1.0)
        d_cs[d_idx] = self.c0 * pow(tmp2, 0.5)

class ContinuityEquation(Equation):
    """Default Continuity equation.

    Args:
        Equation (_type_): _description_
    """
    def initialize(self, d_arho, d_idx):
        d_arho[d_idx] = 0.0
    
    def loop(self, d_arho, d_idx, d_rho, s_rho, s_m, s_idx, VIJ, DWIJ):
        dot = VIJ[0]*DWIJ[0] + VIJ[1]*DWIJ[1]
        Vj = s_m[s_idx]/s_rho[s_idx]
        d_arho[d_idx] += d_rho[d_idx] * Vj * dot

class ContinuityEquationLANS(Equation):
    """Continuity equation for turbulent flows, using the SPH-epsilon
    technique of Monaghan.
    
    Ref:
        Monaghan, J. J. (2011). A turbulence model for smoothed particle hydrodynamics.
        European Journal of Mechanics-B/Fluids, 30(4), 360-370.
    
    """
    def initialize(self, d_arho, d_idx):
        d_arho[d_idx] = 0.0
    
    def loop(self, d_arho, d_idx, d_rho, s_rho, s_m, s_idx,
             d_usph, d_vsph, s_usph, s_vsph,VIJ, DWIJ):
        uij = VIJ[0] + (d_usph[d_idx]-s_usph[s_idx])
        vij = VIJ[1] + (d_vsph[d_idx]-s_vsph[s_idx])
        dot = uij*DWIJ[0] + vij*DWIJ[1]
        Vj = s_m[s_idx]/s_rho[s_idx]
        d_arho[d_idx] += s_m[s_idx]*dot#d_rho[d_idx] * Vj * dot

class VelGrad2D(Equation):
    """Computes the 2D velocity gradient with the SPH formalism
    """
    def initialize(self, d_idx, d_v00, d_v01, d_v10, d_v11):
        d_v00[d_idx] = 0.0
        d_v01[d_idx] = 0.0
        d_v10[d_idx] = 0.0
        d_v11[d_idx] = 0.0
    
    def loop(self, d_idx, s_idx, s_m, s_rho,
             d_v00, d_v01, d_v10, d_v11, VIJ, DWIJ):
        Vj = s_m[s_idx]/s_rho[s_idx]
        d_v00[d_idx] += -VIJ[0]*DWIJ[0]*Vj
        d_v01[d_idx] += -VIJ[0]*DWIJ[1]*Vj
        d_v10[d_idx] += -VIJ[1]*DWIJ[0]*Vj
        d_v11[d_idx] += -VIJ[1]*DWIJ[1]*Vj

class ViscosityBinghamRho(Equation):
    def __init__(self, dest, sources, rho0=500.0, mp=1000.0):
        self.rho0 = rho0
        self.mp = mp
        super(ViscosityBinghamRho, self).__init__(dest, sources)
    
    def initialize(self, d_idx, d_mu, d_v00, d_v01, d_v10, d_v11,
                   d_g00, d_g01, d_g10, d_g11, d_gmag):
        # computes the shear rate tensor
        d_g00[d_idx] = d_v00[d_idx]
        d_g01[d_idx] = 0.5*(d_v01[d_idx] + d_v10[d_idx])
        d_g10[d_idx] = 0.5*(d_v10[d_idx] + d_v01[d_idx])
        d_g11[d_idx] = d_v11[d_idx]
        
        IIe = d_g00[d_idx]**2.0 + d_g01[d_idx]**2.0 + d_g10[d_idx]**2.0 + d_g11[d_idx]**2.0
        a = 2.88e-4
        b = 1.42e-2
        c = 36.72e-3
        d = 4.92e-3
        mu_0 = a*exp(b*self.rho0)
        tau_0 = c*exp(d*self.rho0)
        if IIe>1e-8:
            d_gmag[d_idx] = ((2.0*IIe)**0.5)
            d_mu[d_idx] = mu_0 + tau_0*(1.0 - exp(-self.mp*d_gmag[d_idx]))/d_gmag[d_idx]
        else:
            d_gmag[d_idx] = 0.0
            d_mu[d_idx] = mu_0 + tau_0*self.mp

class ViscosityCross(Equation):
    def __init__(self, dest, sources,
                 nu0=2.1,nu1=2.7e-3,kc=1.1,rho0=500.0):
        self.rho0 = rho0
        self.nu0 = nu0
        self.nu1=nu1
        self.kc=kc
        super(ViscosityCross, self).__init__(dest, sources)
    
    def initialize(self, d_idx, d_mu, d_v00, d_v01, d_v10, d_v11,
                   d_g00, d_g01, d_g10, d_g11, d_gmag):
        # computes the shear rate tensor
        d_g00[d_idx] = d_v00[d_idx]
        d_g01[d_idx] = 0.5*(d_v01[d_idx] + d_v10[d_idx])
        d_g10[d_idx] = 0.5*(d_v10[d_idx] + d_v01[d_idx])
        d_g11[d_idx] = d_v11[d_idx]
        
        IIe = d_g00[d_idx]**2.0 + d_g01[d_idx]**2.0 + d_g10[d_idx]**2.0 + d_g11[d_idx]**2.0
        d_gmag[d_idx] = ((2.0*IIe)**0.5)
        d_mu[d_idx] = ((self.nu0 + self.nu1*self.kc*d_gmag[d_idx])/
                       (1.0+self.kc*d_gmag[d_idx]))*self.rho0

class ViscositySnowMix(Equation):
    def __init__(self, dest, sources, rho0=500.0, mp=1000.0,
                 nu0=2.1,nu1=2.7e-3,kc=1.1,
                 xc=2.5e3,lx=10.0):
        self.rho0 = rho0
        self.mp = mp
        self.nu0 = nu0
        self.nu1=nu1
        self.kc=kc
        self.xc=xc
        self.lx=lx
        super(ViscositySnowMix, self).__init__(dest, sources)
    
    def initialize(self, d_idx, d_mu, d_v00, d_v01, d_v10, d_v11,
                   d_g00, d_g01, d_g10, d_g11, d_gmag,d_x):
        # computes the shear rate tensor
        d_g00[d_idx] = d_v00[d_idx]
        d_g01[d_idx] = 0.5*(d_v01[d_idx] + d_v10[d_idx])
        d_g10[d_idx] = 0.5*(d_v10[d_idx] + d_v01[d_idx])
        d_g11[d_idx] = d_v11[d_idx]
        
        IIe = d_g00[d_idx]**2.0 + d_g01[d_idx]**2.0 + d_g10[d_idx]**2.0 + d_g11[d_idx]**2.0
        a = 2.88e-4
        b = 1.42e-2
        c = 36.72e-3
        d = 4.92e-3
        mu_0 = a*exp(b*self.rho0)
        tau_0 = c*exp(d*self.rho0)
        cond1 = (d_x[d_idx] < (self.xc-self.lx/2.0))
        cond2 = (d_x[d_idx] > (self.xc+self.lx/2.0))
        if not(cond1 + cond2):
            # Fluidized
            if IIe>1e-8:
                d_gmag[d_idx] = ((2.0*IIe)**0.5)
                d_mu[d_idx] = mu_0 + tau_0*(1.0 - exp(-self.mp*d_gmag[d_idx]))/d_gmag[d_idx]
            else:
                d_gmag[d_idx] = 0.0
                d_mu[d_idx] = mu_0 + tau_0*self.mp
        else:
            # Not Fluidized
            d_gmag[d_idx] = ((2.0*IIe)**0.5)
            d_mu[d_idx] = ((self.nu0 + self.nu1*self.kc*d_gmag[d_idx])/
                        (1.0+self.kc*d_gmag[d_idx]))*self.rho0

class TensorsTauEps2D(Equation):
    """Computes Epsilon and Tau tensors, and computes the momentum out of it
    """
    def __init__(self, dest, sources, mu0=1e-3, mp=100.0, t1=0.0, ty=0.0):
        self.mu0 = mu0
        self.ty = ty
        self.t1 = t1
        self.mp = mp
        super(TensorsTauEps2D, self).__init__(dest, sources)

    def initialize(self, d_idx,
                   d_v00, d_v01, d_v10, d_v11,
                   d_g00, d_g01, d_g10, d_g11, d_gmag,
                   d_t00, d_t01, d_t10, d_t11, d_tmag,
                   d_mu):
        
        # computes the shear rate tensor
        d_g00[d_idx] = d_v00[d_idx]
        d_g01[d_idx] = 0.5*(d_v01[d_idx] + d_v10[d_idx])
        d_g10[d_idx] = 0.5*(d_v10[d_idx] + d_v01[d_idx])
        d_g11[d_idx] = d_v11[d_idx]
        
        IIe = d_g00[d_idx]**2.0 + d_g01[d_idx]**2.0 + d_g10[d_idx]**2.0 + d_g11[d_idx]**2.0
        
        IIe = d_g00[d_idx]**2+2*d_g01[d_idx]**2.0+d_g11[d_idx]**2.0
        # |g|
        muo = self.mu0
        mp = self.mp
        ty = self.ty
        t1 = self.t1
        # computes the shear rate magnitude (second invariant)
        # Gamma -> 0 => mu_app -> mu0 + mp*tau0
        if IIe>1e-8:
            d_gmag[d_idx] = ((2.0*IIe)**0.5)
            d_mu[d_idx]   = (
                muo*exp(-t1*d_gmag[d_idx]) +
                ty*(1-exp(-mp*d_gmag[d_idx]))/d_gmag[d_idx]
            )
        else:
            d_gmag[d_idx] = 0.0
            d_mu[d_idx] = muo+mp*ty

class ViscousNNewtNaCl2D(Equation):
    def __init__(self, dest, sources, mu0=1e-3, mp=1000.0):
        """Computes Bingham Viscosity.
        <!> Supposes a Costa Viscosity is applied first!
        """
        self.mu0 = mu0
        self.mp = mp
        super(ViscousNNewtNaCl2D, self).__init__(dest, sources)

    def initialize(self, d_idx,
                   d_v00, d_v01, d_v10, d_v11,
                   d_g00, d_g01, d_g10, d_g11, d_gmag,
                   d_mu, d_mu_c, d_qstate):
        
        # computes the shear rate tensor
        d_g00[d_idx] = d_v00[d_idx]
        d_g01[d_idx] = 0.5*(d_v01[d_idx] + d_v10[d_idx])
        d_g10[d_idx] = 0.5*(d_v10[d_idx] + d_v01[d_idx])
        d_g11[d_idx] = d_v11[d_idx]
        
        IIe = d_g00[d_idx]**2.0 + d_g01[d_idx]**2.0 + d_g10[d_idx]**2.0 + d_g11[d_idx]**2.0
        # |g|
        mp = self.mp
        ty = 0.1341*exp(10.14*d_qstate[d_idx])
        # computes the shear rate magnitude (second invariant)
        # Gamma -> 0 => mu_app -> mu0 + mp*tau0
        if IIe>1e-8:
            d_gmag[d_idx] = ((2.0*IIe)**0.5)
            d_mu[d_idx]   = (
                self.mu0*d_mu_c[d_idx] +
                ty*(1.0 - exp(-mp*d_gmag[d_idx]))/d_gmag[d_idx]
            ) 
        else:
            d_gmag[d_idx] = 0.0
            d_mu[d_idx] = (self.mu0*d_mu_c[d_idx]+mp*ty)

class MomentumMorris(Equation):
    r"""Computes the Momentum equation by using a dynamic viscosity
    with the method proposed by Morris et al.
    WARNING: Default to Laminar flow! If turbulence, use `ContinuityEquationLANS`,
    `LANScorrection` and `XSPHCorrectionForLeapFrog` with `eps=0.8`.
    
    Ref:
        Morris, J. P., Fox, P. J., & Zhu, Y. (1997). Modeling low Reynolds number
        incompressible flows using SPH. Journal of computational physics,
        136(1), 214-226.

    Args:
        gx (float, optional): Acceleration due to gravity /x. Default to 0.0.
        gy (float, optional): Acceleration due to gravity /y. Default to -9.81.
    """
    def __init__(self, dest, sources, gx=0.0, gy=-9.81):
        self.gx = gx
        self.gy = gy
        super(MomentumMorris, self).__init__(dest, sources)
    
    def initialize(self, d_au, d_av, d_idx):
        d_au[d_idx] = self.gx
        d_av[d_idx] = self.gy
    
    def loop(self, d_au, d_av, d_p, d_rho, d_idx, s_m, s_p, s_rho, s_idx, d_mu, s_mu,
             VIJ, XIJ, R2IJ, EPS, DWIJ):
        Pressure = (
            (d_p[d_idx] + s_p[s_idx])/(d_rho[d_idx]*s_rho[s_idx])
        )
        # to remove if doesn't work
        vijdotrij = VIJ[0]*XIJ[0] + VIJ[1]*XIJ[1]
        mu_ij = 0.0
        if vijdotrij <0.0:
            mu_ij = 4*(d_mu[d_idx]*s_mu[s_idx])/(d_mu[d_idx]+s_mu[s_idx])
        V_jr  = s_m[s_idx]/(d_rho[d_idx]*s_rho[s_idx])
        dot   = XIJ[0]*DWIJ[0] + XIJ[1]*DWIJ[1]
        F_ij  = dot/(R2IJ + EPS)
        
        Visc = V_jr*mu_ij*F_ij
        
        d_au[d_idx] += -s_m[s_idx] * Pressure * DWIJ[0] + Visc*VIJ[0]
        d_av[d_idx] += -s_m[s_idx] * Pressure * DWIJ[1] + Visc*VIJ[1]

class MomentumAdami(Equation):
    r"""Computes the Momentum equation by using a dynamic viscosity
    with the transport velocity method proposed by Adami et al.
    
    Args:
        gx (float, optional): Acceleration due to gravity /x. Default to 0.0.
        gy (float, optional): Acceleration due to gravity /y. Default to -9.81.
    Ref:
        Adami, S., Hu, X. Y., & Adams, N. A. (2013).
            A transport-velocity formulation for smoothed particle hydrodynamics.
            Journal of Computational Physics, 241, 292-307.
    """
    def __init__(self, dest, sources, pb=0.0, gx=0.0, gy=-9.81):
        self.gx = gx
        self.gy = gy
        self.pb = pb
        super(MomentumAdami, self).__init__(dest, sources)
    
    def initialize(self, d_au, d_av, d_auhat, d_avhat, d_idx):
        d_au[d_idx] = self.gx
        d_av[d_idx] = self.gy
        d_auhat[d_idx] = 0.0
        d_avhat[d_idx] = 0.0
        
    
    def loop(self, d_au, d_av, d_p, d_rho, d_idx, d_m, s_m, s_p, s_rho, s_idx,
             d_mu, s_mu, d_u, d_v, d_uhat, d_vhat,s_u, s_v, s_uhat, s_vhat,
             d_auhat, d_avhat,
             VIJ, XIJ, R2IJ, EPS, DWIJ):
        rhoi = d_rho[d_idx]
        ui = d_u[d_idx]
        uhati = d_uhat[d_idx]
        vi = d_v[d_idx]
        vhati = d_vhat[d_idx]
        
        rhoj = s_rho[s_idx]
        uj = s_u[s_idx]
        uhatj = s_uhat[s_idx]
        vj = s_v[s_idx]
        vhatj = s_vhat[s_idx]
        
        
        Vi = d_m[d_idx]/d_rho[d_idx]
        Vj = s_m[s_idx]/s_rho[s_idx]
        
        Vi2 = Vi*Vi
        Vj2 = Vj*Vj
        
        
        Axxi = rhoi*ui*(uhati - ui)
        Axyi = rhoi*ui*(vhati - vi)
        Ayxi = rhoi*vi*(uhati - ui)
        Ayyi = rhoi*vi*(vhati - vi)
        
        Axxj = rhoj*uj*(uhatj - uj)
        Axyj = rhoj*uj*(vhatj - vj)
        Ayxj = rhoj*vj*(uhatj - uj)
        Ayyj = rhoj*vj*(vhatj - vj)
        
        Ax = 0.5*(
            (Axxi + Axxj)*DWIJ[0] +
            (Axyi + Axyj)*DWIJ[1]
        )
        Ay = 0.5*(
            (Ayxi + Ayxj)*DWIJ[0] +
            (Ayyi + Ayyj)*DWIJ[1]
        )
        
        pij = (
            (d_p[d_idx]*s_rho[s_idx] + s_p[s_idx]*d_rho[d_idx])/
            (d_rho[d_idx] + s_rho[s_idx])
        )
        mu_ij = 2.0*d_mu[d_idx]*s_mu[s_idx]/(d_mu[d_idx] + s_mu[s_idx])
        mi1 = 1.0/d_m[d_idx]
        
        Fij = DWIJ[0]*XIJ[0] + DWIJ[1]*XIJ[1]
        
        
        
        Pressure = -mi1*(Vi*Vi + Vj*Vj)*pij
        Visc = mi1*(Vi*Vi + Vj*Vj)*mu_ij*Fij/(R2IJ + EPS)
        AdvecX = mi1*(Vi2 + Vj2) * Ax
        AdvecY = mi1*(Vi2 + Vj2) * Ay
        
        d_au[d_idx] += Pressure * DWIJ[0] + Visc*VIJ[0] + AdvecX
        d_av[d_idx] += Pressure * DWIJ[1] + Visc*VIJ[1] + AdvecY
        d_auhat[d_idx] += - self.pb * mi1 * (Vi*Vi + Vj*Vj) * DWIJ[0]
        d_avhat[d_idx] += - self.pb * mi1 * (Vi*Vi + Vj*Vj) * DWIJ[1]

class MomentumMonaghanDynamic(Equation):
    """Computes the Momentum equation by using a dynamic viscosity
    with the method proposed by Cleary.
    Ref:
        Cleary, P. W. (1998). Modelling confined multi-material heat and mass
            flows using SPH. Applied Mathematical Modelling, 22(12), 981-993.

    Args:
        gx (float, optional): Acceleration due to gravity /x. Default to 0.0.
        gy (float, optional): Acceleration due to gravity /y. Default to -9.81.
    """
    def __init__(self, dest, sources, gx=0.0, gy=-9.81):
        self.gx = gx
        self.gy = gy
        super(MomentumMonaghanDynamic, self).__init__(dest, sources)
    
    def initialize(self, d_au, d_av, d_idx):
        d_au[d_idx] = self.gx
        d_av[d_idx] = self.gy
    
    def loop(self, d_au, d_av, d_p, d_rho, d_idx, s_m, s_p, s_rho, s_idx, d_mu, s_mu,
             VIJ, XIJ, R2IJ, EPS, DWIJ):
        Pressure = (
            (d_p[d_idx] + s_p[s_idx])/(d_rho[d_idx]*s_rho[s_idx])
        )
        # to remove if doesn't work
        vijdotrij = VIJ[0]*XIJ[0] + VIJ[1]*XIJ[1]
        pi_ij = 0.0
        if vijdotrij <0.0:
            mu_ij = (
                (4.0 * d_mu[d_idx] * s_mu[s_idx])/
                (d_mu[d_idx] + s_mu[s_idx])
            )
            dot = vijdotrij/(R2IJ+EPS)
            pi_ij = -4.0 * mu_ij/(d_rho[d_idx]*s_rho[s_idx])*dot
        
        d_au[d_idx] += -s_m[s_idx] * (Pressure+pi_ij) * DWIJ[0]
        d_av[d_idx] += -s_m[s_idx] * (Pressure+pi_ij) * DWIJ[1]


class ViscosityCosta(Equation):
    r"""Computes the dynamic viscosity related to the ice fraction using a
    Costa model [Mader et al, 2013]

    Args:
        dest (str): Destination particles
        sources (list): List of sources particles
        mu0 (float, optional): Dynamic viscosity of the fully liquid
        flow.
        Defaults to 1e-3.
        B (float, optional): Einstein coefficient. Defaults to 2.5.
        beta (list, optional): List = (\phi_*, alpha, delta, chi).
        Defaults to (0.6, 5.0, 8.0, 1.0e-4).
    """
    def __init__(self, dest, sources, mu0:float=1e-3, B:float=2.5,
                 phi_s:float=0.6, alpha:float=5.0, delta:float=8.0,
                 chi:float=1.0e-4):
        
        self.mu0=mu0
        self.B = B
        self.phi_s = phi_s
        self.alpha = alpha
        self.delta = delta
        self.chi = chi
        super(ViscosityCosta, self).__init__(dest, sources)
        
    def initialize(self, d_mu, d_qstate, d_idx):
        phi = d_qstate[d_idx]
        B = self.B
        phi_s = self.phi_s
        alpha = self.alpha
        delta = self.delta
        chi = self.chi
        chi1 = 1.0 - chi
        phi_r = phi/phi_s
        F = chi1 * erf(
            sqrt(pi)/(2.0*chi1)*phi_r*(1.0+pow(phi_r, alpha))
        )
        mu_r = (1.0 + pow(phi_r, delta))/pow(1.0 - F, B*phi_s)
        d_mu[d_idx] = self.mu0 * mu_r
        
class ViscosityCostaNNewt(Equation):
    r"""Computes the relative viscosity to the ice fraction using a
    Costa model [Mader et al, 2013]
    
    <!> Used for Non Newtonian computations only.

    Args:
        dest (str): Destination particles
        sources (list): List of sources particles
        mu0 (float, optional): Dynamic viscosity of the fully liquid
        flow.
        Defaults to 1e-3.
        B (float, optional): Einstein coefficient. Defaults to 2.5.
        beta (list, optional): List = (\phi_*, alpha, delta, chi).
        Defaults to (0.6, 5.0, 8.0, 1.0e-4).
    """
    def __init__(self, dest, sources, mu0:float=1e-3, B:float=2.5,
                 phi_s:float=0.6, alpha:float=5.0, delta:float=8.0,
                 chi:float=1.0e-4):
        self.B = B
        self.phi_s = phi_s
        self.alpha = alpha
        self.delta = delta
        self.chi = chi
        super(ViscosityCostaNNewt, self).__init__(dest, sources)
        
    def initialize(self, d_mu_c, d_qstate, d_idx):
        phi = d_qstate[d_idx]
        B = self.B
        phi_s = self.phi_s
        alpha = self.alpha
        delta = self.delta
        chi = self.chi
        chi1 = 1.0 - chi
        phi_r = phi/phi_s
        F = chi1 * erf(
            sqrt(pi)/(2.0*chi1)*phi_r*(1.0+pow(phi_r, alpha))
        )
        mu_r = (1.0 + pow(phi_r, delta))/pow(1.0 - F, B*phi_s)
        d_mu_c[d_idx] = mu_r

class ViscosityMean(Equation):
    """Computes an averaged viscosity from the neighbouring
    particles
    """
    def loop_all(self, d_idx, d_mu, s_mu, NBRS, N_NBRS):
        s_idx = declare('long')
        my_sum = 0.0
        if N_NBRS != 0:
            for i in range(N_NBRS):
                s_idx = NBRS[i]
                my_sum += s_mu[s_idx]
            d_mu[d_idx] = my_sum/N_NBRS
        
        
class XSPHCorrectionForLeapFrog(Equation):
    """Computes the smoothed velocity term for
    XSPH correction.
    Ref:
        Monaghan, J.J.: On the Problem of Penetration
            in Particle Methods, Journal of Computational
            Physics 82, 1-15 (1989)

    Args:
        eps (float, optional): Epsilon factor for the XSPH correction.
        Default to 0.5.
    """
    def __init__(self, dest, sources, eps=0.5):
        self.eps = eps
        super(XSPHCorrectionForLeapFrog, self).__init__(dest, sources)
    
    def initialize(self, d_idx, d_usph, d_vsph):
        d_usph[d_idx] = 0.0
        d_vsph[d_idx] = 0.0

    def loop(self, s_idx, d_idx, s_m, d_usph, d_vsph, WIJ, s_rho, VIJ):
        tmp = -self.eps * s_m[s_idx]*WIJ/s_rho[s_idx]

        d_usph[d_idx] += tmp * VIJ[0]
        d_vsph[d_idx] += tmp * VIJ[1]
        
class LANScorrection(Equation):
    def __init__(self, dest, sources, eps=0.5, rho0=1000.0):
        self.eps=eps
        self.rho0 = rho0
        super(LANScorrection, self).__init__(dest, sources)
    
    def loop(self, d_au, d_av, s_idx, s_m,d_idx,RHOIJ1, VIJ, DWIJ):
        vij2 = VIJ[0]*VIJ[0] + VIJ[1]*VIJ[1]
        fact = s_m[s_idx]*0.5*self.eps*vij2*RHOIJ1
        
        d_au[d_idx] += fact*DWIJ[0]
        d_av[d_idx] += fact*DWIJ[1]
        
class InletVelocity(Equation):
    def __init__(self, dest, sources, vmax:float=0.1, time:float=20.0, Xc:float=20.0, R:float=5.0):
        self.tmax = time
        self.vmax = vmax
        self.Xc = Xc
        self.R = R
        super(InletVelocity, self).__init__(dest, sources)
    
    def initialize(self, d_idx, d_v, d_x, t):
        vmax = self.vmax
        R = self.R
        if t <self.tmax:
            d_v[d_idx] = vmax*(1.0-((d_x[d_idx]-self.Xc)/R)**2.0)
        else:
            d_v[d_idx] = 0.0