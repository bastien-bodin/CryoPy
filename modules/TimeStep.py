from pysph.sph.equation import Equation

class AdaptiveTSMorris(Equation):
    """Computes the adaptive time step for fluid dynamics using the dynamic viscosity
    presented by Morris et al.
    
    Ref:
        Morris, J. P., Fox, P. J., & Zhu, Y. (1997). Modeling low Reynolds number
        incompressible flows using SPH. Journal of computational physics, 136(1), 214-226.
    """
    def initialize(self, d_idx, d_h, d_cs, d_au, d_av, d_dt_adapt, d_dtf, d_dtc,
                   d_dtv, d_rho, d_mu):
        C1, C2, C3 = 0.3, 0.3, 0.125
        f = (d_au[d_idx]**2.0 + d_av[d_idx]**2.0)**0.5
        nu = d_mu[d_idx]/d_rho[d_idx]
        
        d_dtf[d_idx] = C1 * (d_h[d_idx]/f)**0.5
        d_dtc[d_idx] = C2 * d_h[d_idx] / d_cs[d_idx]
        d_dtv[d_idx] = C3 * d_h[d_idx]**2.0/nu
        d_dt_adapt[d_idx] = min(d_dtf[d_idx], d_dtc[d_idx], d_dtv[d_idx])

class AdaptiveTSThermo(Equation):
    """Computes the adaptive time step for thermodynamics using the method of
    Cleary and Monaghan.
    
    Ref:
        Cleary, P. W., & Monaghan, J. J. (1999). Conduction modelling using smoothed
        particle hydrodynamics. Journal of Computational Physics, 148(1), 227-264.
    """
    def initialize(self, d_idx, d_h, d_rho, d_cp, d_kappa, d_dttherm, d_dt_adapt):
        C4 = 0.1
        
        alpha = d_kappa[d_idx]/(d_rho[d_idx]*d_cp[d_idx])
        d_dttherm[d_idx] = C4 * d_h**2.0/alpha
        d_dt_adapt[d_idx] = d_dttherm[d_idx]

class AdaptiveTSCoupled(Equation):
    """Computes the adaptive time step for BOTH fluid dynamics AND thermodynamics
    using the dynamic viscosity presented by Morris et al and Cleary and Monaghan.

    Refs:
        Morris, J. P., Fox, P. J., & Zhu, Y. (1997). Modeling low Reynolds number
        incompressible flows using SPH. Journal of computational physics, 136(1), 214-226.
        
        Cleary, P. W., & Monaghan, J. J. (1999). Conduction modelling using smoothed
        particle hydrodynamics. Journal of Computational Physics, 148(1), 227-264.
    """
    def initialize(self, d_idx, d_h, d_cs, d_au, d_av, d_dt_adapt, d_dtf, d_dtc,
                   d_dtv, d_rho, d_mu, d_cp, d_kappa, d_dttherm):
        C1, C2, C3, C4 = 0.3, 0.3, 0.125, 0.1
        f = (d_au[d_idx]**2.0 + d_av[d_idx]**2.0)**0.5
        nu = d_mu[d_idx]/d_rho[d_idx]

        d_dtf[d_idx] = C1 * (d_h[d_idx]/f)**0.5
        d_dtc[d_idx] = C2 * d_h[d_idx] / d_cs[d_idx]
        d_dtv[d_idx] = C3 * d_h[d_idx]**2.0/nu
        alpha = d_kappa[d_idx]/(d_rho[d_idx]*d_cp[d_idx])
        d_dttherm[d_idx] = C4 * d_h[d_idx]**2.0/alpha
        d_dt_adapt[d_idx] = min(
            d_dtf[d_idx], d_dtc[d_idx], d_dtv[d_idx], d_dttherm[d_idx]
        )
