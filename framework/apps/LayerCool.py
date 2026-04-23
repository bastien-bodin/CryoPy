from modules.Geometries import MyBlock2D

from modules.HeatTransfer import Conduction, Radiation, HeatTransfer
from modules.Integrators import (
    VerletIntegrator, VerletStepThermal
)
from modules.TimeStep import AdaptiveTSThermo
from math import exp, erf, sqrt, pi

from pysph.sph.equation import Equation
from compyle.api import declare

# PySPH imports
from pysph.sph.equation import Group
from pysph.base.utils import get_particle_array_wcsph
from pysph.solver.application import Application
from pysph.solver.solver import Solver
from pysph.base.kernels import CubicSpline
from pysph.base.nnps import DomainManager

# dp = 0.5
# L = 10
# W = 5
# D = 3*W
# G = 5*dp

# x_snow, y_snow = MyBlock2D(x0=0.0, y0=0.0, Lx=L,Ly=W,dp=dp)
# mask = (y_snow > (W-dp))

# x_ice, y_ice = MyBlock2D(x0=0, y0=-D, Lx=L,Ly=D,dp=dp)
# x_ground, y_ground = MyBlock2D(x0=0, y0=-D-G, Lx=L,Ly=G,dp=dp)

# plt.figure()
# plt.plot(x_snow[mask],y_snow[mask],'.')
# plt.plot(x_snow[~mask],y_snow[~mask],'.')
# plt.plot(x_ice,y_ice,'.')
# plt.plot(x_ground,y_ground,'.')
# plt.show()
def kappa_ice_f(T:float=0.0) -> float:
    """Computes the thermal conductivity of ice from the temperature
    (in degC). From Ratcliffe (1962)

    Args:
        T (float, optional): Temperature of ice in degC. Defaults to 0.0.

    Returns:
        float: Thermal conductivity in W/(m K)
    """
    return 612.0/T

def cp_ice_f(T:float=100.0) -> float:
    return (0.183 + 0.689e-2*T)*1.0e3

def kappa_snow_f(T, rho) -> float:
    """Computes the thermal conductivity of snow from porous model
    The porosity is defined as: P = 1-rho_snow/rho_ice

    The thermal conductivity of ice is deduced from the one of
    Wolfenbarger et al (2021)

    Reference: Smith et al (2013), Wolfenbarger et al (2021)
    """
    rho_ice = 917.0

    k_s = 612.0/ T

    C1 =  33.9729025; C2 = -164.702679; C3 =  262.108546
    C4 = -21.5346955; C5 = -443.455815; C6 =  607.339582
    C7 = -368.790121; C8 =  111.296674; C9 = -13.4122465
    Tc = 132.52
    Kappa = 4.358e-3
    Tr = T/Tc
    k_v = (
        C1/Tr   + C2 * Tr**(-2.0/3.0) + C3 * Tr**(-1.0/3.0) +
        C4      + C5 * Tr**(1.0/3.0)  + C6 * Tr**(2.0/3.0)  +
        C7 * Tr + C8 * Tr**(4.0/3.0)  + C9 *Tr**(5.0/3.0)
    )* Kappa

    P = 1.0 - rho/rho_ice
    if P < 0.15:
        # Maxwell-Eucken
        num = k_v + 2.0 * k_s + 2.0*P*(k_v - k_s)
        den = k_v + 2.0 * k_s -     P*(k_v - k_s)
        return k_s*num/den
    elif P < 0.65:
        # Landauer (1952)
        return 0.25*(
            k_v*(3*P-1) + k_s*(2-3*P) +
            (
                (k_v*(3*P-1) + k_s*(2-3*P))**2.0 +
                8*k_s*k_v
            )**0.5
        )
    else:
        # Russell (1935)
        return k_s*(
            k_s+P**(2.0/3.0)*(k_v-k_s)
        )/(
            k_s+(k_v-k_s)*(P**(2.0/3.0) - P)
        )

class kappa_ice(Equation):
    """Computes the thermal conductivity of ice
    Reference: Wolfenbarger et al (2021)
    """
    def initialize(self, d_idx, d_kappa, d_Temp):
        d_kappa[d_idx] = 612.0/d_Temp[d_idx]

class kappa_snow(Equation):
    """Computes the thermal conductivity of snow from porous model
    The porosity is defined as: P = 1-rho_snow/rho_ice

    The thermal conductivity of ice is deduced from the one of
    Wolfenbarger et al (2021)

    Reference: Smith et al (2013), Wolfenbarger et al (2021)
    """
    def initialize(self, d_idx, d_kappa, d_rho, d_Temp):
        rho_ice = 917.0

        k_s = 612.0/d_Temp[d_idx]

        C1 =  33.9729025; C2 = -164.702679; C3 =  262.108546
        C4 = -21.5346955; C5 = -443.455815; C6 =  607.339582
        C7 = -368.790121; C8 =  111.296674; C9 = -13.4122465
        Tc = 132.52
        Kappa = 4.358e-3
        Tr = d_Temp[d_idx]/Tc
        k_v = (
            C1/Tr   + C2 * Tr**(-2.0/3.0) + C3 * Tr**(-1.0/3.0) +
            C4      + C5 * Tr**(1.0/3.0)  + C6 * Tr**(2.0/3.0)  +
            C7 * Tr + C8 * Tr**(4.0/3.0)  + C9 *Tr**(5.0/3.0)
        )* Kappa

        P = 1.0 - d_rho[d_idx]/rho_ice
        if P < 0.15:
            # Maxwell-Eucken
            num = k_v + 2.0 * k_s + 2.0*P*(k_v - k_s)
            den = k_v + 2.0 * k_s -     P*(k_v - k_s)
            d_kappa[d_idx] = k_s*num/den
        elif P < 0.65:
            # Landauer (1952)
            d_kappa[d_idx] = 0.25*(
                k_v*(3*P-1) + k_s*(1-3*P) +
                (
                    (k_v*(3*P-1) + k_s*(1-3*P))**2.0 +
                    8*k_s*k_v
                )**0.5
            )
        else:
            # Russell (1935)
            d_kappa[d_idx] = k_s*(
                k_s+P**(2.0/3.0)*(k_v-k_s)
            )/(
                k_s+(k_v-k_s)*(P**(2.0/3.0) - P)
            )

class cp_ice(Equation):
    def initialize(self, d_idx, d_cp, d_Temp):
        d_cp[d_idx] = (0.183 + 0.689e-2*d_Temp[d_idx])*1.0e3

class ThermoSnowEuropa(Application):
    def __init__(self, fname=None, output_dir=None, domain=None,
                 dp:float=0.05, rho_snow:float=917.0,
                 W:float=1.0, tf_days:float=500.0):
        self.dp         : float = dp
        self.W          : float = W + dp
        self.rho_snow   : float = rho_snow
        self.tf         : float = tf_days*3600.0*24.0
        super(ThermoSnowEuropa, self).__init__(fname, output_dir, domain)

    def initialize(self):
        # self.dp         : float = 0.05
        self.hdp        : float = 1.3
        # self.rho_snow   : float = 917.0
        self.rho_ice    : float = 917.0
        self.T0         : float = 273.15
        self.Tf         : float = 96.0      # Ashkenazy (2019)
        self.gravity    : float = 1.315
        self.c0         : float = 10*(5.0*self.gravity)**0.5
        self.emissivity : float = 0.94      # Ashkenazy (2019)
        self.Tenv       : float = 96.0      # Ashkenazy (2019)


        h = self.dp*self.hdp
        # self.tf = 3600.0*24.0*20.0

        self.k_snow  : float = kappa_snow_f(self.T0, self.rho_snow)
        self.k_ice_0 : float = kappa_ice_f(self.Tf)
        self.cp_ice  : float = cp_ice_f(self.Tf)

        alpha = self.k_snow/(self.rho_snow*self.cp_ice)
        alpha_ice = self.k_ice_0/(self.rho_ice*self.cp_ice)



        # self.W = 1.0+self.dp # width of snow
        self.L = 1.0 # length of the semi-1d domain

        # depth of the ground (variable temperature)
        self.D = sqrt(10.0*alpha_ice*self.tf) # Based on the Fourier Number

        self.G = 5.0*self.dp # depth of the ground (fixed temperature)
        C4 = 0.1
        self.dt_snow : float = C4 * h*h/alpha
        self.dt_ice : float = C4 * h*h/alpha_ice

        Npoints = 200
        self.output_times = [i*self.tf/Npoints for i in range(0,Npoints+1)]

    def create_particles(self):
        L = self.L
        W = self.W
        D = self.D
        G = self.G

        m_snow = self.rho_snow*self.dp**2.0
        m_ice = self.rho_ice*self.dp**2.0
        h = self.dp*self.hdp

        x_snow, y_snow = MyBlock2D(
            x0=0.0, y0=0.0,
            Lx=L,Ly=W,
            dp=self.dp
        )
        mask = (y_snow > (W-self.dp))

        x_ice, y_ice = MyBlock2D(
            x0=0, y0=0.0,
            Lx=L,Ly=D,
            dp=self.dp
        )
        x_ground, y_ground = MyBlock2D(
            x0=0, y0=y_ice.max()+self.dp*0.5,
            Lx=L,Ly=G,
            dp=self.dp
        )

        # surf = get_particle_array_wcsph(
        #     name='surf', x=x_snow[mask], y=y_snow[mask], m=m_snow, h=h,
        #     rho=self.rho_snow, cs=self.c0
        # )
        snow = get_particle_array_wcsph(
            name='snow', x=x_snow, y=y_snow, m=m_snow, h=h,
            rho=self.rho_snow, cs=self.c0
        )
        ice = get_particle_array_wcsph(
            name='ice', x=x_ice, y=-y_ice, m=m_ice, h=h,
            rho=self.rho_ice, cs=self.c0
        )
        ground = get_particle_array_wcsph(
            name='ground', x=x_ground, y=-y_ground, m=m_ice, h=h,
            rho=self.rho_ice, cs=self.c0
        )

        PropsThermo = [
            'kappa', 'cp', 'Temp','Temp0', 'aTemp', 'Cond', 'Conv', 'Rad'
        ]
        PropsDt = ['dt_adapt', 'dtf', 'dtc', 'dtv', 'dttherm', 'Surface']

        Props = PropsThermo + PropsDt
        for pa in (
            # surf,
            snow,
            ice,
            ground
        ):
            for prop in Props:
                pa.add_property(prop, data=0.0)
            pa.cp[:] = self.cp_ice
            pa.dt_adapt[:] = 1.0e-7
            pa.Temp[:] = self.Tf

        snow.Surface[mask] = 1.0
        # surf.Temp[:] = self.T0
        # surf.kappa[:]   = self.k_snow

        snow.dt_adapt[:] = self.dt_snow
        snow.Temp[:] = self.T0
        snow.kappa[:]   = self.k_snow

        ice.kappa[:]    = self.k_ice_0
        ground.kappa[:] = self.k_ice_0
        ground.dt_adapt[:] = self.dt_ice

        for pa in (
            # surf,
            snow,
            ice,
            ground
        ):
            pa.add_output_arrays(Props)

        return [
            # surf,
            snow,
            ice,
            ground
        ]

    def create_equations(self):
        equations = [
            Group(equations=[
                kappa_ice(dest='ice', sources=None),
                kappa_snow(dest='snow', sources=None),
            #     kappa_snow(dest='surf', sources=None),
                cp_ice(dest='ice', sources=None),
                cp_ice(dest='snow', sources=None),
            #     cp_ice(dest='surf', sources=None),
            ]),
            Group(equations=[
                Conduction(
                    dest='snow', sources=['snow', 'ice']
                ),
                Conduction(
                    dest='ice', sources=['snow', 'ice', 'ground']
                ),
                Radiation(
                    dest='snow', sources=None, Tenv=self.Tenv,
                    eps=self.emissivity
                )
            ]),
            Group(equations=[
                HeatTransfer(dest='snow', sources=None, dp=self.dp),
                HeatTransfer(dest='ice', sources=None, dp=self.dp),
                # HeatTransfer(dest='surf', sources=None, dp=self.dp),
            ]),
            Group(equations=[
                AdaptiveTSThermo(dest='ice', sources=None),
                AdaptiveTSThermo(dest='snow', sources=None),
                # AdaptiveTSThermo(dest='surf', sources=None),
            ])
        ]
        return equations

    def create_domain(self):
        return DomainManager(
            xmin=0.0, xmax=self.L, ymin=-self.D-self.G, ymax=self.W,
            periodic_in_x=True, periodic_in_y=False
        )

    def create_solver(self):
        kernel = CubicSpline(dim=2)
        integrator = VerletIntegrator(
            snow = VerletStepThermal(),
            ice = VerletStepThermal(),
            # surf = VerletStepThermal(),
        )
        tf = self.tf
        solver = Solver(
            dim=2, integrator=integrator, kernel=kernel,
            tf = tf, adaptive_timestep=True, fixed_h=True,
            output_at_times=self.output_times, pfreq=1e14
        )
        return solver

if __name__=='__main__':
    app = ThermoSnowEuropa(fname='EuropaThermoRad_e-09_W-1m_rho_917-tf-20')
    app.run()
