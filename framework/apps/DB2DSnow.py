# Personal modules derived from PySPH Library
from modules.Geometries import MyFissure2D, MyBlock2D, MyPit2D, MyCuve2D
from modules.FluidDynamics import (
    TaitEoS, ContinuityEquation, ContinuityEquationLANS,
    GradientCorrection, MomentumMorris, ViscosityCosta,
    ViscosityMean, XSPHCorrectionForLeapFrog, LANScorrection,
    InletVelocity, ViscosityCross, VelGrad2D
)
from modules.HeatTransfer import (
    CenterOfMass, UpdateThermalProps, MassCluster,
    Conduction, Radiation, HeatTransfer, Convection
)
from modules.Integrators import (
    VerletIntegrator,
    VerletStepFluid, VerletStepDBC, VerletStepInlet,
    VerletStepCoupledSolidificationAndMeltPure
)
from modules.TimeStep import AdaptiveTSMorris, AdaptiveTSCoupled

# PySPH imports
from pysph.sph.equation import Group
from pysph.base.utils import get_particle_array_wcsph
from pysph.solver.application import Application
from pysph.solver.solver import Solver
from pysph.base.kernels import CubicSpline

# PySPH modules for Inlet
from pysph.sph.bc.donothing.simple_inlet_outlet import (
    SimpleInletOutlet
)
from pysph.sph.bc.inlet_outlet_manager import InletInfo

class DB2D(Application):
    def __init__(self, fname=None, output_dir=None, domain=None,
                 x0:float=0.0, y0:float=0.0,
                 Lx:float=5.0e3, Ly:float=100.0,
                 lx_fluid:float=80.0, ly_fluid:float=50.0,
                 dp:float=0.05, vmax=1.0,
                 qstate:float=0.0, rho0:float=1000.0, mu0:float=1.0e-3,
                 cp_l:float=4187.0, cp_s:float=2050.0,
                 k_l:float=0.598, k_s:float=2.22,
                 Lh:float=3.34e5, Tsol:float=273.15, Tenv:float=93.7,
                 tf:float=10.0, t_inlet:float=None):
        """SPH Application for a Dam Break where both Hydrodynamics
           and Thermodynamics are at play

        Args:
            fname (str, optional): Name of the simulation. Defaults to None.
            output_dir (str, optional): Path of the directory where the outputs will
            be placed. Defaults to None.
            domain (`pysph.base.nnps_base.DomainManager`, optional): _description_. Defaults to None.
            x0 (float, optional): Starting position /x. Defaults to 0.0.
            y0 (float, optional): Starting position /y. Defaults to 0.0.
            Lx (float, optional): Length of the cuve. Defaults to 4.0.
            Ly (float, optional): Height of the cuve. Defaults to 4.0.
            lx_fluid (float, optional): Length of fluid. Defaults to 1.0.
            ly_fluid (float, optional): Height of fluid. Defaults to 2.0.
            dp (float, optional): Spacing between particles. Defaults to 0.05.
            qstate (float, optional): Initial ice fraction. Defaults to 0.0.
            rho0 (float, optional): Initial density of the fluid. Defaults to 1000.0.
            mu0 (float, optional): Initial dynamic viscosity of the fluid. Defaults to 1.0e-3.
            cp_l (float, optional): Specific heat capacity of the liquid fluid. Defaults to 4187.0.
            cp_s (float, optional): Specific heat capacity of the solidified fluid. Defaults to 2050.0.
            k_l (float, optional): Thermal conductivity of the fluid. Defaults to 0.598.
            k_s (float, optional): Thermal conductivity of the solidified fluid. Defaults to 2.22.
            Lh (float, optional): Enthalpy of fusion. Defaults to 3.34e5.
            Tsol (float, optional): Solidification temperature. Defaults to 273.15.
            Tenv (float, optional): Environment temperature. Defaults to 105.0.
            tf (float, optional): Final computed time. Defaults to 2.0.
        """
        self.dp = dp; self.vmax = vmax
        self.x0, self.y0 = x0, y0
        self.Lx, self.Ly = Lx, Ly
        self.rho0 = rho0
        self.mu0 = mu0
        self.qstate=qstate
        self.tf = tf
        if t_inlet == None:
            self.t_inlet = self.tf
        else:
            self.t_inlet = t_inlet
        # fluid geometry
        self.lx, self.ly = lx_fluid, ly_fluid
        super(DB2D, self).__init__(fname, output_dir, domain)
        
    def initialize(self):
        self.hdp = 1.5
        self.h = self.dp*self.hdp
        self.gx, self.gy = 0.0, -1.36
        self.eps = 0.8
        
        self.mass = self.rho0 * self.dp**2.0
        self.c0 = 10*(self.Ly*(-self.gy))**0.5
        
        # Parameters for Costa rheological model
        self.B = 2.5
        self.phi_s = 0.6
        self.alpha = 5.0
        self.delta = 8.0
        self.chi = 1.0e-4
        
        Npoints = 2000
        self.output_times = [i*self.tf/Npoints for i in range(0,Npoints+1)]
        
    def create_particles(self):
        Xb, Yb = MyCuve2D(
            self.x0, self.y0, self.Lx, self.Ly, self.dp, nl=2
        )
        Xf, Yf = MyBlock2D(
            x0=self.x0+0.5*(self.Lx-self.lx), y0=self.y0,
            Lx=self.lx, Ly=self.ly, dp=self.dp
        )
        
        
        boundaries = get_particle_array_wcsph(
            name='bound', x=Xb, y=Yb, h=self.h, m=self.mass,
            rho=self.rho0, cs=self.c0
        )
        fluid = get_particle_array_wcsph(
            name='fluid', x=Xf, y=Yf, h=self.h, m=self.mass,
            rho=self.rho0, cs=self.c0
        )
        
        PA = [fluid, boundaries]
        
        PropsMu = ['mu']
        PropsGrad = ['v00', 'v01', 'v10', 'v11',
                     'g00', 'g01', 'g10', 'g11', 'gmag',
                     't00', 't01', 't10', 't11', 'tmag']
        PropsXSPH = ['usph','vsph']
        PropsDt = ['dt_adapt', 'dtf', 'dtc', 'dtv', 'dttherm']
        PropsVerlet = ['rho0']
        # PropsCoM = ['Surface', 'x_com', 'y_com', 'MassCluster']
        # PropsSol = ['qstate', 'Temp', 'Temp0', 'kappa', 'cp',
        #             'Cond', 'Conv', 'Rad', 'aTemp', 'Lh', 'Tsol']
        PropsIO = ['ioid', 'disp']
        PropsGradCorr = ['m_mat']
        
        Props = (PropsMu + PropsGrad +
                 PropsDt + PropsVerlet +
                 # PropsSol + PropsCoM +
                 PropsXSPH + PropsIO)
        
        for pa in PA:
            for prop in Props:
                pa.add_property(prop, data=0.0)
            for prop in PropsGradCorr:
                pa.add_property(prop, data=0.0, stride=3)
        
        for pa in PA:
            pa.mu[:] = self.mu0
            # pa.qstate[:] = self.qstate
            pa.dt_adapt[:] = 1e-7
            # pa.Lh[:] = self.Lh
            # pa.Tsol[:] = self.Tsol
        PA[1].dt_adapt[:] = 1e7
        
        # inlet
        # PA[2].Temp[:] = self.Tsol
        # PA[2].kappa[:] = self.k_l
        # PA[2].cp[:]=self.cp_l
        
        # fluid
        # PA[0].Temp[:] = self.Tsol
        # PA[0].kappa[:] = self.k_l
        # PA[0].cp[:]=self.cp_l
        
        # boundaries
        # PA[1].Temp[:] = self.Tenv
        # PA[1].kappa[:] = self.k_s
        # PA[1].cp[:]=self.cp_s
        # PA[1].qstate[:] = 1.0
        
        
        
        for pa in PA:
            pa.add_output_arrays(Props)
        
        return PA
    
    def create_equations(self):
        xc = self.x0 + 0.5*self.Lx
        equations = [
            #Group(equations=[
            #    GradientCorrection(
            #        dest='fluid', sources=['fluid', 'bound', 'inlet']
            #    ),
            #    GradientCorrection(
            #        dest='bound', sources=['fluid', 'bound', 'inlet']
            #    ),
            #    GradientCorrection(
            #        dest='inlet', sources=['fluid', 'bound', 'inlet']
            #    ),
            #]),
            # Group(equations=[
            #     InletVelocity(
            #         dest='inlet', sources=None,
            #         Xc=xc, R = self.lx/2.0,
            #         time=self.t_inlet, vmax=self.vmax
            #     ),
            # ]),
            Group(equations=[
                VelGrad2D(dest='fluid', sources=['fluid', 'bound']),
                VelGrad2D(dest='bound', sources=['fluid']),
            ]),
            Group(equations=[
                ViscosityCross(dest='fluid', sources=None, rho0=self.rho0),
            ]),
            Group(equations=[
                ViscosityMean(dest='bound', sources=['fluid']),
            ]),
            Group(equations=[
                TaitEoS(dest='fluid',sources=None, rho0=self.rho0,
                        c0=self.c0),
                TaitEoS(dest='bound',sources=None, rho0=self.rho0,
                        c0=self.c0),
            ]),
            Group(equations=[
                ContinuityEquationLANS(
                    dest='fluid',sources=['fluid','bound'],
                ),
                ContinuityEquation(
                    dest='bound',sources=['fluid'],
                ),
                MomentumMorris(
                    dest='fluid', sources=['fluid', 'bound'],
                    gx = self.gx, gy=self.gy
                ),
            ]),
            Group(equations=[
                LANScorrection(
                    dest='fluid', sources=['fluid'],
                    eps=self.eps
                ),
                XSPHCorrectionForLeapFrog(
                    dest='fluid', sources=['fluid'],
                    eps=self.eps
                )
            ]),
            Group(equations=[
                AdaptiveTSMorris(dest='fluid', sources=None),
            ])
        ]
        return equations
    
    # def _create_inlet_outlet_manager(self):
    #     from pysph.sph.bc.donothing.inlet import Inlet
    #     props_to_copy = ['x', 'y', 'z', 'u', 'v', 'w', 'm', 'h', 'rho', 'p', 'ioid']
    #     x0pit = self.x0 + 0.5*(self.Lx-self.lx)
    #     inlet_info = InletInfo(
    #         pa_name='inlet',
    #         normal=[0.0, -1.0, 0.0],
    #         refpoint=[x0pit,
    #                   0.0, # self.y0-self.ly+0.5*self.dp,
    #                   0.0],
    #         has_ghost=False,
    #         update_cls=Inlet
    #     )
        
    #     iom = SimpleInletOutlet(
    #         fluid_arrays=['fluid'], inletinfo=[inlet_info],
    #         outletinfo=[]
    #     )
    #     return iom
    
    # def create_inlet_outlet(self, particle_arrays):
    #     iom = self.iom
    #     io = iom.get_inlet_outlet(particle_arrays)
    #     return io
    
    def create_solver(self):
        # self.iom = self._create_inlet_outlet_manager()
        kernel = CubicSpline(dim=2)
        integrator = VerletIntegrator(
            fluid = VerletStepFluid(self.rho0),
            bound = VerletStepDBC(self.rho0),
            # inlet = VerletStepInlet(self.rho0)
        )
        tf = self.tf
        # self.iom.active_stages = [2]
        # self.iom.setup_iom(dim=2, kernel=kernel)
        # self.iom.update_dx(dx=self.dp)
        solver = Solver(dim=2, integrator=integrator, kernel=kernel,
                        tf = tf, adaptive_timestep=True, fixed_h=True,
                        output_at_times=self.output_times, pfreq=1e14)
        return solver
