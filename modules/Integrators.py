from pysph.sph.integrator import Integrator
from pysph.sph.integrator_step import IntegratorStep

"""
- Verlet Integrator
- Integrator Steps, Fluid Dynamics only
- Integrator Steps, Thermal aspect only
- Integrator Steps, Coupled Hydro/Thermo
"""
#######################################################################
# Verlet Integrator
#######################################################################
class VerletIntegrator(Integrator):
    def one_timestep(self, t, dt):
        self.initialize()
        self.compute_accelerations()
        self.stage1()
        self.update_domain()
        self.do_post_stage(0.5*dt,1)
        self.compute_accelerations()
        self.stage2()
        self.update_domain()
        self.do_post_stage(dt,2)

#######################################################################
# Integrator Steps, Fluid Dynamics only
#######################################################################

class VerletStepFluid(IntegratorStep):
    def __init__(self, rho0=1000.0) -> None:
        self.rho0 = rho0
        super(VerletStepFluid, self).__init__()
    
    def initialize(self, d_idx, d_x, d_x0, d_y, d_y0, d_u, d_u0, d_v, d_v0, d_rho, d_rho0):
        d_u0[d_idx] = d_u[d_idx]
        d_v0[d_idx] = d_v[d_idx]
        d_x0[d_idx] = d_x[d_idx]
        d_y0[d_idx] = d_y[d_idx]
        d_rho0[d_idx] = d_rho[d_idx]
    
    def stage1(self, d_idx, d_x, d_y, d_u, d_v, d_rho, d_au, d_av, d_arho,
               d_usph, d_vsph, dt):
        dt2 = dt/2.0
        d_x[d_idx] += (d_u[d_idx] + d_usph[d_idx])*dt2
        d_y[d_idx] += (d_v[d_idx] + d_vsph[d_idx])*dt2
        d_u[d_idx] += d_au[d_idx]*dt2
        d_v[d_idx] += d_av[d_idx]*dt2
        d_rho[d_idx] += d_arho[d_idx]*dt2
    
    def stage2(self, d_idx, d_x, d_x0, d_y, d_y0, d_u, d_u0, d_v, d_v0,
               d_rho, d_rho0, d_usph, d_vsph, d_au, d_av, d_arho,  dt):
        d_u[d_idx] = d_u0[d_idx] + dt*d_au[d_idx] 
        d_v[d_idx] = d_v0[d_idx] + dt*d_av[d_idx]
        u2 = 0.5*(d_u0[d_idx]+d_u[d_idx])
        v2 = 0.5*(d_v0[d_idx]+d_v[d_idx])
        d_x[d_idx] = d_x0[d_idx] + dt*(u2 + d_usph[d_idx])
        d_y[d_idx] = d_y0[d_idx] + dt*(v2 + d_vsph[d_idx])
        
        eps = -(d_arho[d_idx]/d_rho[d_idx])*dt
        rho_new = d_rho0[d_idx]*(2.0 - eps)/(2.0 + eps)
        d_rho[d_idx] = rho_new #max(rho_new, self.rho0)

class VerletStepDBC(IntegratorStep):
    def __init__(self, rho0=1000.0) -> None:
        self.rho0 = rho0
        super(VerletStepDBC, self).__init__()
    
    def initialize(self, d_idx, d_rho, d_rho0):
        d_rho0[d_idx] = d_rho[d_idx]
    
    def stage1(self, d_idx, d_rho, d_arho, dt):
        dt2 = dt/2.0
        d_rho[d_idx] += d_arho[d_idx]*dt2
    
    def stage2(self, d_idx, d_rho, d_rho0, d_arho, dt):
        eps = -(d_arho[d_idx]/d_rho[d_idx])*dt
        # rho_new = d_rho0[d_idx] + dt * d_arho[d_idx]
        rho_new = d_rho0[d_idx]*(2.0 - eps)/(2.0 + eps)
        d_rho[d_idx] =  max(rho_new, self.rho0)

class VerletStepInlet(IntegratorStep):
    def __init__(self, rho0=1000.0) -> None:
        self.rho0 = rho0
        super(VerletStepInlet, self).__init__()
        
    def initialize(self, d_idx, d_rho, d_rho0, d_y, d_y0):
        d_y0[d_idx] = d_y[d_idx]
        d_rho0[d_idx] = d_rho[d_idx]
        
    def stage1(self, d_idx, d_v, d_rho, d_arho, d_y, dt):
        dt2 = dt/2.0
        d_rho[d_idx] += d_arho[d_idx]*dt2
        d_y[d_idx] += d_v[d_idx]*dt2
    
    def stage2(self, d_idx, d_rho, d_rho0, d_y, d_y0, d_arho, d_v, dt):
        eps = -(d_arho[d_idx]/d_rho[d_idx])*dt
        rho_new = d_rho0[d_idx]*(2.0 - eps)/(2.0 + eps)
        d_rho[d_idx] = max(rho_new, self.rho0)
        d_y[d_idx] = d_y0[d_idx] + dt*d_v[d_idx]

#######################################################################
# Integrator Steps, Thermal aspect only
#######################################################################

class VerletStepThermal(IntegratorStep):
    def initialize(self, d_idx, d_Temp, d_Temp0):
        d_Temp0[d_idx] = d_Temp[d_idx]
    
    def stage1(self, d_idx, d_Temp, d_aTemp, dt):
        dt2 = dt/2.0
        d_Temp[d_idx] += d_aTemp[d_idx]*dt2
    
    def stage2(self, d_idx, d_Temp, d_Temp0, d_aTemp, dt):
        d_Temp[d_idx] = d_Temp0[d_idx] + d_aTemp[d_idx]*dt

class VerletStepThermalSolidificationPure(IntegratorStep):
    def initialize(self, d_idx, d_Temp, d_Temp0):
        d_Temp0[d_idx] = d_Temp[d_idx]
        
    def stage1(self, d_idx, d_Temp, d_aTemp, dt):
        dt2 = dt*0.5
        d_Temp[d_idx] += dt2*d_aTemp[d_idx]
    
    def stage2(self, d_idx, d_Temp, d_aTemp, d_Temp0, d_qstate,
               d_cp, d_Lh, d_Tsol, dt):
        Tnew = d_Temp0[d_idx] + dt*d_aTemp[d_idx]
        if (Tnew < d_Tsol[d_idx]) and (d_qstate[d_idx] < 1.0):
            St = (
                d_cp[d_idx]*(d_Tsol[d_idx]-Tnew)/
                d_Lh[d_idx]
            )
            qnew = d_qstate[d_idx] + St
            if qnew >= 1.0:
                d_qstate[d_idx] = 1.0
                d_Temp[d_idx] = d_Tsol[d_idx] - d_Lh[d_idx]*(
                    qnew - 1.0
                )/d_cp[d_idx]
            else:
                d_qstate[d_idx] = qnew
                d_Temp[d_idx] = d_Tsol[d_idx]
        else:
            d_Temp[d_idx] = Tnew

class VerletStepThermalSolidificationAndMeltPure(IntegratorStep):
    def initialize(self, d_idx, d_Temp, d_Temp0):
        d_Temp0[d_idx] = d_Temp[d_idx]
        
    def stage1(self, d_idx, d_Temp, d_aTemp, dt):
        dt2 = dt*0.5
        d_Temp[d_idx] += dt2*d_aTemp[d_idx]
    
    def stage2(self, d_idx, d_Temp, d_aTemp, d_Temp0, d_qstate,
               d_cp, d_Lh, d_Tsol, dt):
        Tnew = d_Temp0[d_idx] + dt*d_aTemp[d_idx]
        # Solidification
        if (Tnew < d_Tsol[d_idx]) and (d_qstate[d_idx] < 1.0):
            St = (
                d_cp[d_idx]*(d_Tsol[d_idx]-Tnew)/
                d_Lh[d_idx]
            )
            qnew = d_qstate[d_idx] + St
            if qnew >= 1.0:
                d_qstate[d_idx] = 1.0
                d_Temp[d_idx] = d_Tsol[d_idx] - d_Lh[d_idx]*(
                    qnew - 1.0
                )/d_cp[d_idx]
            else:
                d_qstate[d_idx] = qnew
                d_Temp[d_idx] = d_Tsol[d_idx]
        # Melt
        elif (Tnew > d_Tsol[d_idx]) and (d_qstate[d_idx] > 0.0):
            St = (
                d_cp[d_idx]*(d_Tsol[d_idx]-Tnew)/
                d_Lh[d_idx]
            )
            qnew = d_qstate[d_idx] + St
            if qnew <= 0.0:
                d_qstate[d_idx] = 0.0
                d_Temp[d_idx] = d_Tsol[d_idx] - d_Lh[d_idx]*qnew/d_cp[d_idx]
            else:
                d_qstate[d_idx] = qnew
                d_Temp[d_idx] = d_Tsol[d_idx]
        else:
            d_Temp[d_idx] = Tnew

#######################################################################
# Integrator Steps, Thermal and FluidMech aspects coupled
#######################################################################

class VerletStepCoupled(IntegratorStep):
    """Updates density, velocity, position, temperature.
    Use Symplectic Verlet algorithm.
    """
    def __init__(self, rho0=1000.0) -> None:
        self.rho0 = rho0
        super(VerletStepCoupled, self).__init__()
        
    def initialize(self, d_idx, d_Temp, d_Temp0, d_u, d_u0, d_v, d_v0,
                   d_x, d_x0, d_y, d_y0, d_rho, d_rho0):
        d_u0[d_idx] = d_u[d_idx]
        d_v0[d_idx] = d_v[d_idx]
        d_x0[d_idx] = d_x[d_idx]
        d_y0[d_idx] = d_y[d_idx]
        d_rho0[d_idx] = d_rho[d_idx]
        d_Temp0[d_idx] = d_Temp[d_idx]
    
    def stage1(self, d_idx, d_Temp, d_aTemp, d_u, d_usph, d_v, d_vsph,
                   d_x, d_y, d_au, d_av, d_rho, d_arho, dt):
        dt2 = dt/2.0
        d_x[d_idx] += (d_u[d_idx] + d_usph[d_idx])*dt2
        d_y[d_idx] += (d_v[d_idx] + d_vsph[d_idx])*dt2
        d_u[d_idx] += d_au[d_idx]*dt2
        d_v[d_idx] += d_av[d_idx]*dt2
        d_rho[d_idx] += d_arho[d_idx]*dt2
        d_Temp[d_idx] += d_aTemp[d_idx]*dt2
    
    def stage2(self, d_idx, d_Temp, d_Temp0, d_aTemp, d_u, d_u0, d_v,
               d_v0, d_x, d_x0, d_y, d_y0, d_rho, d_rho0, d_au, d_av,
               d_usph, d_vsph, d_arho, dt):
        d_u[d_idx] = d_u0[d_idx] + dt*d_au[d_idx] 
        d_v[d_idx] = d_v0[d_idx] + dt*d_av[d_idx]
        u2 = 0.5*(d_u0[d_idx]+d_u[d_idx])
        v2 = 0.5*(d_v0[d_idx]+d_v[d_idx])
        d_x[d_idx] = d_x0[d_idx] + dt*(u2 + d_usph[d_idx])
        d_y[d_idx] = d_y0[d_idx] + dt*(v2 + d_vsph[d_idx])
        
        eps = -(d_arho[d_idx]/d_rho[d_idx])*dt
        rho_new = d_rho0[d_idx]*(2.0 - eps)/(2.0 + eps)
        d_rho[d_idx] = max(rho_new, self.rho0)
        d_Temp[d_idx] = d_Temp0[d_idx] + d_aTemp[d_idx]*dt

class VerletStepCoupledSolidificationPure(IntegratorStep):
    """Updates density, velocity, position, temperature and
    solid fraction by considering solidification only of
    a pure substance. Use Symplectic Verlet algorithm.
    """
    def __init__(self, rho0=1000.0) -> None:
        self.rho0 = rho0
        super(VerletStepCoupledSolidificationPure, self).__init__()
        
    def initialize(self, d_idx, d_Temp, d_Temp0, d_u, d_u0, d_v, d_v0,
                   d_x, d_x0, d_y, d_y0, d_rho, d_rho0):
        d_u0[d_idx] = d_u[d_idx]
        d_v0[d_idx] = d_v[d_idx]
        d_x0[d_idx] = d_x[d_idx]
        d_y0[d_idx] = d_y[d_idx]
        d_rho0[d_idx] = d_rho[d_idx]
        d_Temp0[d_idx] = d_Temp[d_idx]
        
    def stage1(self, d_idx, d_Temp, d_aTemp, d_u, d_usph, d_v, d_vsph,
                   d_x, d_y, d_au, d_av, d_rho, d_arho, dt):
        dt2 = dt/2.0
        d_x[d_idx] += (d_u[d_idx] + d_usph[d_idx])*dt2
        d_y[d_idx] += (d_v[d_idx] + d_vsph[d_idx])*dt2
        d_u[d_idx] += d_au[d_idx]*dt2
        d_v[d_idx] += d_av[d_idx]*dt2
        d_rho[d_idx] += d_arho[d_idx]*dt2
        d_Temp[d_idx] += d_aTemp[d_idx]*dt2
    
    def stage2(self, d_idx, d_Temp, d_Temp0, d_aTemp, d_u, d_u0, d_v,
               d_v0, d_x, d_x0, d_y, d_y0, d_rho, d_rho0, d_au, d_av,
               d_usph, d_vsph, d_arho, d_Tsol, d_qstate, d_cp, d_Lh,
               dt):
        d_u[d_idx] = d_u0[d_idx] + dt*d_au[d_idx] 
        d_v[d_idx] = d_v0[d_idx] + dt*d_av[d_idx]
        u2 = 0.5*(d_u0[d_idx]+d_u[d_idx])
        v2 = 0.5*(d_v0[d_idx]+d_v[d_idx])
        d_x[d_idx] = d_x0[d_idx] + dt*(u2 + d_usph[d_idx])
        d_y[d_idx] = d_y0[d_idx] + dt*(v2 + d_vsph[d_idx])
        
        eps = -(d_arho[d_idx]/d_rho[d_idx])*dt
        rho_new = d_rho0[d_idx]*(2.0 - eps)/(2.0 + eps)
        d_rho[d_idx] = max(rho_new, self.rho0)
        Tnew = d_Temp0[d_idx] + dt*d_aTemp[d_idx]
        if (Tnew < d_Tsol[d_idx]) and (d_qstate[d_idx] < 1.0):
            St = (
                d_cp[d_idx]*(d_Tsol[d_idx]-Tnew)/
                d_Lh[d_idx]
            )
            qnew = d_qstate[d_idx] + St
            if qnew >= 1.0:
                d_qstate[d_idx] = 1.0
                d_Temp[d_idx] = d_Tsol[d_idx] - d_Lh[d_idx]*(
                    qnew - 1.0
                )/d_cp[d_idx]
            else:
                d_qstate[d_idx] = qnew
                d_Temp[d_idx] = d_Tsol[d_idx]
        else:
            d_Temp[d_idx] = Tnew

class VerletStepCoupledSolidificationAndMeltPure(IntegratorStep):
    """Updates density, velocity, position, temperature and
    solid fraction by considering solidification/melting of
    a pure substance. Use Symplectic Verlet algorithm.
    """
    def __init__(self, rho0=1000.0) -> None:
        self.rho0 = rho0
        super(VerletStepCoupledSolidificationAndMeltPure, self).__init__()
        
    def initialize(self, d_idx, d_Temp, d_Temp0, d_u, d_u0, d_v, d_v0,
                   d_x, d_x0, d_y, d_y0, d_rho, d_rho0):
        d_u0[d_idx] = d_u[d_idx]
        d_v0[d_idx] = d_v[d_idx]
        d_x0[d_idx] = d_x[d_idx]
        d_y0[d_idx] = d_y[d_idx]
        d_rho0[d_idx] = d_rho[d_idx]
        d_Temp0[d_idx] = d_Temp[d_idx]
        
    def stage1(self, d_idx, d_Temp, d_aTemp, d_u, d_usph, d_v, d_vsph,
                   d_x, d_y, d_au, d_av, d_rho, d_arho, dt):
        dt2 = dt/2.0
        d_x[d_idx] += (d_u[d_idx] + d_usph[d_idx])*dt2
        d_y[d_idx] += (d_v[d_idx] + d_vsph[d_idx])*dt2
        d_u[d_idx] += d_au[d_idx]*dt2
        d_v[d_idx] += d_av[d_idx]*dt2
        d_rho[d_idx] += d_arho[d_idx]*dt2
        d_Temp[d_idx] += d_aTemp[d_idx]*dt2
    
    def stage2(self, d_idx, d_Temp, d_Temp0, d_aTemp, d_u, d_u0, d_v,
               d_v0, d_x, d_x0, d_y, d_y0, d_rho, d_rho0, d_au, d_av,
               d_usph, d_vsph, d_arho, d_Tsol, d_qstate, d_cp, d_Lh,
               dt):
        d_u[d_idx] = d_u0[d_idx] + dt*d_au[d_idx] 
        d_v[d_idx] = d_v0[d_idx] + dt*d_av[d_idx]
        u2 = 0.5*(d_u0[d_idx]+d_u[d_idx])
        v2 = 0.5*(d_v0[d_idx]+d_v[d_idx])
        d_x[d_idx] = d_x0[d_idx] + dt*(u2 + d_usph[d_idx])
        d_y[d_idx] = d_y0[d_idx] + dt*(v2 + d_vsph[d_idx])
        
        eps = -(d_arho[d_idx]/d_rho[d_idx])*dt
        rho_new = d_rho0[d_idx]*(2.0 - eps)/(2.0 + eps)
        d_rho[d_idx] = rho_new#max(rho_new, self.rho0)
        Tnew = d_Temp0[d_idx] + dt*d_aTemp[d_idx]
        # Solidification
        if (Tnew < d_Tsol[d_idx]) and (d_qstate[d_idx] < 1.0):
            St = (
                d_cp[d_idx]*(d_Tsol[d_idx]-Tnew)/
                d_Lh[d_idx]
            )
            qnew = d_qstate[d_idx] + St
            if qnew >= 1.0:
                d_qstate[d_idx] = 1.0
                d_Temp[d_idx] = d_Tsol[d_idx] - d_Lh[d_idx]*(
                    qnew - 1.0
                )/d_cp[d_idx]
            else:
                d_qstate[d_idx] = qnew
                d_Temp[d_idx] = d_Tsol[d_idx]
        # Melt
        elif (Tnew > d_Tsol[d_idx]) and (d_qstate[d_idx] > 0.0):
            St = (
                d_cp[d_idx]*(d_Tsol[d_idx]-Tnew)/
                d_Lh[d_idx]
            )
            qnew = d_qstate[d_idx] + St
            if qnew <= 0.0:
                d_qstate[d_idx] = 0.0
                d_Temp[d_idx] = d_Tsol[d_idx] - d_Lh[d_idx]*qnew/d_cp[d_idx]
            else:
                d_qstate[d_idx] = qnew
                d_Temp[d_idx] = d_Tsol[d_idx]
        else:
            d_Temp[d_idx] = Tnew
            
class VerletStepCoupledSolidificationAndMeltPureTVF(IntegratorStep):
    """Updates density, velocity, position, temperature and
    solid fraction by considering solidification/melting of
    a pure substance. Use Symplectic Verlet algorithm.
    """
    def __init__(self, rho0=1000.0) -> None:
        self.rho0 = rho0
        super(VerletStepCoupledSolidificationAndMeltPureTVF, self).__init__()
        
    def initialize(self, d_idx, d_Temp, d_Temp0, d_u, d_u0, d_v, d_v0,
                   d_x, d_x0, d_y, d_y0, d_rho, d_rho0):
        d_u0[d_idx] = d_u[d_idx]
        d_v0[d_idx] = d_v[d_idx]
        d_x0[d_idx] = d_x[d_idx]
        d_y0[d_idx] = d_y[d_idx]
        d_rho0[d_idx] = d_rho[d_idx]
        d_Temp0[d_idx] = d_Temp[d_idx]
        
    def stage1(self, d_idx, d_Temp, d_aTemp, d_u, d_v,
                   d_x, d_y, d_au, d_av, d_rho, d_arho,
                   d_uhat, d_vhat, d_auhat, d_avhat,
                   dt):
        dt2 = dt/2.0
        d_u[d_idx] += d_au[d_idx]*dt2
        d_v[d_idx] += d_av[d_idx]*dt2
        d_uhat[d_idx] = d_u[d_idx] + d_auhat[d_idx]*dt2
        d_vhat[d_idx] = d_v[d_idx] + d_avhat[d_idx]*dt2
        d_x[d_idx] += (d_uhat[d_idx])*dt2
        d_y[d_idx] += (d_vhat[d_idx])*dt2
        d_rho[d_idx] += d_arho[d_idx]*dt2
        d_Temp[d_idx] += d_aTemp[d_idx]*dt2
    
    def stage2(self, d_idx, d_Temp, d_Temp0, d_aTemp, d_u, d_u0, d_v,
               d_v0, d_x, d_x0, d_y, d_y0, d_rho, d_rho0, d_au, d_av,
               d_uhat, d_vhat, d_auhat, d_avhat,
               d_arho, d_Tsol, d_qstate, d_cp, d_Lh,
               dt):
        d_u[d_idx] = d_u0[d_idx] + dt*d_au[d_idx] 
        d_v[d_idx] = d_v0[d_idx] + dt*d_av[d_idx]
        u2 = 0.5*(d_u0[d_idx]+d_u[d_idx])
        v2 = 0.5*(d_v0[d_idx]+d_v[d_idx])
        d_uhat[d_idx] = u2 + d_auhat[d_idx]*dt
        d_vhat[d_idx] = v2 + d_avhat[d_idx]*dt
        d_x[d_idx] = d_x0[d_idx] + dt*(d_uhat[d_idx])
        d_y[d_idx] = d_y0[d_idx] + dt*(d_vhat[d_idx])
        
        eps = -(d_arho[d_idx]/d_rho[d_idx])*dt
        rho_new = d_rho0[d_idx]*(2.0 - eps)/(2.0 + eps)
        d_rho[d_idx] = max(rho_new, self.rho0)
        Tnew = d_Temp0[d_idx] + dt*d_aTemp[d_idx]
        # Solidification
        if (Tnew < d_Tsol[d_idx]) and (d_qstate[d_idx] < 1.0):
            St = (
                d_cp[d_idx]*(d_Tsol[d_idx]-Tnew)/
                d_Lh[d_idx]
            )
            qnew = d_qstate[d_idx] + St
            if qnew >= 1.0:
                d_qstate[d_idx] = 1.0
                d_Temp[d_idx] = d_Tsol[d_idx] - d_Lh[d_idx]*(
                    qnew - 1.0
                )/d_cp[d_idx]
            else:
                d_qstate[d_idx] = qnew
                d_Temp[d_idx] = d_Tsol[d_idx]
        # Melt
        elif (Tnew > d_Tsol[d_idx]) and (d_qstate[d_idx] > 0.0):
            St = (
                d_cp[d_idx]*(d_Tsol[d_idx]-Tnew)/
                d_Lh[d_idx]
            )
            qnew = d_qstate[d_idx] + St
            if qnew <= 0.0:
                d_qstate[d_idx] = 0.0
                d_Temp[d_idx] = d_Tsol[d_idx] - d_Lh[d_idx]*qnew/d_cp[d_idx]
            else:
                d_qstate[d_idx] = qnew
                d_Temp[d_idx] = d_Tsol[d_idx]
        else:
            d_Temp[d_idx] = Tnew


class VerletStepFluidTVF(IntegratorStep):
    """Updates density, velocity, position, temperature and
    solid fraction by considering solidification/melting of
    a pure substance. Use Symplectic Verlet algorithm.
    """
    def __init__(self, rho0=1000.0) -> None:
        self.rho0 = rho0
        super(VerletStepFluidTVF, self).__init__()
        
    def initialize(self, d_idx, d_u, d_u0, d_v, d_v0,
                   d_x, d_x0, d_y, d_y0, d_rho, d_rho0):
        d_u0[d_idx] = d_u[d_idx]
        d_v0[d_idx] = d_v[d_idx]
        d_x0[d_idx] = d_x[d_idx]
        d_y0[d_idx] = d_y[d_idx]
        d_rho0[d_idx] = d_rho[d_idx]
        
    def stage1(self, d_idx, d_u, d_v,
                   d_x, d_y, d_au, d_av, d_rho, d_arho,
                   d_uhat, d_vhat, d_auhat, d_avhat,
                   dt):
        dt2 = dt/2.0
        d_u[d_idx] += d_au[d_idx]*dt2
        d_v[d_idx] += d_av[d_idx]*dt2
        d_uhat[d_idx] = d_u[d_idx] + d_auhat[d_idx]*dt2
        d_vhat[d_idx] = d_v[d_idx] + d_avhat[d_idx]*dt2
        d_x[d_idx] += (d_uhat[d_idx])*dt2
        d_y[d_idx] += (d_vhat[d_idx])*dt2
        d_rho[d_idx] += d_arho[d_idx]*dt2
    
    def stage2(self, d_idx, d_u, d_u0, d_v,
               d_v0, d_x, d_x0, d_y, d_y0, d_rho, d_rho0, d_au, d_av,
               d_uhat, d_vhat, d_auhat, d_avhat,
               d_arho,
               dt):
        d_u[d_idx] = d_u0[d_idx] + dt*d_au[d_idx] 
        d_v[d_idx] = d_v0[d_idx] + dt*d_av[d_idx]
        u2 = 0.5*(d_u0[d_idx]+d_u[d_idx])
        v2 = 0.5*(d_v0[d_idx]+d_v[d_idx])
        d_uhat[d_idx] = u2 + d_auhat[d_idx]*dt
        d_vhat[d_idx] = v2 + d_avhat[d_idx]*dt
        d_x[d_idx] = d_x0[d_idx] + dt*(d_uhat[d_idx])
        d_y[d_idx] = d_y0[d_idx] + dt*(d_vhat[d_idx])
        
        eps = -(d_arho[d_idx]/d_rho[d_idx])*dt
        rho_new = d_rho0[d_idx]*(2.0 - eps)/(2.0 + eps)
        d_rho[d_idx] = max(rho_new, self.rho0)