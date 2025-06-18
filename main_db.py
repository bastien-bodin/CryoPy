# Main program for Doom Mons simulation with reservoir instead of inlet, Water only

# User inputs
dp = 0.01/(2.0**0.5)                   # Spacing between particles
x0, y0 = 0.0, 0.0           # Starting coords for making the geometry
Lx, Ly = 4.0, 3.0        # Length and Height of the geometry
lx, ly = 1.0, 2.0           # Width and Depth of the Inlet
rho0 = 1000.0    # Density and viscosity of the fluid
vmax = 0.5                  # Maximum velocity
t_inlet = 2.0            # Max time were the inlet adds particles. If None, then t_inlet = tf
tf = t_inlet               # Final physical time of the simulation
# sources: Engineering toolbox

path2dir = '../outputs/' # You can also set as absolute path
prefix_output = f'Water_tf_{int(tf/60)}'
saving_directory = f'{path2dir}{prefix_output}_lx_{int(lx)}_ly_{int(ly)}_rho{rho0}_output/'

def main():
    from apps.DB2D import DB2D
    
    # Defining the simulation
    MyApp = DB2D(
        fname=prefix_output, output_dir=saving_directory,
        dp = dp,
        x0=x0, y0=y0,
        Lx=Lx, Ly=Ly, lx_fluid=lx, ly_fluid=ly,
        rho0=rho0, vmax=vmax,
        tf=tf, t_inlet=t_inlet
    )
    MyApp.run() # Running the simulation
    return None

if __name__=='__main__':
    main()
