# Presentation
***CryoPy*** is a set of Python scripts for numerical simulation of **fluid dynamics** and **heat transfer**, primarily based on the **SPH** (Smoothed Particle Hydrodynamics) **method**. This project is organized to allow the modeling of various **cryolava scenarios**, based on the physics of avalanches, fluidized snow, and ice/water mixtures.

This README is not finished yet
What's next:
- Installation markdown
- Usage cases

# Project Structure

```text
CryoPy/
│   .gitignore
│   main_db.py
│
├── apps/
│   ├── DB2D.py
│   ├── DB2DFluidized.py
│   ├── DB2DSnow.py
│   ├── FissureInletAvalanche.py
│   ├── FissureInletMix.py
│   ├── FissureInletSnow.py
│
└── modules/
    ├── FluidDynamics.py
    ├── Geometries.py
    ├── HeatTransfer.py
    ├── Integrators.py
    └── TimeStep.py
```

# Folder Descriptions

- `main_db.py`: Main script to launch a simulation (e.g., Doom Mons with a reservoir).
- `apps/`: Contains various simulation scenarios, each implemented as a class inheriting from `Application`.
    - `DB2D.py`, `DB2DFluidized.py`, `DB2DSnow.py`: 2D simulations of fluids, fluidized materials, or snow.
    - `FissureInletAvalanche.py`, `FissureInletMix.py`, `FissureInletSnow.py`: Simulations of flow through a fissure (avalanches, mixtures, snow).
- `modules/`: Core modules for physics and geometry.
    - `FluidDynamics.py`: Equations and models for fluid dynamics.
    - `Geometries.py`: Geometry generation (fissures, blocks, reservoirs).
    - `HeatTransfer.py`: Heat transfer models.
    - `Integrators.py`: Time integrators for the simulation.
    - `TimeStep.py`: Adaptive time step management.

# Prerequisites
- Python 3.8
- PySPH (and other scientific dependencies: numpy, matplotlib, etc.)

For more information about the installation procedure, please see the INSTALLATION.md

# Example case
An example is given in the project. It simulates a dam break in 2D. To launch it, just run:
```bash
python main_db.py
```
You can modify the physical and geometrical parameters directly from the `main_db.py` file.

# Personalized applications

To create a new scenario, add a file into the `apps/` repository and implement a class inheriting from `Application`. The basic structure is the following:
```python
class MyOwnApp(Application):
    def intitialize(self):
        # define init parameters
        # ...
    
    def create_particles(self):
        # create the particles based on the wanted geometry
        # Some are available in modules/Geometries.py
        # ...
    
    def create_equations(self):
        # Use the equations defined in modules/FluidDynamics.py
        #    or in modules/HeatTransfer.py
        # The Time Step computation is also an equation to
        #    include here.
        # Some are available in modules/TimeSteps.py
    
    def create_solver(self):
        # define which solver to use on which particle group
        # Some integrators and integrator steps are located in
        # modules/Integrators
```
Feel free to use the modules available in the `modules/` repository to define the equations, the geometries and the integrators. Don't hesitate to add your own equations for your own specific cases!


Author
Bastien Bodin

