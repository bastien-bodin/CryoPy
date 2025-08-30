# CryoPy: A Python Framework for Cryolava Simulation using SPH

**CryoPy** is a Python-based numerical simulation toolkit for **fluid dynamics** and **heat transfer**, primarily using the **Smoothed Particle Hydrodynamics (SPH)** method. It is designed to model **cryolava scenarios**, including avalanches, fluidized snow, and ice/water mixtures.

---

## 📂 Project Structure

```text
CryoPy/
│
├── .gitignore
├── main_db.py                # Main script to launch simulations
├── INSTALLATION.md           # Detailed installation instructions
│
├── apps/                     # Simulation scenarios
│   ├── DB2D.py               # 2D fluid simulation
│   ├── DB2DFluidized.py      # 2D fluidized material simulation
│   ├── DB2DSnow.py           # 2D snow simulation
│   ├── FissureInletAvalanche.py  # Avalanche simulation through a fissure
│   ├── FissureInletMix.py    # Mixture simulation through a fissure
│   └── FissureInletSnow.py   # Snow simulation through a fissure
│
└── modules/                  # Core physics and geometry modules
    ├── FluidDynamics.py      # Fluid dynamics equations and models
    ├── Geometries.py         # Geometry generation (fissures, blocks, reservoirs)
    ├── HeatTransfer.py       # Heat transfer models
    ├── Integrators.py        # Time integrators
    └── TimeStep.py           # Adaptive time step management
```

---

## 📦 Prerequisites

- **Python 3.8.10**
- **PySPH** (and scientific dependencies: `numpy`, `matplotlib`, etc.)

### ⚙️ Installation Scripts
To simplify setup, use one of the provided scripts:
- `install_pysph.sh`: Local installation with parallel support.
- `install_pysph_slurm_nocray.sh`: Installation on non-Cray clusters with openMPI.
- `install_pysph_slurm_cray_mpich.sh`: Installation on Cray clusters with MPICH.

> ⚠️ **Note**
> If installing **PySPH** on a cluster, ensure the correct **Python** and **MPI** modules are loaded.

For detailed instructions, see **[INSTALLATION.md](./INSTALLATION.md)**.

---

## 🌊 Example Usage: 2D Dam-Break Flow Simulation
An available example simulates the classical dam-break problem in 2D using the Smoothed Particle Hydrodynamics (SPH) method. A column of fluid is initially held behind a virtual dam. When the simulation starts, the dam is removed, and the fluid collapses under gravity, spreading downstream.

To run it, just run:
```bash
(env_PySPH) $ python main_db.py
```
You can customize physical and geometric parameters directly in `main_db.py`.

---

## 🔧 Creating Custom Scenarios

To add a new scenario:
1. Create a file in the `apps/` directory.
2. Implement a class inheriting from `Application`:
   ```python
   class MyOwnApp(Application):
       def initialize(self):
           # Define initial parameters
           pass

       def create_particles(self):
           # Create particles using geometries from `modules/Geometries.py`
           pass

       def create_equations(self):
           # Use equations from:
           # - `modules/FluidDynamics.py` (fluid dynamics)
           # - `modules/HeatTransfer.py` (heat transfer)
           # - `modules/TimeSteps.py` (time step computation)
           pass

       def create_solver(self):
           # Define solvers for particle groups
           # Integrators are available in `modules/Integrators.py`
           pass
   ```
Feel free to extend or modify existing modules for your specific needs!

---

## 📝 Author
**Bastien Bodin**

## Cite
```bibtex
@software{cryopy,
    title={CryoPy: A Python Framework for Cryolava Simulation using SPH},
    author={Bastien Bodin},
    year={2024},
    url={https://github.com/username/CryoPy}
}
```

## License
[BSD-3](./LICENSE.txt)

## 🆘 Support
For questions, bug reports, or feature requests, please:
- Open an issue on GitHub
- Contact: [bastien.bodin@proton.me](mailto:bastien.bodin@proton.me)