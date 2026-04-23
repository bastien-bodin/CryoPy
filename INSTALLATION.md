# PySPH Model Installation Guide

This guide provides step-by-step instructions for installing the PySPH model on a Linux-based system
(tested on Ubuntu 20.04). To avoid compatibility issues, **strictly follow the specified versions**
and use a dedicated virtual environment.

If your distro is too recent (Ubuntu 24.04 or above, for instance), we recommend that you directly
refer to the Docker Installation Part

---

# Direct Installation

Here are the instructions for setting-up the python environment and all necessary dependencies
manually and directly to your distro/cluster. If you encounter any issue, or if your working
environment is too recent, please skip this part and go directly to the Docker Installation

## 📋 Prerequisites

- **Operating System:** Linux (Ubuntu 20.04 or higher recommended)
- **C/C++ Compiler:** Required for building dependencies
- **Python Version:** 3.8.10 ([Download Python 3.8.10](https://www.python.org/downloads/release/python-3810/))
- **Virtual Environment:** Required to isolate dependencies

---

## 🛠️ Setting Up the Virtual Environment

### 1. Create the Virtual Environment
```bash
cd /path/to/your/project
python3.8 -m venv env_PySPH
```

### 2. Activate the Environment
```bash
source env_PySPH/bin/activate
```
> **✅ Best Practice:** Always activate the environment before installing packages or running the model.

---

## 📦 Installing Core Dependencies

Install the following packages **in the specified order** to minimize compatibility issues:

```bash
pip install numpy==1.23.1 cython==0.29.32 mako==1.2.3 pytest==7.2.0
pip install carray==1.1 --no-build-isolation
pip install compyle==0.8.1 --no-build-isolation
```

> **📄 Note:** These dependencies are also listed in the [`requirements_bases.txt`](./requirements_bases.txt) file.

---

## ⚡ Parallel Simulations and MPI

PySPH supports parallel simulations using MPI. Below are the installation steps for both **local machines** and **clusters**.

### 🖥️ Local Installation (OpenMPI + Zoltan)

#### 1. Install OpenMPI
```bash
sudo apt-get install openmpi-bin libopenmpi-dev
```

#### 2. Install Zoltan
```bash
sudo apt-get install libtrilinos-zoltan-dev
```

#### 3. Set Environment Variables
Add the following lines to your `~/.bashrc` file:
```bash
export ZOLTAN_INCLUDE=/usr/include/trilinos
export ZOLTAN_LIBRARY=/usr/lib/x86_64-linux-gnu
export USE_TRILINOS=1
```
Reload the `.bashrc` file:
```bash
source ~/.bashrc
```

---

### 🏢 Cluster Installation (Zoltan from Source)

#### 1. Download and Extract Zoltan
```bash
wget https://github.com/sandialabs/Zoltan/archive/refs/tags/v3.901.tar.gz
tar xvfz v3.901.tar.gz
```

#### 2. Configure Zoltan
```bash
mkdir Zoltan && cd Zoltan
../Zoltan-3.901/configure --with-cflags=-fPIC --enable-mpi \
    --with-mpi=/path/to/mpi --prefix=/absolute/path/to/Zoltan
```
> **🔍 Tip:** Find your MPI path using `mpicc --showme` (OpenMPI) or `mpicc -show` (MPICH).

#### 3. Build and Install Zoltan
```bash
make everything
make install
```
> **⚠️ Note:** Warnings during compilation are normal.

#### 4. Set Environment Variables
```bash
export ZOLTAN_INCLUDE=/path/to/Zoltan/include
export ZOLTAN_LIBRARY=/path/to/Zoltan/lib
```

---

## 🐍 Installing Python Dependencies

Install the following packages **after** setting up Zoltan:

| Package       | Version       | Command                                      |
|---------------|---------------|----------------------------------------------|
| pip           | upgrade       | `pip install --upgrade pip`                  |
| numpy         | 1.23.1        | `pip install numpy==1.23.1`                  |
| matplotlib    | latest        | `pip install matplotlib`                     |
| jupyter       | latest        | `pip install jupyter`                        |
| beaker        | latest        | `pip install beaker`                         |
| cython        | 0.29.32       | `pip install cython==0.29.32`                |
| mako          | 1.2.3         | `pip install mako==1.2.3`                    |
| pytest        | 7.2.0         | `pip install pytest==7.2.0`                  |
| carray        | 1.1           | `pip install cyarray==1.1 --no-build-isolation` |
| compyle       | 0.8.1         | `pip install compyle==0.8.1 --no-build-isolation` |
| meshio        | 5.3.4         | `pip install meshio==5.3.4`                  |
| mpi4py        | 3.1.4         | `pip install mpi4py==3.1.4 --no-build-isolation` |
| scipy         | 1.9.3         | `pip install scipy==1.9.3`                   |
| h5py          | 3.7.0         | `pip install h5py==3.7.0`                    |

> **🔧 Cray Architectures:** Create a `config.py` file in `~/.compyle` with MPI flags and paths.

---

## 📥 Installing PyZoltan

#### 1. Download PyZoltan 1.0.1
```bash
git clone --branch v1.0.1 https://github.com/pypr/pyzoltan.git
cd pyzoltan
```

#### 2. Modify `setup.py` (if needed)
- Replace `compiler = 'gcc'` with `compiler = 'cray'` for Cray architectures.
- Add `elif compiler == 'cray'` and define `link_args` and `compile_args` as per your `config.py`.

#### 3. Install PyZoltan
```bash
pip install .
```

---

## 🎉 Installing PySPH

#### 1. Clone PySPH 1.0b1
```bash
git clone --branch 1.0b1 https://github.com/pypr/pysph.git
cd pysph
```

#### 2. Install PySPH
```bash
pip install .
```

---

## 📝 Summary of Commands

Here is a summary of all the commands seen above. Please use what you want to fit your needs! (avoid copy/paste)

```bash
# Virtual Environment
python3.8 -m venv env_PySPH
source env_PySPH/bin/activate

# Core Dependencies
pip install numpy==1.23.1 cython==0.29.32 mako==1.2.3 pytest==7.2.0
pip install carray==1.1 --no-build-isolation
pip install compyle==0.8.1 --no-build-isolation

# MPI + Zoltan (Local)
sudo apt-get install openmpi-bin libopenmpi-dev libtrilinos-zoltan-dev
echo "export ZOLTAN_INCLUDE=/usr/include/trilinos" >> ~/.bashrc
echo "export ZOLTAN_LIBRARY=/usr/lib/x86_64-linux-gnu" >> ~/.bashrc
echo "export USE_TRILINOS=1" >> ~/.bashrc
source ~/.bashrc

# MPI + Zoltan (Cluster)
wget https://github.com/sandialabs/Zoltan/archive/refs/tags/v3.901.tar.gz
tar xvfz v3.901.tar.gz
mkdir Zoltan && cd Zoltan
../Zoltan-3.901/configure --with-cflags=-fPIC --enable-mpi --with-mpi=/path/to/mpi --prefix=/absolute/path/to/Zoltan
make everything && make install
export ZOLTAN_INCLUDE=/path/to/Zoltan/include
export ZOLTAN_LIBRARY=/path/to/Zoltan/lib

# Python Dependencies
pip install --upgrade pip
pip install numpy==1.23.1 matplotlib jupyter beaker cython==0.29.32 mako==1.2.3 pytest==7.2.0
pip install carray==1.1 --no-build-isolation
pip install compyle==0.8.1 --no-build-isolation
pip install mpi4py==3.1.4 --no-build-isolation
pip install meshio==5.3.4 scipy==1.9.3 h5py==3.7.0

# PyZoltan
git clone --branch v1.0.1 https://github.com/pypr/pyzoltan.git
cd pyzoltan
# Modify setup.py if needed
pip install .

# PySPH
git clone --branch 1.0b1 https://github.com/pypr/pysph.git
cd pysph
pip install .
```

---

## 🚀 Next Steps
- Verify the installation by running the PySPH test suite:
  ```bash
  pysph test -v
  ```
- Explore the [PySPH documentation](https://pysph.readthedocs.io/en/main/index.html) for usage examples.


# Docker Installation (Recommended for Local Use)
 
## What is Docker?
 
Docker is a tool that packages a program and *all its dependencies* into a single,
isolated unit called a **container**. Think of it as a self-contained virtual computer
running inside your machine, pre-configured with the exact libraries CryoPy needs
(Python, PySPH, OpenMPI, etc.). You do not need to install anything else manually,
and it will not interfere with your existing software.
 
> **Why does it matter for CryoPy?**
> PySPH has a complex dependency chain that can be difficult to install on macOS or
> Windows. Docker handles all of that for you, regardless of your operating system.
 
---
 
## Prerequisites
 
You need **Docker Desktop**, the official graphical application that runs Docker on
your machine. Download and install it from the official website:
 
- [macOS / Windows / Linux → docker.com/products/docker-desktop](https://www.docker.com/products/docker-desktop/)
Once installed, open Docker Desktop and leave it running in the background.
You can verify everything is working by opening a terminal and typing:
 
```bash
docker --version
docker compose version
```
 
Both commands should print a version number without errors.
 
---
 
## First-Time Setup: Building the Image
 
An **image** is the blueprint Docker uses to create containers. You only need to
build it once (or again if you modify the `Dockerfile`).
 
From the root of the repository (the folder containing `docker-compose.yml`), run:
 
```bash
docker compose up -d --build
```
 
> The first build takes **10–20 minutes** because it compiles PySPH and its
> dependencies from source. Subsequent starts are instant.
 
What this command does, step by step:
 
| Step | What happens |
|------|--------------|
| `--build` | Builds the Docker image from the `Dockerfile` |
| `up` | Creates and starts the container |
| `-d` | Runs it in the background (detached mode) |
 
---
 
## Running a Simulation
 
Once the container is running, open an interactive shell inside it:
 
```bash
docker compose exec pysph bash
```
 
You are now *inside* the container, in the `/framework` directory — which is a
direct mirror of the `framework/` folder on your local machine. Any file you
create or modify here is immediately visible on your computer, and vice versa.
 
Run a simulation:
 
```bash
python main_db.py
```
 
Output files will be written to `framework/outputs/` and will persist on your
local disk after the container is stopped.
 
---
 
## Stopping the Container
 
When you are done, stop the container with:
 
```bash
docker compose down
```
 
Your files in `framework/` are never deleted by this command.
 
---
 
## Quick Reference
 
| Action | Command |
|--------|---------|
| Build and start | `docker compose up -d --build` |
| Start (after first build) | `docker compose up -d` |
| Open a shell inside | `docker compose exec pysph bash` |
| Stop | `docker compose down` |
| Check running containers | `docker ps` |
 
---
 
## Notes
 
- **CPU only**: the current image runs on CPU. MPI parallelism is supported
  (OpenMPI is installed); see the *Parallel Runs* section for details.
- **No GPU support**: CUDA is not configured in this image.
- **Rebuilding the image**: if you add a new dependency to the `Dockerfile`,
  re-run `docker compose up -d --build`.
---
 
## Using CryoPy on a Computing Cluster (HPC)
 
Most HPC clusters do **not** support Docker for security reasons. The standard
alternative is **Singularity** (also called Apptainer), which can convert a Docker
image into a format safe for shared systems.
 
If you need to run CryoPy on a cluster, refer to the shell scripts provided at the
root of the repository:
 
- `install_pysph_slurm_nocray.sh` — for Slurm clusters with OpenMPI
- `install_pysph_slurm_cray_mpich.sh` — for Cray systems with MPICH
Singularity support and a conversion guide are planned for a future release.