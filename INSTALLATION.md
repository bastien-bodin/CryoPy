Installation procedure for PySPH Model
------
This document outlines the dependencies required to run the PySPH model, based on the [official PySPH installation guide](https://pysph.readthedocs.io/en/main/index.html). The following setup was tested on a Dell Precision 3561 running Ubuntu 20.04. We strongly recommend using a Linux-based OS for this installation.

> ⚠️ **Important**:
> The installation process can be sensitive to package versions. To avoid compatibility issues, strictly follow the versions specified below and use a dedicated virtual environment.


# Prerequisites
These are the prerequisites we recommend as it was installed on thesespecific versions and generated no incompatibilities
- Operating System: Linux (Ubuntu 20.04 or higher)
- C/C++ compiler
- Python Version: 3.8.10 [(Download here)](https://www.python.org/downloads/release/python-3810/)
- Virtual Environment: Required to isolate dependencies


# Setting Up the Virtual Environment
## Create the Virtual Environment
We recommend creating a dedicated virtual environment named `env_PySPH` to avoid conflicts with other projects.
1. Navigate to your preferred directory
```bash
cd /path/to/your/project
```
2. Create the virtual environment
```bash
python3.8 -m venv env_PySPH
```
3. Activate the environment
```bash
source env_PySPH/bin/activate
```
>✅ **Best Practice**:
> Always activate the environment before installing packages or running the model:
> ```bash
> source env_PySPH/bin/activate
> ```

# Installing Dependencies
## Core Dependencies
Install the following packages in the specified order to minimize compatibility issues:
- Numpy 1.23.1
- Cython 0.29.32
- Mako 1.2.3
- Cyarray 1.1
- compyle 0.8.1
- pytest 7.2.0

All these packages are available in the [requirement_bases](./requirements_bases.txt) file.
The installation can be done as follows:
```bash
$ pip install -r requirements_bases.txt
```

## Parallel simulations and MPI
As the SPH method is suited for CFD problems, running simulations can become quite computationally intensive. Therefore, we recommend to use this model with the parallel installation.
For that, a couple more dependencies are needed.
- MPI, with libraries such as openMPI or MPICH
- Zoltan, a library for dynamic load balancing
- A lot more Python packages

We develop the procedures in the following

### Installing the Zoltan library (local)
On a local machine, we recommend using openMPI, as the installation procedure is easier. To do so, just run on a terminal:
```bash
$ sudo apt-get install openmpi-bin libopenmpi-dev
```
Then, to install Zoltan, just run:
```bash
$ sudo apt-get install libtrilinos-zoltan-dev
```
Then, you just need to set some variables in the `~/.bashrc` file
```bash
export ZOLTAN_INCLUDE=/usr/include/trilinos
export ZOLTAN_LIBRARY=/usr/lib/x86_64-linux-gnu
export USE_TRILINOS=1
```
> ✅ Don't forget to reload the `.bashrc` with:
> ```bash
> $ source ~/.bashrc
> ```

### Installing the Zoltan library on a cluster
If you want more computational power, you would probably want this model to run on a supercomputer. As you probably won't have the sudo rights on this machine, the procedure is a bit different.
1. Getting and extracting Zoltan 3.901:
```bash
$ wget https://github.com/sandialabs/Zoltan/archive/refs/tags/v3.901.tar.gz
$ tar xvfz v3.901.tar.gz
```
A `Zoltan-3.901/` directory should have appeared.

2. Configure a Zoltan directory
```bash
$ mkdir Zoltan; cd Zoltan
$ ../Zoltan-3.901/configure --with-cflags=fPIC --enable-mpi \
    --with-mpi=/path/to/mpi --prefix=/absolute/path/to/Zoltan
```
> You can check for the location of the mpi library with `mpicc --showme` or `mpicc -show` depending on if you are using `openMPI` or `MPICH` respectively.

3. Build the Zoltan library
```bash
$ make everything
$ make install
```
> A *lot* of warnings will appear, that's normal
4. Set the environment variables
```bash
export ZOLTAN_INCLUDE=/path/to/Zoltan/include
export ZOLTAN_LIBRARY=/path/to/Zoltan/lib
```
And *voilà!*

## Installing the overall Python dependencies
Once Zoltan is installed, you can *finally* install ***more*** dependencies.
Create and activate a virtual environment, as described above in the serial installation, then install the following packages:
 | Package      | Version |
 |--------------|---------|
 | pip          | upgrade |
 | numpy        | 1.23.1  |
 | matplotlib   | -       |
 | jupyter      | -       |
 | Beaker       | -       |
 | pyproject.toml | -     |
 | cython       | 0.29.32 |
 | mako         | 1.2.3   |
 | pytest       | 7.2.0   |
 | carray       | 1.1*    |
 | compyle      | 0.8.1*  |
 | meshio       | 5.3.4   |
 | mpi4py       | 3.1.4*  |
 | scipy        | 1.9.3   |
 | h5py         | 3.7.0   |

`*` : Installation avec `--no-build-isolation`

For `Cray` architectures, a small modification is necessary.
Create a `config.py` file in `~/.compyle` with the following content (replace to what suits your environment):
```python
import os

os.environ['CC'] = 'cc'
os.environ['CXX'] = 'CC'

USE_ZOLTAN = 1
ZOLTAN = "/path/to/Zoltan"

MPI_CFLAGS = ['...']
MPI_LINK = ['...']
```

### Installing pyZoltan
To install pyZoltan, we must get the 1.0.1 version. To do so, we must obtain the archive onth [Github page](https://github.com/pypr/pyzoltan/tree/v1.0.1). If the the cluster is not on a Cray architecture, just run the `setup.py`.

If not, you might have to modify the `setup.py` file with the following informations:
- Replace compiler = ’gcc’ with compiler = ’cray’;
- add an `elif compiler == 'cray'` statement and define the `link_args` and `compile_args` variables as defined in the `config.py` file

# Installing PySPH (finally)
Once everything is *hopefully* installed, you can install pysph. The compatible version is [1.0b1](https://github.com/pypr/pysph/tree/1.0b1)

The procedure is quite similar to PyZoltan