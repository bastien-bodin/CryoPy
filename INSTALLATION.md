# Dependencies for PySPH Model
This document outlines the dependencies required to run the PySPH model, based on the [official PySPH installation guide](https://pysph.readthedocs.io/en/main/index.html). The following setup was tested on a Dell Precision 3561 running Ubuntu 20.04. We strongly recommend using a Linux-based OS for this installation.

>⚠️ Important:
The installation process can be sensitive to package versions. To avoid compatibility issues, strictly follow the versions specified below and use a dedicated virtual environment.


# Prerequisites
These are the prerequisites we recommend as it was installed on thesespecific versions and generated no incompatibilities
- Operating System: Linux (Ubuntu 20.04 or higher)
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
>✅ Best Practice:
Always activate the environment before installing packages or running the model

# Installing Dependencies
3.1 Core Dependencies
Install the following packages in the specified order to minimize compatibility issues:


## Core dependencies
The core dependencies are the packages/software needed only if you want to "play" with the model on small cases in serial on a personnal laptop or desktop.
These are:
- Numpy 1.23.1
- Cython 0.29.32
- Mako 1.2.3
- Cyarray 1.1
- compyle
- pytest 7.2.0

All these packages are available in the [requirement_bases](./requirements_bases.txt) file.
The installation can be done as follows:
```bash
pip install -r requirements_bases.txt
```