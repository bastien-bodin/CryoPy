#!/bin/bash

# Installing system dependencies
echo "[1/8] Installing system dependencies..."
sudo apt-get update
sudo apt-get install -y git curl wget build-essential python3.8 python3.8-venv python3-pip openmpi-bin libopenmpi-dev libtrilinos-zoltan-dev

# Configure Zoltan's environment variables
echo "[2/8] Setting environment variables..."
echo "" >> ~/.bashrc
echo "export ZOLTAN_INCLUDE=/usr/include/trilinos" >> ~/.bashrc
echo "export ZOLTAN_LIBRARY=/usr/lib/x86_64-linux-gnu" >> ~/.bashrc
echo "export USE_TRILINOS=1" >> ~/.bashrc
echo "" >> ~/.bashrc
source ~/.bashrc

# Create and activate virtual environment
echo "[3/8] Creating Python virtual environment..."
python3.8 -m venv env_pysph
source env_pysph/bin/activate

# Upgrade pip
echo "[4/8] Upgrading pip..."
pip install --upgrade pip

# Install core dependencies
echo "[5/8] Installing core dependencies..."
pip install numpy==1.23.1 cython==0.29.32 mako==1.2.3 pytest==7.2.0
pip install compyle==0.8.1 carray==1.1 --no-build-isolation

# Install parallel and additional dependencies
echo "[6/8] Installing MPI-related dependencies..."
pip install mpi4py==3.1.4 --no-build-isolation
pip install meshio==5.3.4 scipy==1.9.3 h5py==3.7.0

# Installing PyZoltan
echo "[7/8] Installing PyZoltan..."
pip install git+https://github.com/pypr/pyzoltan.git@v1.0.1 --no-build-isolation

# Installing PySPH
echo "[8/8] Installing PySPH..."
pip install git+https://github.com/pypr/pysph.git@1.0b1 --no-build-isolation