#!/bin/bash
WORKDIR="$HOME"
echo "Working directory: ${WORKDIR}"
cd "$WORKDIR" || { echo "Error: Failed to cd to $WORKDIR"; exit 1; }

# Load modules
echo "[1/8] Loading modules on the cluster"
module purge
module load openmpi
module load python/3.8.10

# Check for MPI
MPI_PATH=${MPI_HOME:-$OMPI_HOME}
if [ -z "$MPI_PATH" ]; then
    echo "Error: MPI_HOME or OMPI_HOME undefined. Check that OpenMPI is correctly loaded."
    exit 1
fi
echo "OpenMPI found at: $MPI_PATH"

# Check for Python
if ! command -v python &> /dev/null; then
    echo "Error: Python is not available. Check that the python-3.8 module is correctly loaded."
    exit 1
fi
PYTHON_VERSION=$(python --version 2>&1)
if [[ "$PYTHON_VERSION" != *"3.8"* ]]; then
    echo "Error: Python 3.8 is required. Current version: $PYTHON_VERSION"
    exit 1
fi

# Install Zoltan
echo "[2/8] Installing Zoltan..."
# Define the expected Zoltan paths
ZOLTAN_PATH="${WORKDIR}/Zoltan"

# Check if Zoltan is already installed
if [ -f "${ZOLTAN_PATH}/include/zoltan.h" ] && [ -f "${ZOLTAN_PATH}/lib/libzoltan.a" ]; then
    echo "Zoltan is already installed at ${ZOLTAN_PATH}. Skipping installation."
else
    wget https://github.com/sandialabs/Zoltan/archive/refs/tags/v3.901.tar.gz || {
        echo "Error: Failed to download Zoltan."
        exit 1
    }
    tar xvfz v3.901.tar.gz || {
        echo "Error: Failed to extract Zoltan."
        exit 1
    }
    mkdir -p "${ZOLTAN_PATH}" && cd "${ZOLTAN_PATH}" || {
        echo "Error: Failed to create or enter Zoltan directory."
        exit 1
    }
    ../Zoltan-3.901/configure --with-cflags=-fPIC --enable-mpi --with-mpi="$MPI_PATH" --prefix="${ZOLTAN_PATH}" || {
        echo "Error: Zoltan configure failed."
        exit 1
    }
    make everything || {
        echo "Error: Zoltan make failed."
        exit 1
    }
    make install || {
        echo "Error: Zoltan make install failed."
        exit 1
    }
    echo "Zoltan installed successfully at ${ZOLTAN_PATH}."
fi

# Configure Zoltan's environment variables (only if not already set)
echo "[3/8] Setting Zoltan environment variables..."
if ! grep -q "export ZOLTAN_INCLUDE=${ZOLTAN_PATH}/include" ~/.bashrc; then
    echo "" >> ~/.bashrc
    echo "export ZOLTAN_INCLUDE=${ZOLTAN_PATH}/include" >> ~/.bashrc
    echo "export ZOLTAN_LIBRARY=${ZOLTAN_PATH}/lib" >> ~/.bashrc
    echo "" >> ~/.bashrc
    source ~/.bashrc
    echo "Zoltan environment variables added to ~/.bashrc."
else
    echo "Zoltan environment variables already set in ~/.bashrc."
fi

# Return to working directory
cd "$WORKDIR" || { echo "Error: Failed to return to $WORKDIR"; exit 1; }

# Create and activate virtual environment
echo "[4/8] Creating Python virtual environment..."
python -m venv env_pysph || {
    echo "Error: Failed to create virtual environment."
    exit 1
}
source env_pysph/bin/activate || {
    echo "Error: Failed to activate virtual environment."
    exit 1
}
echo "Using Python from: $(which python)"

# Upgrade pip
echo "[5/8] Upgrading pip..."
pip install --upgrade pip || {
    echo "Error: Failed to upgrade pip."
    exit 1
}

# Install core dependencies
echo "[6/8] Installing dependencies..."
pip install numpy==1.23.1 cython==0.29.32 mako==1.2.3 pytest==7.2.0 || {
    echo "Error: Failed to install core dependencies."
    exit 1
}
pip install compyle==0.8.1 --no-build-isolation || {
    echo "Error: Failed to install compyle."
    exit 1
}
pip install carray==1.1 --no-build-isolation || {
    echo "Error: Failed to install carray."
    exit 1
}
pip install mpi4py==3.1.4 --no-build-isolation || {
    echo "Error: Failed to install mpi4py."
    exit 1
}
pip install meshio==5.3.4 scipy==1.9.3 h5py==3.7.0 || {
    echo "Error: Failed to install meshio, scipy, or h5py."
    exit 1
}

# Install PyZoltan and PySPH
echo "[7/8] Installing PyZoltan..."
pip install git+https://github.com/pypr/pyzoltan.git@v1.0.1 --no-build-isolation || {
    echo "Error: Failed to install PyZoltan."
    exit 1
}

echo "[8/8] Installing PySPH..."
pip install git+https://github.com/pypr/pysph.git@1.0b1 --no-build-isolation || {
    echo "Error: Failed to install PySPH."
    exit 1
}

echo "Installation completed successfully!"