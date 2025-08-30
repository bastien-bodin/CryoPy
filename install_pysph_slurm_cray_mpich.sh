#!/bin/bash
WORKDIR="$HOME"
echo "Working directory: ${WORKDIR}"
cd "$WORKDIR" || { echo "Error: Failed to cd to $WORKDIR"; exit 1; }

# Load modules for Cray system
echo "[1/8] Loading modules on the Cray cluster"
module purge
module load cray-mpich
module load python/3.8.10  # or similar version available on your Cray system

# Check for MPI - Cray systems use different environment variables
MPI_PATH=${CRAY_MPICH_DIR:-$MPICH_DIR}
if [ -z "$MPI_PATH" ]; then
    echo "Error: CRAY_MPICH_DIR or MPICH_DIR undefined. Check that cray-mpich module is correctly loaded."
    exit 1
fi
echo "Cray MPICH found at: $MPI_PATH"

# Check for Python
if ! command -v python &> /dev/null; then
    echo "Error: Python is not available. Check that the python module is correctly loaded."
    exit 1
fi

PYTHON_VERSION=$(python --version 2>&1)
if [[ "$PYTHON_VERSION" != *"3.8"* ]]; then
    echo "Warning: Python 3.8 preferred. Current version: $PYTHON_VERSION"
    echo "Continuing with available Python version..."
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
    
    # Configure Zoltan with Cray-specific settings
    ../Zoltan-3.901/configure \
        --with-cflags="-fPIC" \
        --enable-mpi \
        --with-mpi="$MPI_PATH" \
        --prefix="${ZOLTAN_PATH}" || {
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

# Set up Zoltan environment variables for current session
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

# Create compyle configuration for Cray
echo "[3.5/8] Creating compyle configuration for Cray..."
mkdir -p ~/.compyle
cat > ~/.compyle/config.py << EOF
import os
os.environ['CC'] = 'cc'
os.environ['CXX'] = 'CC'
USE_ZOLTAN = 1
ZOLTAN = "${ZOLTAN_PATH}"
MPI_CFLAGS = ['-I${MPI_PATH}/include']
MPI_LINK = ['-L${MPI_PATH}/lib', '-lmpich']
EOF

echo "compyle configuration created at ~/.compyle/config.py"

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

# Install dependencies
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

# Install mpi4py with Cray-specific settings
echo "Installing mpi4py with Cray MPICH..."
MPICC=cc pip install mpi4py==3.1.4 --no-build-isolation || {
    echo "Error: Failed to install mpi4py."
    exit 1
}

pip install meshio==5.3.4 scipy==1.9.3 h5py==3.7.0 || {
    echo "Error: Failed to install meshio, scipy, or h5py."
    exit 1
}

# Install PyZoltan and PySPH with Cray environment
echo "[7/8] Installing PyZoltan..."

# Function to modify setup.py for Cray compatibility
modify_setup_for_cray() {
    local setup_file="$1"
    if [ -f "$setup_file" ]; then
        echo "Modifying $setup_file for Cray compatibility..."
        # Create backup
        cp "$setup_file" "${setup_file}.backup"

        # Replace gcc with cray
        sed -i "s/compiler = 'gcc'/compiler = 'cray'/g" "$setup_file"

        # Find the start of the get_mpi_flags function
        local func_start=$(grep -n "def get_mpi_flags()" "$setup_file" | cut -d: -f1)

        # Find the end of the get_mpi_flags function by looking for the next 'def' or end of file
        local func_end=$(sed -n "${func_start},\$p" "$setup_file" | grep -n "^[[:space:]]*def " | head -1 | cut -d: -f1)
        if [ -z "$func_end" ]; then
            func_end=$(wc -l < "$setup_file")
        else
            func_end=$((func_start + func_end - 2))
        fi

        # Find the line with "if compiler == 'intel':" inside get_mpi_flags
        local intel_line=$(sed -n "${func_start},${func_end}p" "$setup_file" | grep -n "if compiler == 'intel':" | cut -d: -f1)
        if [ -n "$intel_line" ]; then
            intel_line=$((func_start + intel_line - 1))
            # Find the next "else:" after the intel block
            local else_line=$(sed -n "${intel_line},${func_end}p" "$setup_file" | grep -n "^[[:space:]]*else:" | head -1 | cut -d: -f1)
            if [ -n "$else_line" ]; then
                else_line=$((intel_line + else_line - 1))
                # Check if cray block already exists
                if ! sed -n "${intel_line},${else_line}p" "$setup_file" | grep -q "elif compiler == 'cray':"; then
                    # Insert cray block before the else
                    sed -i "${else_line}i\\
        elif compiler == 'cray':\\
            # Cray-specific MPI flags\\
            link_args = '-L\${MPI_PATH}/lib -lmpich -dynamic'\\
            compile_args = '-I\${MPI_PATH}/include'" "$setup_file"
                    echo "Cray block added to get_mpi_flags."
                else
                    echo "Cray block already present in get_mpi_flags."
                fi
            else
                echo "Could not find 'else:' after 'if compiler == \"intel\":' in get_mpi_flags."
            fi
        else
            echo "Could not find 'if compiler == \"intel\":' in get_mpi_flags."
        fi
        echo "setup.py modified for Cray compatibility."
    else
        echo "File $setup_file not found."
    fi
}

# First attempt to install PyZoltan
echo "Attempting PyZoltan installation (first try)..."
if ! CC=cc CXX=CC pip install git+https://github.com/pypr/pyzoltan.git@v1.0.1 --no-build-isolation; then
    echo "First attempt failed. Trying manual setup.py modification..."
    
    # Clone the repository to modify setup.py
    cd /tmp || { echo "Error: Cannot access /tmp"; exit 1; }
    rm -rf pyzoltan_temp
    git clone https://github.com/pypr/pyzoltan.git pyzoltan_temp || {
        echo "Error: Failed to clone PyZoltan repository."
        exit 1
    }
    cd pyzoltan_temp || { echo "Error: Failed to enter PyZoltan directory."; exit 1; }
    git checkout v1.0.1 || { echo "Error: Failed to checkout v1.0.1 tag."; exit 1; }
    
    # Modify setup.py files
    find . -name "setup.py" -exec bash -c 'modify_setup_for_cray "$0"' {} \;
    
    # Try installation from modified source
    echo "Installing PyZoltan from modified source..."
    CC=cc CXX=CC pip install . --no-build-isolation || {
        echo "Error: Failed to install PyZoltan even after setup.py modification."
        cd "$WORKDIR"
        rm -rf /tmp/pyzoltan_temp
        exit 1
    }
    
    # Cleanup
    cd "$WORKDIR"
    rm -rf /tmp/pyzoltan_temp
    echo "PyZoltan installed successfully after setup.py modification."
else
    echo "PyZoltan installed successfully on first attempt."
fi

echo "[8/8] Installing PySPH..."

# First attempt to install PySPH
echo "Attempting PySPH installation (first try)..."
if ! CC=cc CXX=CC pip install git+https://github.com/pypr/pysph.git@1.0b1 --no-build-isolation; then
    echo "First attempt failed. Trying manual setup.py modification..."
    
    # Clone the repository to modify setup.py
    cd /tmp || { echo "Error: Cannot access /tmp"; exit 1; }
    rm -rf pysph_temp
    git clone https://github.com/pypr/pysph.git pysph_temp || {
        echo "Error: Failed to clone PySPH repository."
        exit 1
    }
    cd pysph_temp || { echo "Error: Failed to enter PySPH directory."; exit 1; }
    git checkout 1.0b1 || { echo "Error: Failed to checkout 1.0b1 tag."; exit 1; }
    
    # Modify setup.py files
    find . -name "setup.py" -exec bash -c 'modify_setup_for_cray "$0"' {} \;
    
    # Try installation from modified source
    echo "Installing PySPH from modified source..."
    CC=cc CXX=CC pip install . --no-build-isolation || {
        echo "Error: Failed to install PySPH even after setup.py modification."
        cd "$WORKDIR"
        rm -rf /tmp/pysph_temp
        exit 1
    }
    
    # Cleanup
    cd "$WORKDIR"
    rm -rf /tmp/pysph_temp
    echo "PySPH installed successfully after setup.py modification."
else
    echo "PySPH installed successfully on first attempt."
fi

echo "Installation completed successfully!"