FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    python3.10 \
    python3-pip \
    python3-dev \
    build-essential \
    gcc \
    g++ \
    gfortran \
    git \
    libopenmpi-dev \
    openmpi-bin \
    libtrilinos-zoltan-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

RUN ln -sf /usr/lib/x86_64-linux-gnu/libtrilinos_zoltan.so \
           /usr/lib/x86_64-linux-gnu/libzoltan.so && \
    ln -sf /usr/lib/x86_64-linux-gnu/libtrilinos_zoltan.a \
           /usr/lib/x86_64-linux-gnu/libzoltan.a

RUN ln -sf /usr/bin/python3.10 /usr/bin/python

RUN pip3 install --upgrade pip&& \
    pip3 install "setuptools==59.8.0" "wheel==0.38.4"

ENV SETUPTOOLS_USE_DISTUTILS=stdlib

ENV ZOLTAN_INCLUDE=/usr/include/trilinos
ENV ZOLTAN_LIBRARY=/usr/lib/x86_64-linux-gnu

RUN pip3 install "numpy==1.23.1"

RUN pip install matplotlib jupyter Beaker

RUN pip install pyproject.toml

RUN pip3 install "cython==0.29.32"
RUN pip3 install "mako==1.2.3"
RUN pip3 install "pytest==7.2.0"

RUN pip3 install --no-build-isolation "cyarray==1.1"

RUN pip3 install --no-build-isolation "compyle==0.8.1"

RUN pip install "meshio==5.3.4"

RUN pip3 install --no-build-isolation "mpi4py==3.1.4"

RUN pip install "scipy==1.9.3"

RUN pip install "h5py==3.7.0"

RUN pip3 install --no-build-isolation "pyzoltan==1.0.1"

RUN pip3 install --no-build-isolation "pysph==v1.0b1"


ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
WORKDIR /framework