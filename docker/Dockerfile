# Build image for both Python and R
FROM python:3.11-slim-bullseye
# Create USER
ARG USER=docker
# set environment debian
ENV DEBIAN_FRONTEND=noninteractive

# Essential updates
RUN pip3 install --upgrade pip
RUN apt-get update && apt-get install -y \
    r-base

# Copy requirements
#COPY requirements.txt .
COPY requirements.R .

# Python and R dependencies
RUN pip3 install pandas numpy scipy
RUN Rscript requirements.R

# Setup cutadapt
RUN pip3 install cutadapt

# Setup biom
RUN pip3 install biom-format h5py

# Setup multiqc
RUN pip3 install multiqc

RUN useradd -ms /bin/bash ${USER}

# Set added user as default user
USER ${USER}