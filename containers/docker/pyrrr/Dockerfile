#------------------------------------------------------------------------------------#
#
#   Created by Alem Gusinac, last modified at 30-09-2024
# 
#   This Dockerfile creates a container suitable for both Python:3.11 and r-base:latest
# 
#   Current available setups by section:
#       1. Linux package installation
#       2. Python:3.11 (see requirements/python_req.txt for overview of packages)
#       3. r-base:latest (current version: 4.4.0)
#       4. PhantomJS (dependent on r-base)
#       5. PEAR
#       6. Seqkit
# 
#   Additional setups can be included after the latest section
#   It is important to include lib/pkgs for your setup in the first section
#
#   NOTE:   - Addition of python packages can be included to the 'python_req.txt' file
#           - Addition of R packages can be included in section 3; under CRAN, BiocManager or Github
#
#------------------------------------------------------------------------------------#

FROM ubuntu:20.04

ARG USER=docker

# set environment without graphical interface
ENV DEBIAN_FRONTEND=noninteractive

#------------------------------------------------------------------------------------#
# 1. Installation of essential linux packages
#------------------------------------------------------------------------------------#
RUN apt-get update && apt-get install -y \
    autoconf automake \
    make \
    pandoc \
    cargo \
    curl wget pgp \ 
    libtool \
    libcurl4-openssl-dev \ 
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    libncurses5-dev \
    libgdbm-dev \
    libnss3-dev \
    libreadline-dev \
    libffi-dev \
    libbz2-dev \
    liblzma-dev \
    fastqc \
    parallel \
    default-jre \
    && rm -rf /var/lib/apt/lists/*

#------------------------------------------------------------------------------------#
# 2. Python:3.11 setup
#------------------------------------------------------------------------------------#
# Getting python 3.11 from source
RUN wget https://www.python.org/ftp/python/3.11.3/Python-3.11.3.tgz \
    && tar -xvf Python-3.11.3.tgz \
    && cd Python-3.11.3 \
    && ./configure --enable-optimizations \
    && make altinstall \
    && rm -rf /tmp/* /var/lib/apt/lists/*

# Install pip3 for Python 3.11
RUN wget https://bootstrap.pypa.io/get-pip.py \
    && python3.11 get-pip.py

# Copy requirements
COPY python_req.txt .

# Python dependencies
RUN pip3 install --upgrade pip \
    && pip3 install -r python_req.txt \
    && rm -rf /tmp/* /var/lib/apt/lists/*

#------------------------------------------------------------------------------------#
# 3. R-base:latest setup
#------------------------------------------------------------------------------------#
# Setting up R-base:latest, includes public key setup and addition of r-base version to apt depository
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" >> /etc/apt/sources.list
# apt update, upgrade and clean after above set-up is required!
# Otherwise it gives an error of broken pipes.
RUN apt-get update && apt-get upgrade -y && apt-get clean
RUN apt-get update && apt-get install -y r-base

# Copy requirements
COPY install2.r .
COPY installBioc.r .
COPY installGithub.r .

# Required package for install2.r
RUN R -e "install.packages('docopt',dependencies=TRUE)"

# R dependencies from CRAN
RUN Rscript install2.r --error --skipinstalled \
    tibble \
    networkD3 \
    webshot \
    htmlwidgets \
    && rm -rf /tmp/downloaded_packages

#------------------------------------------------------------------------------------#
# 4. Setup PhantomJS within PATH
#------------------------------------------------------------------------------------#
RUN mkdir /tmp/phantomjs \
  && curl -L https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-linux-x86_64.tar.bz2 \
          | tar -xj --strip-components=1 -C /tmp/phantomjs \
  && cd /tmp/phantomjs \
  && mv bin/phantomjs /usr/local/bin \
  && cd \
  && apt-get purge --auto-remove -y \
      curl \
  && apt-get clean \
  && rm -rf /tmp/* /var/lib/apt/lists/*

# PhantomJS requires an OpenSSL config even if it's an empty one
ENV OPENSSL_CONF=/opt/openssl.cnf

#------------------------------------------------------------------------------------#
# 5. Setup PEAR (paired-end read merger)
#------------------------------------------------------------------------------------#
COPY pear-0.9.11-linux-x86_64.tar.gz .
RUN tar -xzf pear-0.9.11-linux-x86_64.tar.gz \
    && mv pear-0.9.11-linux-x86_64/bin/pear /usr/local/bin \
    && rm pear-0.9.11-linux-x86_64.tar.gz

#------------------------------------------------------------------------------------#
# 6. Setup Seqkit
#------------------------------------------------------------------------------------#
COPY seqkit_linux_386.tar.gz .
RUN tar -xzf seqkit_linux_386.tar.gz \
    && mv seqkit /usr/local/bin \
    && rm seqkit_linux_386.tar.gz

#------------------------------------------------------------------------------------#
# Final steps
#------------------------------------------------------------------------------------#

RUN useradd -ms /bin/bash ${USER}
USER ${USER}