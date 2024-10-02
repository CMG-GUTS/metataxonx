#------------------------------------------------------------------------------------#
#
#   Created by Alem Gusinac, last modified at 26-09-2024
# 
#   This Dockerfile creates a container suitable for fastqc and parallel
#   It builds on top of a slim ubuntu (70MB) base
# 
#   Current available setups by section:
#       1. Fastqc and GNU parallel
# 
#   Additional setups can be included after the latest section
#   It is important to include lib/pkgs for your setup in the first section
#
#------------------------------------------------------------------------------------#

FROM ubuntu:20.04

ARG USER=docker

# set environment without graphical interface
ENV DEBIAN_FRONTEND=noninteractive

#------------------------------------------------------------------------------------#
# 1. Installation of essential linux packages
#------------------------------------------------------------------------------------#
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
    wget \
    unzip \
    fastqc \
    parallel \
    default-jre \
    && apt-get clean

#------------------------------------------------------------------------------------#
# Final steps
#------------------------------------------------------------------------------------#

RUN useradd -ms /bin/bash ${USER}
USER ${USER}