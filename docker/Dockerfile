# Build image for both Python and R
FROM python:3.11-slim-bullseye
# Create USER
ARG USER=docker
# set environment debian
ENV DEBIAN_FRONTEND=noninteractive

# Essential updates
RUN pip3 install --upgrade pip
RUN apt-get update && apt-get install -y \
    r-base \
    build-essential autoconf automake libtool \
    pandoc wget curl

# Copy requirements
COPY requirements/python_req.txt .
COPY requirements/R_req.R .
COPY requirements/pear-0.9.11-linux-x86_64.tar.gz .

# Python and R dependencies
RUN pip3 install -r python_req.txt
RUN Rscript R_req.R

# Setup PhantomJS within PATH
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
      curl \
  && mkdir /tmp/phantomjs \
  && curl -L https://bitbucket.org/ariya/phantomjs/downloads/phantomjs-2.1.1-linux-x86_64.tar.bz2 \
          | tar -xj --strip-components=1 -C /tmp/phantomjs \
  && cd /tmp/phantomjs \
  && mv bin/phantomjs /usr/local/bin \
  && cd \
  && apt-get purge --auto-remove -y \
      curl \
  && apt-get clean \
  && rm -rf /tmp/* /var/lib/apt/lists/*

# Setup PEAR (paired-end read merger)
RUN tar -xzf pear-0.9.11-linux-x86_64.tar.gz \
    && mv pear-0.9.11-linux-x86_64/bin/pear /usr/local/bin \
    && rm pear-0.9.11-linux-x86_64.tar.gz

RUN useradd -ms /bin/bash ${USER}

# Set added user as default user
USER ${USER}
# PhantomJS requires  an OpenSSL config even if it's an empty one
ENV OPENSSL_CONF=/opt/openssl.cnf