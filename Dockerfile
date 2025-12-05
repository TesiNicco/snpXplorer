# syntax=docker/dockerfile:1

FROM python:3.8-slim-bullseye

WORKDIR /

# System packages: Python tooling, tabix, Redis, R deps
RUN apt-get update \
    && apt-get install -y \
        wget \
        tabix \
        redis-server \
        sed \
        procps \
        git \
        gnupg2 \
        software-properties-common \
        libcurl4-openssl-dev \
        libxml2-dev \
        libssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libfribidi-dev \
        libharfbuzz-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Python deps
COPY requirements.txt requirements.txt
RUN pip3 install --no-cache-dir -r requirements.txt

# Copy app
COPY Exploration /Exploration

# Add CRAN key & repo for Debian bullseye and install R
RUN apt-get update \
    && apt-get install -y gnupg2 ca-certificates \
    && gpg --keyserver keyserver.ubuntu.com \
           --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' \
    && gpg --armor --export '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' \
         > /etc/apt/trusted.gpg.d/cran_debian_key.asc \
    && echo "deb http://cloud.r-project.org/bin/linux/debian bullseye-cran40/" \
         > /etc/apt/sources.list.d/cran.list \
    && apt-get update \
    && apt-get install -y r-base \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# R packages
RUN Rscript -e "install.packages(c('dynamicTreeCut','argparse','viridis','dendextend', 'pheatmap', 'data.table','RColorBrewer','circlize','plotrix'), repos='https://cloud.r-project.org', quiet=TRUE)"

# Expose Redis port
EXPOSE 6379

# Env
ENV FLASK_APP=Exploration/app.py

# Start Redis + Flask
# Specify the command to run when the container starts
CMD service redis-server start && python3 -m flask run --host=82.165.237.220 -p 8007
#CMD python3 -m flask run --host=82.165.237.220 -p 8001

