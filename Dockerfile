# syntax=docker/dockerfile:1

FROM python:3.10-slim-bullseye

WORKDIR /

# -------------------------------------------------------------------
# System packages: Python tooling, tabix, Redis, build deps for R
# -------------------------------------------------------------------
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        wget \
        tabix \
        redis-server \
        sed \
        procps \
        git \
        gnupg2 \
        ca-certificates \
        software-properties-common \
        build-essential \
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

# -------------------------------------------------------------------
# Python deps
# -------------------------------------------------------------------
COPY requirements.txt requirements.txt
RUN pip3 install --no-cache-dir -r requirements.txt

# -------------------------------------------------------------------
# App code
# -------------------------------------------------------------------
COPY Exploration /Exploration

# -------------------------------------------------------------------
# R from CRAN (for Debian bullseye)
# -------------------------------------------------------------------
RUN apt-get update \
    && apt-get install -y --no-install-recommends gnupg2 ca-certificates \
    && apt-key adv --keyserver keyserver.ubuntu.com \
                   --recv-keys '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' \
    && echo "deb http://cloud.r-project.org/bin/linux/debian bullseye-cran40/" \
         > /etc/apt/sources.list.d/cran.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends r-base \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# -------------------------------------------------------------------
# R packages from CRAN
# -------------------------------------------------------------------
RUN Rscript -e "install.packages( \
        c('dynamicTreeCut','argparse','viridis','dendextend','pheatmap', \
          'data.table','RColorBrewer','circlize','plotrix'), \
        repos='https://cloud.r-project.org', \
        quiet=TRUE \
    )"

# -------------------------------------------------------------------
# Expose Redis port & Flask env
# -------------------------------------------------------------------
EXPOSE 6379

ENV FLASK_APP=Exploration/app.py

# Change directory to Exploration
WORKDIR /Exploration

# -------------------------------------------------------------------
# Start Redis + Flask
# -------------------------------------------------------------------
#CMD service redis-server start && python3 -m flask run --host=82.165.237.220 -p 8007
# Use the following line for deployment with gunicorn
CMD service redis-server start && gunicorn -w 1 --threads 4 -b 0.0.0.0:8007 --timeout 300 app:app

# If you ever run this locally, usually you'd want:
# CMD service redis-server start && python3 -m flask run --host=0.0.0.0 -p 8007