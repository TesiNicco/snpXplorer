# syntax=docker/dockerfile:1

FROM eclipse-temurin:17-jre-jammy AS java17

FROM python:3.10-slim-bookworm

# The Google Cloud Storage connector installed below requires Java 17.
COPY --from=java17 /opt/java/openjdk /opt/java/openjdk
ENV JAVA_HOME=/opt/java/openjdk
ENV PATH="${JAVA_HOME}/bin:${PATH}"

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
        libpq-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff-dev \
        libjpeg-dev \
        libfribidi-dev \
        libharfbuzz-dev \
        libopenblas0-pthread \
        liblapack3 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# -------------------------------------------------------------------
# Python deps
# -------------------------------------------------------------------
COPY requirements.txt requirements.txt
RUN pip3 install --no-cache-dir -r requirements.txt

# Keep the original requirements file untouched for the standard image.
# Hail requires NumPy 2; replace the NumPy 1-only Matplotlib wheel from requirements.txt.
RUN pip3 install --no-cache-dir gnomad-toolbox==0.0.1 'matplotlib>=3.9,<3.11' \
    && pip3 check \
    && java -version \
    && python3 -c "import hail, gnomad_toolbox, matplotlib, wordcloud"

# Access the public gnomAD Hail tables without Google Cloud authentication.
RUN wget -q https://broad.io/install-gcs-connector -O /tmp/install-gcs-connector.py \
    && python3 /tmp/install-gcs-connector.py --auth-type UNAUTHENTICATED \
    && rm /tmp/install-gcs-connector.py

ENV GNOMAD_SOURCE=hail

# -------------------------------------------------------------------
# App code
# -------------------------------------------------------------------
COPY Exploration /Exploration

# -------------------------------------------------------------------
# R from CRAN (for Debian)
# -------------------------------------------------------------------
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        r-base \
        r-base-dev \
        r-cran-circlize \
        r-cran-ggplot2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# -------------------------------------------------------------------
# R packages from CRAN
# -------------------------------------------------------------------
RUN Rscript -e "packages <- c( \
        'dynamicTreeCut','argparse','viridis','dendextend','pheatmap', \
        'data.table','RColorBrewer','circlize','plotrix','ggplot2' \
    ); install.packages(setdiff(packages, rownames(installed.packages())), \
        repos='https://cloud.r-project.org'); \
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]; \
    if (length(missing)) stop('Missing R packages: ', paste(missing, collapse=', '))"

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
# Original execution with dev server (not for production)
#CMD service redis-server start && python3 -m flask run --host=88.99.91.126 -p 8001
CMD service redis-server start && python3 -m flask run --host=0.0.0.0 -p 8001

# Use the following line for deployment with gunicorn
#CMD service redis-server start && gunicorn -w 1 --threads 4 -b 0.0.0.0:8007 --timeout 300 app:app

# Gunicorn plus gthreaded workers -- 2 workers, 4 threads each
#CMD service redis-server start && gunicorn -k gthread -w 2 --threads 4 --timeout 300 -b 0.0.0.0:8007 app:app

# Gunicorn plus gevent for better performances -- 2 workers
#CMD service redis-server start && gunicorn -k gevent -w 2 --worker-connections 1000 --timeout 300 -b 0.0.0.0:8007 app:app

# If you ever run this locally, usually you'd want:
# CMD service redis-server start && python3 -m flask run --host=0.0.0.0 -p 8007

# To build
# docker build -t snpxplorer --name snpxplorer --network=host .
# To run
# docker run --name snpxplorer --restart unless-stopped -p 127.0.0.1:8001:8001 \
# -v /root/100plus_snpxplorer/definite_standalone/snpXplorer/Exploration:/Exploration \
# -v /root/100plus_snpxplorer/definite_standalone/snpXplorer/Annotation:/Annotation \
# -v /root/100plus_snpxplorer/definite_standalone/snpXplorer/Data:/Data \
# -v /home/nicco/snpxplorer-hail.json:/run/secrets/snpxplorer-hail.json:ro \
# -e GOOGLE_APPLICATION_CREDENTIALS=/run/secrets/snpxplorer-hail.json \
# snpxplorer

# To run interactively
# docker exec -it snpxplorer /bin/bash
