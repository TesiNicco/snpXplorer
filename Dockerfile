# syntax=docker/dockerfile:1

FROM python:3.8-slim-buster

WORKDIR /

# Install required system packages by Python 3.9 (including Tabix)
RUN apt-get update \
    && apt-get install -y \
	wget \
	tabix \
	redis-server \
	sed \
	procps \
	git \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

# Link phantomjs
RUN ln -s /Exploration/phantomjs /usr/local/bin

# Install other system packages for R
RUN apt-get update \
    && apt-get install -y \
	gnupg2 \
	software-properties-common \
	libcurl4-openssl-dev \
	libxml2-dev \
	openjdk-11-jdk \
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

# Add CRAN GPG key
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys B8F25A8A73EACF41

# Add CRAN repository
RUN echo "deb https://cloud.r-project.org/bin/linux/debian buster-cran40/" >> /etc/apt/sources.list

# Install R
RUN apt-get update \
    && apt-get install -y r-base

# Set JAVA_HOME environment variable -- uncomment
#RUN R CMD javareconf && Rscript -e "install.packages('rJava', repos='http://cran.rstudio.com/')"

# Install additional R packages -- uncomment
#RUN Rscript -e "install.packages(c('RCurl', 'data.table', 'stringr', 'parallel', 'ggplot2', 'BiocManager', 'mailR', 'bedr', 'grDevices', 'plyr', 'viridis', 'ggplot2', 'plotrix', 'pheatmap', 'dynamicTreeCut', 'dendextend', 'RColorBrewer', 'htmlwidgets', 'wordcloud2', 'tidytext', 'dplyr', 'webshot', 'tibble', 'devtools', 'ggsci', 'circlize'), verbose=F, repos='http://cran.rstudio.com/')"

# Install some packages through bioconductor -- uncomment
#RUN R -e "BiocManager::install(c('GenomicRanges', 'rtracklayer', 'LDlinkR', 'gprofiler2'))"
#RUN R -e "devtools::install_github('JosephCrispell/basicPlotteR')"

# Copy the rest of your application files
COPY Exploration /Exploration

# Expose Redis port
EXPOSE 6379

# Set environment variables
ENV FLASK_APP=Exploration/app.py

# Specify the command to run when the container starts
CMD service redis-server start && python3 -m flask run --host=82.165.237.220 -p 8007
#CMD python3 -m flask run --host=82.165.237.220 -p 8001

