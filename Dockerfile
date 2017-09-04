################################################################################
### 
### [CONTAINER CORE FUNCTIONS]: 
###     install "Tool - Batch Correction" Galaxy tool (and required third part softwares, libraries, ...).
### [NOTE]
###     please refer to README.md and about_docker.md files for further informations
### 
################################################################################

################################################################################
### fix parent containter
FROM ubuntu:16.04

################################################################################
### set author
MAINTAINER Nils Paulhe <nils.paulhe@inra.fr>

################################################################################
### sets the environment variables
ENV TOOL_VERSION = "v2.0.3"
ENV CONTAINER_VERSION = 0.1

LABEL version = "${CONTAINER_VERSION}"
LABEL tool_version = "${TOOL_VERSION}"

################################################################################
### install third part tools 

# install R
RUN apt-get update && \
    apt-get install -y \ 
    r-base \
    libcurl4-openssl-dev \
    libxml2-dev \
    git

# init R env. (Docker)
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

# install R libs
RUN Rscript -e "install.packages('batch', dep=TRUE)"
RUN Rscript -e "install.packages('ade4', dep=TRUE)"
RUN Rscript -e "source('http://www.bioconductor.org/biocLite.R'); biocLite('pcaMethods')"
RUN Rscript -e "source('http://www.bioconductor.org/biocLite.R'); biocLite('ropls')"

################################################################################
### install core scripts

# init. WORKDIR
RUN [ "mkdir", "/scripts" ]

#
# [NOTE] to add scripts, we have two options: get them from GitHub OR copy them from this directory
# 

# get scripts using Git (option 1)
# RUN cd /scripts && \
#    git clone -b release/${TOOL_VERSION} --recursive https://github.com/workflow4metabolomics/batch_correction

# copy scripts files from this directory (option 2)
COPY "." "/scripts/"

## set WORKDIR
# WORKDIR "/scripts"

## set authorizations
# RUN ["chmod", "a+x", "/scripts/batch_correction_wrapper.R"]

# make tool accessible through PATH
ENV PATH = $PATH:/scripts

################################################################################
### clean
RUN apt-get remove -y git && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/ /tmp/* /var/tmp/*

################################################################################
### Define Entry point script
## ENTRYPOINT ["/scripts/batch_correction_wrapper.R"]

### [END]