FROM rocker/rstudio:4.6.0

RUN usermod -l hafez rstudio && \
    usermod -d /home/hafez -m hafez
RUN echo "hafez:hafez" | chpasswd

# System dependencies for compiled R packages (gfortran needed for distutils/ElPiGraph.R)
RUN apt-get update && apt-get install -y --no-install-recommends \
    gfortran \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libglpk40 libglpk-dev \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install hafez and all dependencies (CRAN + GitHub Remotes are resolved automatically)
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org')" \
 && R -e "remotes::install_github('mamouzgar/hafez', build_vignettes=FALSE)"

# Smoke-test: verify the package loads
RUN R -e "library(hafez); cat('hafez loaded successfully\n')"

EXPOSE 8787
CMD ["/init"]
