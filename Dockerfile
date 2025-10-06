# Start from rocker base with RStudio Server + tidyverse
FROM rocker/rstudio:4.4.1

# Install system dependencies for R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-pip python3-venv \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libglpk40 libglpk-dev \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

# Create a virtualenv for Python packages
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install Python dependencies
RUN pip install --upgrade pip && \
    pip install numpy pandas scanpy scikit-learn matplotlib seaborn jupyter

# Install R dependencies (from CRAN + BiocManager)
RUN R -e "install.packages(c('devtools','reticulate','BiocManager'))" \
 && R -e "BiocManager::install(c('SingleCellExperiment','scater','scran'))"

# Install hafez package directly from GitHub
RUN R -e "devtools::install_github('mamouzgar/hafez')"

# Expose RStudio Server port
EXPOSE 8787

# Default command runs RStudio Server
CMD ["/init"]

