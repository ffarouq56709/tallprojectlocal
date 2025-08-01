# Base image with R, Bioconductor, and Shiny
FROM bioconductor/bioconductor_docker:3.20

# Install system dependencies for Python and R
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    curl \
    git \
    zlib1g-dev \
    libncurses5-dev \
    libreadline-dev \
    libsqlite3-dev \
    libssl-dev \
    libbz2-dev \
    liblzma-dev \
    libffi-dev \
    libgdbm-dev \
    libgdbm-compat-dev \
    uuid-dev \
    libexpat1-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libx11-dev \
    libgl1 \
    libtiff-dev \
    gdal-bin \
    libgdal-dev \
    libudunits2-dev \
    lsb-release \
    ca-certificates \
    gfortran && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Build Python 3.10.13 with --enable-shared
RUN wget https://www.python.org/ftp/python/3.10.13/Python-3.10.13.tgz && \
    tar -xzf Python-3.10.13.tgz && \
    cd Python-3.10.13 && \
    ./configure --enable-optimizations --enable-shared && \
    make -j$(nproc) && \
    make altinstall && \
    cd .. && rm -rf Python-3.10.13*

# Ensure shared lib is registered for reticulate
RUN echo "/usr/local/lib" >> /etc/ld.so.conf.d/python3.10.conf && ldconfig

# Create Python virtual environment
RUN /usr/local/bin/python3.10 -m pip install --upgrade pip setuptools virtualenv && \
    /usr/local/bin/python3.10 -m virtualenv /venv

# Set environment variables for reticulate
ENV RETICULATE_PYTHON=/venv/bin/python
ENV PATH="/venv/bin:$PATH"

# Copy app and Python dependencies
COPY shiny-app /srv/shiny-server/
COPY deps /deps

# Add this line to ensure the risk model is copied into the image
COPY deps/TALLSorts/TALLSorts/models /models

# Install Python packages
RUN /venv/bin/pip install -r /deps/TALLSorts-requirements.txt && \
    cd /deps/TALLSorts && /venv/bin/pip install .

# Install required R packages
RUN R -e "install.packages(c('shiny', 'ggplot2', 'DT', 'git2r', 'readr', 'dplyr', 'tidyr', 'reticulate', 'data.table'), repos='https://cloud.r-project.org')"

# Expose Shiny's default port
EXPOSE 3838

# Start the app
CMD ["Rscript", "/srv/shiny-server/app.R"]

