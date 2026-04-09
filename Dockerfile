FROM rocker/r-ver:4.4.2

# System libraries required by R and Python packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Python
    python3 \
    python3-dev \
    python3-pip \
    # R package compilation
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg-dev \
    # Bioconductor / genomics
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libhdf5-dev \
    # General build tools
    build-essential \
    curl \
    git \
    libgmp3-dev \
    libx11-dev \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1

# GCC 13 on Ubuntu 24.04 treats -Wformat-security as an error; some Bioconductor
# packages (e.g. VariantAnnotation) trigger it. Drop -Werror=format-security globally.
RUN mkdir -p /root/.R && \
    echo 'CFLAGS = -g -O2 -fstack-protector-strong -Wformat -Wdate-time -D_FORTIFY_SOURCE=2' > /root/.R/Makevars && \
    echo 'CXXFLAGS = -g -O2 -fstack-protector-strong -Wformat -Wdate-time -D_FORTIFY_SOURCE=2' >> /root/.R/Makevars

WORKDIR /app

# ── Python dependencies ───────────────────────────────────────────────────────
COPY pyproject.toml pyproject.toml
COPY README.md README.md
COPY svig/ svig/
RUN pip install --no-cache-dir --break-system-packages . jupyter

# ── R: CRAN packages ──────────────────────────────────────────────────────────
RUN R -e "install.packages(c( \
    'BiocManager', 'remotes', \
    'plotrix', 'dplyr', 'tidyverse', 'plyr', \
    'lmtest', 'stringr', 'optparse', 'readxl', 'reshape2', \
    'IRkernel' \
  ), repos='https://cloud.r-project.org')"

# ── R: Bioconductor packages ──────────────────────────────────────────────────
RUN R -e "BiocManager::install(c( \
    'BSgenome.Hsapiens.UCSC.hg19', \
    'GenomicRanges', \
    'regioneR', \
    'Rsamtools', \
    'rhdf5', \
    'MutationTimeR' \
  ), ask=FALSE, update=FALSE)"

# ── R: signature.tools.lib from GitHub ───────────────────────────────────────
RUN R -e "remotes::install_github( \
    'Nik-Zainal-Group/signature.tools.lib', \
    ref='0e6b9f5a7cf41ab6083e3b0edf24fc39a595847f', \
    upgrade='never')"

# ── Register R kernel with Jupyter ────────────────────────────────────────────
RUN R -e "IRkernel::installspec(user=FALSE)"

# ── Copy the rest of the project ─────────────────────────────────────────────
COPY . .
