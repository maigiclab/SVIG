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
    git \
    && rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1

WORKDIR /app

# ── R dependencies via renv ───────────────────────────────────────────────────
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')"

COPY src/renv.lock renv.lock
RUN R -e "renv::restore(lockfile='renv.lock', prompt=FALSE)"

# ── Python dependencies via pyproject.toml ────────────────────────────────────
COPY pyproject.toml pyproject.toml
# Copy the svig package so flit can read it during install
COPY svig/ svig/

RUN pip install --no-cache-dir .

# ── Copy the rest of the project ─────────────────────────────────────────────
COPY . .
