FROM rocker/r-ver:4.4.0

# 1. Install BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

# 2. Install Bioconductor packages: ensembldb, EnsDb.Hsapiens.v75, rtracklayer
RUN R -e "BiocManager::install(c('ensembldb', 'EnsDb.Hsapiens.v75', 'rtracklayer'))"

# 3. Install CRAN packages: locuszoomr + data.table
RUN R -e "install.packages(c('locuszoomr', 'data.table'), repos='https://cloud.r-project.org')"

# 4. Sanity check: fail the build if any required package is missing
RUN R -e "library(data.table); library(locuszoomr); library(ensembldb); library(EnsDb.Hsapiens.v75); library(rtracklayer)"

# Copy project into container and set working directory
COPY . /workspace
WORKDIR /workspace

# Default command: run your solution
CMD ["bash", "solution.sh"]
