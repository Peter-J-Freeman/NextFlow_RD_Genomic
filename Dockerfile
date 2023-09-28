# Use the continuumio/miniconda3 base image
FROM continuumio/miniconda3:23.3.1-0

# Update the package list and install necessary packages using apt
# Note: I'm running a Mac M1, and these packages do not seem to be available
# for ARM via conda, so using apt to install
RUN apt update && apt install -y \
    bwa \
    bcftools \
    samtools \
    libxcb-xinerama0 \
    wget

# Create the Conda environment and install dependencies
RUN conda init bash
RUN echo "conda activate base" > ~/.bashrc

# Specify Conda channels
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Install Conda packages
RUN conda install -y \
    python=3.10 \
    gatk4==4.3.0.0 \
    nextflow==23.04.1 \
    rtg-tools==3.12.1 \
    fastqc==0.12.1 \
    picard

# Set entrypoint
WORKDIR /home
ENTRYPOINT []
CMD ["tail", "-f", "/dev/null"]
