# Use the official Miniconda3 image as the base image
FROM continuumio/miniconda3

# Update package list and install necessary tools
RUN apt-get update && \
    apt-get install -y gfortran libblas-dev liblapack-dev

# Set up an example Conda environment and install packages
RUN conda create -n myenv python=3.8 -y && \
    echo "source activate myenv" > ~/.bashrc

# Set up the environment for running commands
SHELL ["/bin/bash", "-c", "source activate myenv"]

# Install packages within the Conda environment
RUN conda install -c conda-forge numpy scipy matplotlib petsc4py -y

# Set the environment variable
ENV EMCLAW_PATH /home/$USER/emclaw

# Clone the repository and navigate to the riemann directory
RUN git clone https://github.com/MaxwellGEMS/emclaw.git $EMCLAW_PATH && \
    cd $EMCLAW_PATH/emclaw/riemann

# Activate the Conda environment and build/install the package
RUN source activate myenv && \
    python setup.py build_ext -i && \
    conda develop -n myenv $EMCLAW_PATH

# Set the default command to activate the Conda environment
CMD ["bash", "-i"]
