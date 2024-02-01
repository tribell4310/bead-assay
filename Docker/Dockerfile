FROM python:3.8-bookworm

# Define working directory
WORKDIR /app

# Update alpine package manager and add programs
RUN apt-get update && apt-get install -y bash rsync tmux htop which wget vim
RUN apt-get install -y make automake gcc g++ subversion python3-dev
RUN apt-get install -y imagemagick
RUN apt-get install -y ffmpeg libsm6 libxext6

# Install base python with cython
RUN python -m pip install --upgrade pip setuptools wheel
RUN python -m pip install Cython
RUN python -m pip install numpy

# Precompile the bead_assay packages
COPY *.py .
COPY *.pyx .
RUN python build.py build_ext --inplace

# Install miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
 	bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3 && \
 	rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=$PATH:/opt/miniconda3/bin
RUN echo "conda init" >> ~/.bashrc

# Create virtual environment with dependencies
RUN conda create --name bead python=3.8
RUN conda install -y --name bead -c anaconda matplotlib
RUN conda install -y --name bead -c anaconda imageio
RUN conda install -y --name bead -c anaconda scikit-image
RUN conda install -y --name bead -c anaconda scipy
RUN conda install -y --name bead -c anaconda Cython
RUN conda install -y --name bead -c conda-forge opencv

# Copy other scripts
COPY *.sh .
COPY *.csv .

# Prep data directory
RUN mkdir /data

# Initiate interactive session
CMD ["bash"]
