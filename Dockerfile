# Use a miniconda base image
FROM continuumio/miniconda3

# Install system dependencies and build tools
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    build-essential \
    libxml2-dev \
    libxslt1-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /app

# Copy the environment file and create the conda environment
COPY environment.yml .
RUN conda config --add channels conda-forge && \
    conda config --set channel_priority strict && \
    (conda env create -f environment.yml || \
     (sleep 5 && conda env create -f environment.yml) || \
     (sleep 10 && conda env create -f environment.yml))

# Make sure the environment is activated
SHELL ["conda", "run", "-n", "metacell_vae_env", "/bin/bash", "-c"]

# Copy and install SEACells
COPY SEACells/ /app/SEACells/
RUN cd SEACells && pip install -e .

# Copy the application code
COPY scripts/ /app/scripts/
COPY run_pipeline.py .
COPY run_pipeline.sh .

# Copy data directory (if exists)
COPY data/ /app/data/

# Create output directories
RUN mkdir -p /app/seacells_output/batch_corrected && \
    mkdir -p /app/seacells_output/clustering_analysis && \
    mkdir -p /app/metacell_vae_output

# Set permissions
RUN chmod +x run_pipeline.sh

# Set the entrypoint to run the pipeline
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "metacell_vae_env", "python", "run_pipeline.py"]
