FROM mambaorg/micromamba:1.5.8

COPY --chown=$MAMBA_USER:$MAMBA_USER docker/environment.yml /tmp/environment.yml

ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV CONDA_PKGS_DIRS=/tmp/conda_pkgs

RUN micromamba create -y -n perturbseq -f /tmp/environment.yml && \
    micromamba clean --all --yes && \
    rm -rf /tmp/conda_pkgs /tmp/environment.yml

WORKDIR /work

# Sanity checks
RUN micromamba run -n perturbseq python -c "import importlib.metadata as md; print('scanpy', md.version('scanpy')); print('anndata', md.version('anndata'))" && \
    micromamba run -n perturbseq STAR --version