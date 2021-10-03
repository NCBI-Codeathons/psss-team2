FROM bitnami/minideb:buster
MAINTAINER Johannes Köster <johannes.koester@tu-dortmund.de>
ADD . /tmp/repo
WORKDIR /tmp/repo
ENV PATH /opt/conda/bin:${PATH}
ENV LANG C.UTF-8
ENV SHELL /bin/bash
RUN install_packages wget curl bzip2 ca-certificates gnupg2 squashfs-tools git time vim
RUN /bin/bash -c "curl -L https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh > miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh"
RUN /bin/bash -c "conda install -n base -c conda-forge mamba && \
    mamba create -c conda-forge -c bioconda -n snakemake snakemake && \
    source activate snakemake && \
    conda clean --all -y"
RUN echo "source activate snakemake" > ~/.bashrc
ENV PATH /opt/conda/envs/snakemake/bin:${PATH}
