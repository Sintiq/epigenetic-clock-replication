FROM condaforge/miniforge3:latest

WORKDIR /work

COPY environment.yml /work/environment.yml
RUN conda env create -f /work/environment.yml && conda clean -afy

SHELL ["conda", "run", "-n", "epi-clock-replication", "/bin/bash", "-c"]

COPY . /work

CMD ["conda", "run", "--no-capture-output", "-n", "epi-clock-replication", "make", "all"]
