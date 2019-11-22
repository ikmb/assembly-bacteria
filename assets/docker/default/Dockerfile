FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for IKMB Bacterial assembly pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN apt-get install -y libgl1-mesa-glx && apt-get clean -y
ENV PATH /opt/conda/envs/assembly-bacteria-1.0/bin:$PATH

