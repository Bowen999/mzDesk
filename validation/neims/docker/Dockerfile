# 原代码的运行环境：Python 3.6 + TensorFlow 1.13.2 + RDKit 2019.03.4（原论文时期的版本），amd64
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8
USER root
RUN micromamba create -y -n neims -c conda-forge python=3.6 rdkit=2019.03.4 pip numpy=1.16 && micromamba clean -a -y
ENV PATH=/opt/conda/envs/neims/bin:$PATH
RUN pip install --no-cache-dir tensorflow==1.13.2 absl-py six "protobuf<3.21" "matplotlib<3.4" scipy pandas
WORKDIR /work
