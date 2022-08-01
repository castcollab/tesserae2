FROM mambaorg/micromamba:0.25.1 AS builder
COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp/tesserae2/

RUN micromamba install -y -f /tmp/tesserae2/env-build.yml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --upgrade build && cd /tmp/tesserae2/ && COMPILATION_ARCH=haswell python3 -m build

FROM mambaorg/micromamba:0.25.1 AS install
COPY --from=builder --chown=$MAMBA_USER:$MAMBA_USER /tmp/tesserae2/ /tmp/tesserae2
RUN micromamba install -y -f /tmp/tesserae2/env.yml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install /tmp/tesserae2/dist/*.whl
