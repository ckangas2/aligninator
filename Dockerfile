# aligninator runtime: STAR 2.7.11b + the Python validation/conversion layer.
#
#   docker build -t ghcr.io/ckangas2/aligninator:0.2.0 .
#
# Versions are pinned in environment.yml, the single source of truth shared
# with `-profile conda`.
FROM mambaorg/micromamba:1.5.10-bookworm-slim

LABEL org.opencontainers.image.source="https://github.com/ckangas2/aligninator" \
      org.opencontainers.image.description="Host + virus single-cell alignment with STARsolo"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean --all --yes \
 && rm /tmp/environment.yml

# Nextflow also puts the pipeline's bin/ on PATH; baking it in makes the
# image usable on its own (docker run ... solo_to_h5ad.py --help).
COPY --chown=$MAMBA_USER:$MAMBA_USER bin/ /opt/aligninator/bin/

# Explicit PATH so the tools resolve whether or not the micromamba entrypoint
# runs (Nextflow and `docker run --entrypoint` both bypass it).
ENV PATH="/opt/conda/bin:/opt/aligninator/bin:${PATH}"

# Fail the build, not the first pipeline run, if anything is missing.
RUN STAR --version | grep -qx 2.7.11b \
 && python -c "import numpy, scipy, pandas, anndata" \
 && command -v ps zcat solo_to_h5ad.py >/dev/null
