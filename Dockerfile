# aligninator runtime: STAR 2.7.11b + the Python libraries the pipeline needs.
#
#   docker build -t ghcr.io/ckangas2/aligninator:0.2.0 .
#
# Versions are pinned in environment.yml, shared with `-profile conda`.
#
# The pipeline's own scripts are NOT copied in. Nextflow puts the repo's bin/
# on PATH at run time, but appends it; a copy baked into the image would sit
# earlier on PATH and silently run stale code after every edit to bin/.
FROM mambaorg/micromamba:1.5.10-bookworm-slim

LABEL org.opencontainers.image.source="https://github.com/ckangas2/aligninator" \
      org.opencontainers.image.description="Host + virus single-cell alignment with STARsolo"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean --all --yes \
 && rm /tmp/environment.yml

# Explicit PATH so the tools resolve whether or not the micromamba entrypoint
# runs (Nextflow bypasses it).
ENV PATH="/opt/conda/bin:${PATH}"

# Fail the build, not the first pipeline run, if anything is missing.
RUN STAR --version | grep -qx 2.7.11b \
 && python -c "import numpy, scipy, pandas, anndata" \
 && command -v ps zcat >/dev/null
