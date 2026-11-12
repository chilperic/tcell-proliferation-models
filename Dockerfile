# Reproducible R runtime for the T-cell proliferation models repository.
# Pinning policy:
#   - R is pinned by the rocker/r-ver image tag.
#   - R packages are installed from a fixed CRAN snapshot.
FROM rocker/r-ver:4.4.3

ENV CRAN_REPO=https://packagemanager.posit.co/cran/2025-04-01
ENV R_ENVIRON_USER=/workspace/tcell-proliferation-models/.Renviron

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace/tcell-proliferation-models
COPY . .

RUN printf 'CRAN_REPO=%s\n' "$CRAN_REPO" > .Renviron \
    && Rscript tools/setup-r-packages.R \
    && R CMD INSTALL . \
    && Rscript tools/check-project.R \
    && Rscript scripts/run_model.R examples/reproducible_run.json /tmp/tcell-repro

CMD ["Rscript", "tools/check-project.R"]
