# Reproducible local setup

The most reproducible route is Docker. Native installation is provided only as a convenience because the R version provided by the operating system may vary.

## Option A — Docker, recommended

Requirements: Docker or Podman.

```bash
cd tcell-proliferation-models
docker build -t tcell-proliferation-models:r-4.4.3 .
docker run --rm tcell-proliferation-models:r-4.4.3
```

This builds an image using:

- R pinned through `rocker/r-ver:4.4.3`
- R package installation from a fixed CRAN snapshot configured by `CRAN_REPO`
- repository smoke tests through `tools/check-project.R`

To open an interactive shell inside the environment:

```bash
docker run --rm -it -v "$PWD":/workspace/tcell-proliferation-models tcell-proliferation-models:r-4.4.3 bash
```

Then run:

```bash
Rscript tools/check-project.R
Rscript scripts/01_simulate_smith_martin.R
Rscript scripts/02_simulate_cyton.R
```

## Option B — Docker Compose

```bash
cd tcell-proliferation-models
docker compose up --build
```

## Option C — Native Ubuntu/Debian setup

This is less reproducible than Docker because `apt` controls the installed R interpreter version.

```bash
cd tcell-proliferation-models
chmod +x tools/install_r_ubuntu_24_04.sh
./tools/install_r_ubuntu_24_04.sh
```

For later runs:

```bash
./tools/run_local.sh
```

## Files that control the environment

- `Dockerfile`: pinned R runtime and package setup.
- `.Rprofile`: fixed CRAN snapshot repository.
- `tools/setup-r-packages.R`: installs required R packages.
- `tools/check-project.R`: smoke-test after setup.
- `docker-compose.yml`: one-command Docker execution.
- `.devcontainer/devcontainer.json`: VS Code / Dev Container setup.
