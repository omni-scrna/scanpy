# Scanpy

Scanpy-backed PCA module for omnibenchmark scRNA pipelines.

## Setup

```sh
pixi install
pixi run check
```

`pixi run check` imports all runtime dependencies and prints `OK`. Run it after install to confirm the environment is healthy.

## Usage

### PCA

```sh
pixi run python pca.py \
  --output_dir <dir> \
  --name <name> \
  --normalized.h5 <normalized.h5> \
  --selected.genes <selected.genes.gz> \
  --solver <arpack|randomized> \
  --n_components <int> \
  --random_seed <int>
```

Output: `<output_dir>/<name>_pcas.tsv` — TSV with header `cell_id<TAB>PC1<TAB>...<TAB>PC{n_components}`, one row per cell. See `Embedding` in `src/schemas.py` for the schema.

## Tests

```sh
pixi run test
```

Runs the pytest suite under `tests/`.

## Conda environment export

```sh
pixi run export-env
```

Exports the resolved environment to `envs/scanpy.yml`. The environment is named after the repo root folder.

## Citation

If you use this module in your research, please cite it using the information in `CITATION.cff`.
