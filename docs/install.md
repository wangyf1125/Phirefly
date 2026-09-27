# Installation

Use Python 3.10 or 3.11. The tested numerical path uses NumPy 1.x and
MindQuantum 0.12; these bounds are recorded in pyproject.toml.

## Conda plus pip

From the cloned repository:

~~~sh
conda env create -f environment.yml
conda activate phirefly
python -m pip install '.[qaia]'
phirefly --version
bash examples/mhc_smoke/run_mhc_smoke.sh
~~~

Mamba can also create the environment from this file. The environment installs
Python, NumPy, SciPy, pysam and htslib. MindQuantum is installed by pip, not
silently fetched by a Conda post-link hook. Use a new environment rather than
updating an analysis environment in place.

The example starts with staged observations. The qaia-marked pytest tests also
exercise small indexed BAM-to-VCF inputs without truth, with one/two extraction
workers. Without Conda, create a Python 3.10/3.11 virtual environment, provide
bgzip and tabix on PATH, then run the same pip installation.

## GPU

~~~sh
python -m pip install '.[gpu]'
PHIREFLY_BACKEND=gpu-float32 PHIREFLY_BATCH_SIZE=100 PHIREFLY_N_ITER=1000 \
    bash examples/mhc_smoke/run_mhc_smoke.sh out/mhc_gpu
~~~

The GPU extra adds PyTorch. A compatible GPU driver and CUDA-capable PyTorch
installation are still required. The tested publication backend was an H100.
The requested backend is not silently replaced by CPU execution.

## Bioconda

On 27 September 2026, the official Anaconda package API returned 404 for
bioconda/phirefly, conda-forge/mindquantum and bioconda/mindquantum.
Do not advertise conda install phirefly as an existing installation route.

The packaging/bioconda/meta.yaml.in file is a full-workflow recipe template,
not a working or submitted recipe. Before opening a Bioconda PR:

1. Publish the tested source commit and a versioned source archive.
2. Replace __RELEASE_SHA256__ with the archive's real SHA256.
3. Package MindQuantum in an accepted channel, or agree on an explicitly
   limited core-only package with maintainers. Do not hide pip installs in hooks.
4. Render/build/test the recipe with Bioconda tooling.
5. Submit the recipe to bioconda/bioconda-recipes; maintainers review it.
6. Advertise the Conda install command only after the package is built and indexed.

Phirefly itself is pure Python and can use noarch: python. This does not
make compiled dependencies or GPU support platform-independent.
See the [Bioconda guidelines](https://bioconda.github.io/contributor/guidelines.html).
