# Installation

Use Python 3.10 or 3.11. The tested numerical path uses NumPy 1.x and
MindQuantum 0.12; these bounds are recorded in pyproject.toml.

## Personal Conda channel

Version 0.2.1 is available from y1fei for Linux x86-64 with Python 3.11.
An empty-environment remote installation passed 32 tests, installed CPU CFC,
the MHC truth evaluation and a synthetic Hi-C workflow on 28 September 2026.

~~~sh
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
conda create -n phirefly --override-channels -c y1fei -c conda-forge \
    -c bioconda python=3.11 y1fei::phirefly=0.2.1
conda activate phirefly
phirefly --version
PHIREFLY_BENCHMARK=1 bash \
    "$CONDA_PREFIX/share/phirefly/examples/mhc_smoke/run_mhc_smoke.sh" out/mhc
~~~

Only Phirefly needs to be requested. Conda automatically installs MindQuantum
0.12.0, its pinned Rich 10.9.0 build, htslib, WhatsHap and other dependencies.
The personal-channel MindQuantum/Rich builds repackage unmodified upstream
wheels. There is no post-link pip installation. GPU libraries are not included
or certified by this CPU validation. This is not a Bioconda release.

## Conda plus pip

From the cloned repository:

~~~sh
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
conda env create -f environment.yml
conda activate phirefly
python -I -m pip install '.[qaia]'
phirefly --version
bash examples/mhc_smoke/run_mhc_smoke.sh
~~~

Mamba can also create the environment from this file. The environment installs
Python, NumPy, SciPy, pysam and htslib. MindQuantum is installed by pip, not
silently fetched by a Conda post-link hook. Use a new environment rather than
updating an analysis environment in place.

On HPC, start without a loaded Python module that injects external packages.
Using a Conda Python executable alone does not neutralize PYTHONPATH. Python's
-I flag ignores Python environment variables and user site-packages; pip's
--isolated flag alone does not isolate the Python interpreter. Never fix shared
package permission errors with sudo or --user. The workflow checks bgzip/tabix
and the MindQuantum installation before creating result directories.

If using a venv instead, activate it before running the installed command:

~~~sh
unset PYTHONPATH PYTHONHOME
python -I -m venv .venv-phirefly
source .venv-phirefly/bin/activate
python -I -m pip install '.[qaia]'
phirefly --version
~~~

The venv still needs htslib executables on PATH. Conda's own plugin errors are
separate from Phirefly; repair the Conda installation if they persist after
removing external Python paths.

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

## Conda channels

The recipes in packaging/conda build the three packages published to y1fei.
See the [build and upload instructions](../packaging/conda/README.md).
The full Conda package includes WhatsHap and the MHC smoke/truth files under
share/phirefly/examples/mhc_smoke; PHIREFLY_BENCHMARK=1 enables its HE/SE test.

## Bioconda

On 27 September 2026, the official Anaconda package API returned 404 for
bioconda/phirefly, conda-forge/mindquantum and bioconda/mindquantum.
The y1fei installation above does not imply that an unqualified installation
from Bioconda is available.

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
