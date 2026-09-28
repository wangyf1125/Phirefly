# Personal Conda channel

Published to y1fei/main on 28 September 2026: Phirefly 0.2.1,
MindQuantum 0.12.0 and Rich 10.9.0. A clean remote-channel installation passed
32 tests and the installed MHC truth smoke; package URLs and hashes were checked.

These are buildable recipes, not Bioconda submissions. Start with Linux x86-64
and Python 3.11. Phirefly is noarch Python; its MindQuantum dependency is a
platform-specific repack of the unmodified official CPython 3.11 wheel. Its URL,
SHA256 and version-specific LICENSE/NOTICE are fixed in the recipe. The recipe
tests CPU CFC. GPU libraries and hardware are not certified by this build.
Rich 10.9.0 is also repackaged: MindQuantum pins that version, but the available
conda-forge builds target older Python versions. The upstream wheel and dependency
constraints are retained. The Phirefly package includes WhatsHap 2.8 and the MHC
example with its truth files under share/phirefly/examples/mhc_smoke/.

Use a separate build environment containing conda-build, conda-index,
anaconda-client and python-build. On HPC, give dependency solving sufficient
memory and set CONDA_PKGS_DIRS outside a quota-limited home directory.

## Build and test

Run from the repository; WORK must be a new absolute directory outside it:

~~~sh
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
export WORK=/absolute/path/to/new-build
mkdir -p "$WORK/source"
python -I -m build --sdist --outdir "$WORK/dist"
sha256sum "$WORK/dist/phirefly-0.2.1.tar.gz"
tar -xzf "$WORK/dist/phirefly-0.2.1.tar.gz" -C "$WORK/source"
export PHIREFLY_SOURCE_DIR="$WORK/source/phirefly-0.2.1"
conda build packaging/conda/rich --override-channels -c conda-forge \
    --croot "$WORK/build" --output-folder "$WORK/channel" --no-anaconda-upload
conda index "$WORK/channel"
conda build packaging/conda/mindquantum --override-channels \
    -c "file://$WORK/channel" -c conda-forge \
    --croot "$WORK/build" --output-folder "$WORK/channel" --no-anaconda-upload
conda index "$WORK/channel"
conda build packaging/conda/phirefly --override-channels \
    -c "file://$WORK/channel" -c conda-forge -c bioconda \
    --croot "$WORK/build" --output-folder "$WORK/channel" --no-anaconda-upload
conda index "$WORK/channel"
conda create -p "$WORK/clean-env" --override-channels \
    -c "file://$WORK/channel" -c conda-forge -c bioconda python=3.11 phirefly=0.2.1
conda run -p "$WORK/clean-env" phirefly --version
conda run -p "$WORK/clean-env" env PHIREFLY_BENCHMARK=1 \
    bash "$WORK/clean-env/share/phirefly/examples/mhc_smoke/run_mhc_smoke.sh" "$WORK/mhc-test"
~~~

The package build runs tests including a real CPU CFC batch, a tiny BAM-to-VCF
workflow, and Hi-C contact-BAM/edge/orientation/apply commands. Run the MHC smoke
again from the clean environment outside the source checkout before uploading.
Keep package hashes, dependency records and logs alongside the source archive.

## Upload

After validation, authenticate in a terminal, never in a committed script:

~~~sh
anaconda --at anaconda.org login
anaconda --at anaconda.org whoami
anaconda upload --user y1fei "$WORK"/channel/noarch/rich-*.conda
anaconda upload --user y1fei "$WORK"/channel/linux-64/mindquantum-*.conda
anaconda upload --user y1fei "$WORK"/channel/noarch/phirefly-*.conda
~~~

Do not use --force to replace a released package; increment its build number.
After all three packages are publicly indexed, verify a second empty installation:

~~~sh
conda create -n phirefly --override-channels -c y1fei -c conda-forge \
    -c bioconda python=3.11 y1fei::phirefly=0.2.1
conda activate phirefly
PHIREFLY_BENCHMARK=1 bash \
    "$CONDA_PREFIX/share/phirefly/examples/mhc_smoke/run_mhc_smoke.sh" out/mhc
~~~

The Anaconda.org channel is y1fei, not the GitHub username.
Uploading to one's own public channel does not require a Bioconda PR;
it does require account authorization and compliance with the hosting service.
Future versions must also pass upload and remote tests before being advertised.
