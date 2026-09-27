# Release checks

The source version is 0.2.0. Do not overwrite the earlier release tag.
Public release and Bioconda availability are separate from local tests.

## Validation

~~~sh
conda env create -f environment.yml
conda activate phirefly
python -m pip install -e '.[qaia,test]' build twine
phirefly --version
phirefly run --help
phirefly benchmark --help
pytest -q
bash examples/mhc_smoke/run_mhc_smoke.sh out/mhc_cpu
python -m build
python -m twine check dist/*
~~~

Repeat installation from the built wheel outside the checkout. CI tests Python
3.10/3.11 and a MindQuantum CPU workflow. GPU validation is a separate test;
a CPU smoke does not verify it. Keep graph identities and VCF comparisons against
the frozen publication source, and record any numerical differences.

## Local validation

Local validation on 27 September 2026 passed in a clean Conda Python 3.11
environment: 27 tests (none skipped), CPU MHC smoke, tiny BAM-to-VCF workflows,
wheel/sdist build, twine checks and wheel import outside the checkout.
Fixed-input MHC comparisons with the frozen single-round benchmark source also
passed on CPU and RTX A5000 at 100 candidates, 1000 iterations and seed 23.
All 22 checks passed per backend: graph matrices, phaselet mapping, raw objective,
VCF records, solution/debug tables, summary fields and configuration column order.
Timing fields were excluded. Comparisons were within each backend, not a claim
of identical CPU and GPU solutions. Tested Python source checksums matched the
final checkout.

The earlier H100 regression preceded the final runner/output simplification;
it is not a test of this exact source snapshot. The final GPU regression used
Python 3.10, MindQuantum 0.12.0 and PyTorch 2.13.0+cu130. These checks are release
regressions, not fresh all-chromosome timing benchmarks. Remote GitHub Actions
and clean Python 3.10 validation still need their CI results; local Python 3.11
success is not a substitute.

## Repository contents

Publish phirefly/, tests/, examples/, .github/, documentation, packaging
metadata and license. Generated figures, jobs, environments and credentials stay
outside the software commit. Research results belong in a versioned evidence
archive with source checksums. The checkout's figures/ is ignored, not deleted.

Review git diff and git diff --cached before committing. A matching tested
commit, tag, dependency record and permanent source/results archive are required
for the paper. A Bioconda PR needs the real release tarball checksum and all
dependencies in accepted channels; the template is not PR-ready.

## Structure

Keep the algorithm path linear: extraction, MSF, phaselets, hyperreads, CFC, VCF;
Hi-C adds contact-based orientation. Avoid a plugin framework or a second solver
interface in the release CLI. Existing installed command aliases are retained
for compatibility. The phaselet runner separates loading, construction, solving
and optional evaluation. Its output module only collects statistics and writes
solution/debug tables; VCF export lives in vcf.py. Shared text and coordinate
helpers live in core/. No additional algorithm modules are needed.
