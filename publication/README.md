# Publication records (v0.2.3)

- `benchmark_runs.tsv`: the 264 BAM-to-indexed-VCF runs behind the main
  comparison (four methods, 22 HG002 autosomes, three repetitions). Average
  chromosomes equally within a repetition; report the mean and sample SD of the
  three repetition means for Core and E2E. Times are seconds; errors are percentages.
- `quantum_errors.tsv`: chromosome-mean HE and SE. The full evaluation includes
  all regions. The matched evaluation omits the same 253 observed chr19 SNPs
  for IBM, D-Wave and CFC.
- `ibm_chr22.tsv`: nine chr22 jobs, with minimum sampled energy, error counts
  and denominators. The Fez repeat 1 job also contributed chr22 to the autosomal
  IBM result.

Run `python publication/replay_chr22.py` from the checkout to reproduce the
selected energy for Fez repeat 1 from its measured bitstrings. The bundled file
contains the original-scale edges and the 512 measurement counts for each of
51 components. Qiskit bitstrings are stored with qubit 0 on the right; bits 0
and 1 map to spins +1 and -1. The script uses NumPy only. It does not run a QPU
or compute HE/SE from truth. Expected total energy: **-3776.312746**.
