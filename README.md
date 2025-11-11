# TE-derive Peptide Search (TIPs)

This repository provides a configurable Python wrapper around a TE-focused
immunopeptidomics discovery workflow inspired by TIPs.  
It is designed to:

1. Extract high-confidence peptide tags directly from MS/MS spectra using
   multiple state-of-the-art de novo sequencing engines.
2. Build a compact, sample-specific database that mixes:
   - Canonical human proteome (PE=0)
   - Candidate TE-derived proteins / tags (PE=1)
   - Common contaminants
   plus automatic decoy generation.
3. Search the MS/MS data against the database using multiple engines
   (Comet, MSFragger, MS-GF+).
4. Integrate the results statistically (PeptideProphet / iProphet),
   control FDR separately for TE vs canonical peptides,
   remove canonical peptides via short-sequence blastp,
   and optionally score HLA binding with MixMHCpred.

**⚠ Important:** This code exposes the orchestration logic and file handling
in Python, but does *not* bundle external binaries such as Comet, MSFragger,
MS-GF+, TPP, Philosopher, ProteoWizard, MixMHCpred, or NCBI BLAST.
You must install them yourself and point to them in the YAML config.

---

## Installation

```bash
git clone https://github.com/your-org/tips-pipeline.git
cd tips-pipeline

# create environment (example using conda/mamba)
mamba create -n tips-pipeline python=3.10 -y
mamba activate tips-pipeline

pip install -r requirements.txt
pip install -e .   # optional, installs tips_pipeline as a package
