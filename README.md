   # TIPs

TIPs CLI is a command-line wrapper around the TIPs workflow. It standardizes the workflow into a reproducible and portable CLI and unifies all dependencies via Docker. You only need to prepare the **raw data** and download the provided Docker images and resources to run the full pipeline step-by-step or with a single command.

We also provide a **beta Nextflow workflow** that can be executed using the same Docker images and resources.

---

## Features Overview

- `tips init`: Initialize the sample directory structure; optionally symlink resources into the sample directory and generate configuration templates.
- `tips prepare`: Convert Thermo `.raw` files to `mzML/mgf` inside Docker.
- `tips denovo`: Run built-in batch scripts to perform de novo sequencing on `mgf` files using **InstaNovo / Casanovo / PepNet**.
- `tips denovo-te`: Aggregate de novo peptides, BLAST against TE libraries, and generate TE FASTA.
- `tips search-te`: Perform TE-only database searches using Comet / MSFragger / MSGF+, followed by TE-only FDR and canonical removal.
- `tips integrate-te`: Integrate multi-engine results (iProphet / Philosopher), perform strict BLAST-based canonical removal, and optionally score with MixMHCpred.
- `tips run-all`: One-click execution from `init` to `integrate-te`.

---

## Quick Start (One-Click Run)

Assume you have:

- `SAMPLE/`: A sample directory (containing `.raw` files or future `mzML/mgf` files)
- `/path/to/resources/`: The resources directory
- Required Docker images already available

Resources and Docker bundles can be downloaded from Zenodo:  
https://zenodo.org/records/18589565

```bash
tips run-all SAMPLE \
  --resources /path/to/resources
```

Key outputs after completion:

- `SAMPLE/denovo/` – de novo sequencing results
- `SAMPLE/denovo/te/Denovo_TE_SoftMerged_Merged.fasta` – TE FASTA
- `SAMPLE/search/` – TE-only search results and intermediates
- `SAMPLE/integrate/final_TE_results.tsv` – final integrated results

---

## Installation

```bash
git clone <REPO_URL>
cd <REPO_DIR>/CLI
pip install -e .
```

Verify installation:

```bash
tips --help
```

---

## Docker Images

TIPs CLI orchestrates the workflow; all core tools run inside Docker containers.

### Required Images

- Core image (RAW conversion, blastp, search engines, TPP / Philosopher):
  - `tips-core:0.2`
- De novo images (used internally by `tips denovo`):
  - `tips-denovo-casanovo:0.2`
  - `tips-denovo-pepnet:0.2`
  - `tips-denovo-instanovo:0.2`

### Load Images

```bash
zstd -dc tips_all_cu121_0.2.tar | docker load
```

Verify:

```bash
docker images | grep tips
```

---

## Resources Directory Structure

Resources include canonical FASTA files, contaminants, TE data, and search engine parameters. Multiple samples may share the same resources directory.

Minimum required structure:

```text
resources/
  db/
    UP000005640_9606_processed.fasta
    Contaminants.fasta
    human_blastpdb/
      human_proteome_UP000005640_noERV_LINE.pin
      human_proteome_UP000005640_noERV_LINE.psq
      human_proteome_UP000005640_noERV_LINE.phr

  te/
    TE_class_dic.npy
    blastdb/
      hg38_rmsk_6_frame_singleLine_done.npy
      hg38_rmsk_6_frame_singleLine_done.fasta*

  search_params/
    comet.params
    msfragger.params
    msgfplus.params
```

`tips init --resources` will symlink these into `SAMPLE/resources/` and generate `config/denovo_te.yaml`.

---

## Sample Directory Layout

```text
SAMPLE/
  raw/
  mzML/
  mgf/
  denovo/
    casanovo/
    pepnet/
    instanovo/
    tags/
    te/
  search/
  integrate/
  config/
    denovo_te.yaml
  resources/
  logs/
  tmp/
```

---

## Step-by-Step Usage

### 1. Initialize

```bash
tips init SAMPLE --resources /path/to/resources
```

Overwrite existing config if needed:

```bash
tips init SAMPLE --resources /path/to/resources --overwrite-config
```

---

### 2. RAW Conversion

```bash
tips prepare SAMPLE --image tips-core:0.2 --formats mzml,mgf --stage move --jobs 4
```

---

### 3. De Novo Sequencing

```bash
tips denovo SAMPLE --tool all --gpu auto
```

Examples:

```bash
tips denovo SAMPLE --tool instanovo --gpu auto
tips denovo SAMPLE --tool casanovo --gpu -1
tips denovo SAMPLE --tool pepnet --gpu 0
```

Casanovo supports parameter files:

```bash
tips denovo SAMPLE --tool casanovo --param-file /path/to/casanovo.yaml
```

---

### 4. Generate TE FASTA

```bash
tips denovo-te SAMPLE --image tips-core:0.2 --tools instanovo,casanovo,pepnet
```

---

### 5. TE-Only Search

```bash
tips search-te SAMPLE --image tips-core:0.2 --engines comet,msfragger,msgfplus --fdr 0.03
```

---

### 6. Integration

```bash
tips integrate-te SAMPLE --image tips-core:0.2 --fdr 0.03 --threads 20
```

Optional MixMHCpred:

```bash
tips integrate-te SAMPLE --hla-alleles A0201,A6801,B1302
```

Final output:

```text
SAMPLE/integrate/final_TE_results.tsv
```

---

## Reproducibility and Logs

Each step records:

- `*_commands.sh` – reproducible command scripts
- `*.log` – execution logs

---

## FAQ

**Docker permission denied**  
Ensure your user has Docker permissions.

**GPU not available in container**  
Check `nvidia-smi` and NVIDIA Container Toolkit. Use `--gpu -1` to force CPU.

**Missing resources**  
`tips init --resources` validates required paths and will fail early.

**Missing canonical BLAST DB**  
Provide it in resources, or allow `search-te` to build one automatically.

---

## Command Reference

```bash
tips --help
tips init --help
tips prepare --help
tips denovo --help
tips denovo-te --help
tips search-te --help
tips integrate-te --help
tips run-all --help
```
