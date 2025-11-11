# TIPs-Style TE-Derived Peptide Discovery Pipeline (TIPsSearch)

A configurable Python orchestrator for immunopeptidomics focused on discovering **transposable element (TE)-derived antigens**, inspired by TIPs. It wires together de novo tag extraction, sample-specific database construction, multi-engine database search, and integrative post-processing.

**Key capabilities**

1. **De novo tag extraction** from MS/MS using multiple engines; high-confidence subsequences are selected by per-residue scores.
2. **Sample-specific database build:** Human proteome (PE=0) + TE candidates (PE=1) + contaminants + decoys. TE candidates are preferably obtained by *de novo tags → BLAST (short)* against a large TE FASTA → collecting matched entries.
3. **Multi-engine DB search:** Comet / MSFragger / MS-GF+ (configurable).
4. **Post-processing:** TPP (RefreshParser → `xinteract`/iProphet with decoy prefix `rev_`) → table export → charge-state–aware FDR filtering → **short-sequence homology filtering** vs the canonical human proteome (BLAST short, 100% identity across full peptide) → optional HLA-binding prediction (MixMHCpred).

> This repository only contains the Python orchestration and I/O layout. It **does not** ship the external binaries (Comet, MSFragger, MS-GF+, TPP/Philosopher, ProteoWizard, MixMHCpred, NCBI BLAST, etc.). Please install them yourself and provide absolute paths in the YAML config.

---

## System requirements

- **OS:** Linux (expects common GNU tools such as `nohup`, `awk`, `ps`, `grep` for simple polling and parallel BLAST handling).
- **Python:** 3.10 (example environment).
- **Python deps:** see `requirements.txt`. If you are not using `lxml` in your environment, you can safely remove it from `requirements.txt` to slim down the environment.

**Heads-up (known quirk):** `fasta_utils.py` calls `subprocess.check_output(...)` when polling BLAST jobs. Please ensure `import subprocess` is present at the top of that file.

---

## Installation

You have two equally valid ways to run this project. Choose one.

### Option A — Install as a package (recommended)

1) Create and activate an environment
```bash
mamba create -n tips-pipeline python=3.10 -y
mamba activate tips-pipeline
pip install -r requirements.txt
```

2) Add minimal packaging files (e.g., `pyproject.toml`) so editable install works, then:
```bash
pip install -e .
```

3) Run the CLI as a module (relative imports require package context):
```bash
python -m tips_pipeline.cli --config /abs/path/to/config.yaml --stage all
```
`--stage` can be one of: `denovo` | `dbsearch` | `postprocess` | `all`.

> **Minimal `pyproject.toml` (example):**
> ```toml
> [project]
> name = "tips-pipeline"
> version = "0.1.0"
> requires-python = ">=3.10"
> dependencies = [
>   "pandas>=2.0",
>   "numpy>=1.24",
>   "biopython>=1.83",
>   "pyyaml>=6.0",
>   "lxml>=5.2",
> ]
> 
> [project.scripts]
> tipssearch = "tips_pipeline.cli:main"
> ```
> Place all Python files inside a package folder named `tips_pipeline/` (with an `__init__.py`).

### Option B — Run from source (without packaging)

Keep a package-like layout (`tips_pipeline/` containing all `.py` files), then:
```bash
mamba activate tips-pipeline
pip install -r requirements.txt
export PYTHONPATH="$(pwd)"
python -m tips_pipeline.cli --config /abs/path/to/config.yaml --stage all
```

---

## Inputs & working directory layout

We refer to `sample.path` from the YAML as the **sample working root**. Expected structure:

```
${sample.path}/
  mgf/        # MGF files for the de novo stage
  mzML/       # mzML files for database search
  Denovo/     # De novo engines' working dir and merged outputs live here
```

- De novo stage reads `mgf/*.mgf`.
- Database-search stage reads `mzML/*.mzML`.

---

## Pipeline stages & outputs

### 1) De novo (`--stage denovo`)
- Runs Casanovo / PepNet / InstaNovo (enable per config).
- Selects high-confidence subsequences by per-residue scores and merges them.
- **Output:** merged FASTA  
  `${sample.path}/Denovo/Denovo_TE_SoftMerged_InstaNovo.fasta`  
  (each entry header is the peptide sequence itself).

### 2) Database build & search (`--stage dbsearch`)
- Builds per-engine FASTAs:
  - **Human (PE=0)** + **TE (PE=1)** + **contaminants (PE=0)** [+ **decoys**].
  - TE subset is preferably derived from *de novo tags → BLAST short vs a large TE FASTA → matched entries*.
- **Outputs:** per-engine DBs, search results under:  
  `${sample.path}/DB_search_iProphet/{COMET,MSFRAGGER,MSGFPLUS}/...`

### 3) Post-processing (`--stage postprocess`)
A typical chain:
1. TPP RefreshParser
2. `xinteract` to iProphet (**decoy prefix fixed to `rev_`**)
3. `pepxml2csv` table export
4. Charge-state–aware FDR filtering
5. **Short-sequence homology filtering**: remove peptides that have **100% identity** full-length hits against the canonical human proteome (BLAST short).
6. Optional HLA binding affinity prediction using **MixMHCpred**.

**Outputs include (examples):**
- iProphet pepXML/TSV: `${sample.path}/DB_search_iProphet/IPROPHET/`
- FDR-filtered results, homology-filtered tables, and integrated binding predictions; final integrated CSV often named like `final_TE_results.csv`

---

## YAML configuration (field quick reference)

Below is an example `config.yaml`. Fill in absolute paths for your environment.

```yaml
sample:
  name: SAMPLE_01
  path: /abs/path/to/SAMPLE_01

gpu:
  id: "0"

denovo:
  min_score: 0.75
  min_length: 8
  casanovo:
    enable: true
    binary: /path/to/casanovo
    model_ckpt: /path/to/model.ckpt
    config_yaml: /path/to/casanovo.yaml
  pepnet:
    enable: false
    python_bin: /path/to/conda/envs/pepnet/bin/python
    script: /path/to/pepnet_denovo.py
    model_h5: /path/to/model.h5
  instanovo:
    enable: true
    python_bin: /path/to/python
    module_call: instanovo.main
    model_ckpt: /path/to/instanovo.ckpt

database_build:
  human_fasta: /path/to/uniprot_human.fasta
  contaminants_fasta: /path/to/contaminants.fasta
  te_db_fasta: /path/to/TE_big.fasta
  add_decoy: true
  te_selection:
    method: blast
    chunks: 15
    threads_per_chunk: 4

blast:
  blastp_bin: /path/to/blastp

search_engines:
  comet:
    enable: true
    binary: /path/to/comet
    params_file: /path/to/comet.params
  msfragger:
    enable: true
    java_bin: /usr/bin/java
    jar: /path/to/MSFragger.jar
    base_params: /path/to/fragger.params
  msgfplus:
    enable: true
    java_bin: /usr/bin/java
    jar: /path/to/MSGFPlus.jar
    params_file: /path/to/msgfplus.params

postprocessing:
  tpp:
    refreshparser: /path/to/RefreshParser
    xinteract: /path/to/xinteract
    pepxml2csv: /path/to/pepxml2csv.py
  fdr_threshold: 0.03
  blastp:
    enable: true
    binary: /path/to/blastp
    human_proteome_blastdb: /path/to/blastdb/human
    threads_per_job: 4

hla_binding:
  enable: true
  mixmhcpred_bin: /path/to/MixMHCpred
  hla_alleles: "HLA-A02:01,HLA-B07:02"
  peptide_length_min: 8
  peptide_length_max: 14
```

> The code expects the decoy prefix to be **`rev_`** in iProphet. Ensure your decoy generation and TPP setup match this.

---

## FAQ

- **Post-processing says “no pepXML found”.**  
  Likely your search stage did not produce pepXML (e.g., only MSGF+ `.mzid`). Convert to pepXML first, then re-run post-processing.

- **iProphet complains about decoys.**  
  The pipeline assumes the decoy prefix is `rev_`. Make sure your decoy generation and TPP flags are aligned.

- **Parallel BLAST seems to “hang”.**  
  Polling relies on `ps`/`grep` on Unix-like systems. Please run on Linux with common GNU userland available.

---

## Versions & acknowledgements

Please document the tested versions of **Comet / MSFragger / MS-GF+ / TPP (or Philosopher) / MixMHCpred / NCBI BLAST / ProteoWizard** you used so others can reproduce your results.
