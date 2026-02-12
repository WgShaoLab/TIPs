## beta Nextflow Workflow

### Setting

**Before running Nextflow, ensure that the definite path of the reference databases and code in nextflow.config are set correctly**

```bash
// Reference databases and files
te_blastdb_dir   = "/path/to/TIPs/resources/te/blastdb"  
te_npy           = "/path/to/TIPs/resources/te/hg38_rmsk_6_frame_singleLine_done.npy"
te_class_dic     = "/path/to/TIPs/resources/te/TE_class_dic.npy"
uniprot_fasta    = "/path/to/TIPs/resources/db/UP000005640_9606_processed.fasta"
contaminants_fasta = "/path/to/TIPs/resources/db/Contaminants.fasta"
blastp_db        = "/path/to/TIPs/resources/db/human_blastpdb"
// Search engine configs
comet_params     = "/path/to/TIPs/resources/search_params/comet.params"
msfragger_params = "/path/to/TIPs/resources/search_params/msfragger.params"
msgf_params      = "/path/to/TIPs/resources/search_params/msgfplus.params"
// Host code base mounted by container options 
code_base = "/path/to/TIPs/Nextflow/code"
// Host de novo runner script
denovo_runner = "/path/to/TIPs/Nextflow/code/run_denovo_tool_batch_final.py"
```

### Usage

Assume you have:

- Docker images already available
- `nextflow.config`: A config file modified according to user definition
- `sample_dir`: A sample directory(Directory containing raw/*.raw)

- `out_dir`:Specify the directory location for all outputs (including intermediate files)

```bash
nextflow run main.nf -c nextflow.config --sample_dir <DIR> --outdir <DIR>
```

### Output

The output of this workflow contains:

- MzML and mgf format files converted from raw file format in **mzML** and **mgf** 
- De novo Results in **denovo**
- Fasta file composed of sequence tags extracted by blastp in **denovo_db**
- Fasta file used for searching databases in **search_db**
- Database search results in **dbsearch**
- MHC binding filtered final results in **integrate**

```bash
.
├── blastp_denovo/
│   └── merged_raw_peptide_Child_fasta/
├── dbsearch/
│   ├── Comet/
│   ├── fdr/
│   │   └── TE/
│   ├── iprophet/
│   │   └── TE/
│   ├── MSfragger/
│   ├── MS_GF/
│   ├── peptideprophet/
│   │   └── TE/
├── denovo/
│   ├── casanovo/
│   ├── instanovo/
│   ├── merged_raw_peptide.fasta
│   └── pepnet/
├── denovo_db/
├── integrate/
│   └── final_TE_results.tsv
├── mgf/
├── mzML/
├── post/
└── search_db/

```

