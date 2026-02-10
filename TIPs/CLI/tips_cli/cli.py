from __future__ import annotations

from pathlib import Path

import typer
from tips_cli.steps.integrate_te import IntegrateTeOutputs, integrate_te_step
from tips_cli.steps.denovo_te import DenovoTeOutputs, denovo_te_step
from tips_cli.steps.prepare import PrepareOutputs, prepare_step
from tips_cli.steps.denovo import DenovoOutputs, denovo_step_batch_script
from tips_cli.steps.search_te import SearchTeOutputs, search_te_step


app = typer.Typer(help="TIPs dockerized CLI")

def _safe_symlink(src: Path, dst: Path) -> None:
    """Create a symlink if dst does not exist; overwrite only if dst is a broken link."""
    dst.parent.mkdir(parents=True, exist_ok=True)

    # If dst exists as a real path (or a valid symlink), keep it.
    if dst.exists():
        return

    # If dst is a broken symlink, remove it and recreate.
    if dst.is_symlink():
        try:
            if not dst.resolve().exists():
                dst.unlink()
        except Exception:
            # If resolve fails for any reason, treat it as broken and recreate safely.
            dst.unlink()

    dst.symlink_to(src)



def _pick_first_existing(*candidates: Path) -> Path | None:
    """Return the first existing candidate, otherwise None."""
    for p in candidates:
        if p.exists():
            return p
    return None


def _materialize_resources(sample_path: Path, resources_path: Path) -> dict:
    """Symlink required resources into sample_path/resources and return resolved relative paths."""
    resources_path = resources_path.resolve()
    sample_path = sample_path.resolve()

    # Required resource roots
    db_root = resources_path / "db"
    te_root = resources_path / "te"

    # Only search_params is supported; resources/params has been deprecated and removed.
    search_params_root = resources_path / "search_params"

    if not db_root.exists():
        raise FileNotFoundError(f"Missing resources/db under: {resources_path}")
    if not te_root.exists():
        raise FileNotFoundError(f"Missing resources/te under: {resources_path}")
    if not search_params_root.exists():
        raise FileNotFoundError(f"Missing resources/search_params under: {resources_path}")

    # Canonical DB
    human_fasta = _pick_first_existing(
        db_root / "UP000005640_9606_processed.fasta",
        db_root / "UP000005640_9606_processed_tag_0.fasta",
    )
    contaminants_fasta = _pick_first_existing(
        db_root / "Contaminants.fasta",
        db_root / "Contaminants_tag_0.fasta",
    )
    if human_fasta is None:
        raise FileNotFoundError(f"Missing human fasta under: {db_root}")
    if contaminants_fasta is None:
        raise FileNotFoundError(f"Missing contaminants fasta under: {db_root}")

    _safe_symlink(human_fasta, sample_path / "resources" / "db" / human_fasta.name)
    _safe_symlink(contaminants_fasta, sample_path / "resources" / "db" / contaminants_fasta.name)

    # Optional canonical BLAST DB (directory symlink)
    human_blastdb_dir = db_root / "human_blastpdb"
    if human_blastdb_dir.exists():
        _safe_symlink(human_blastdb_dir, sample_path / "resources" / "db" / "human_blastpdb")

    # TE resources
    te_protein_npy = _pick_first_existing(
        te_root / "hg38_rmsk_6_frame_singleLine_done.npy",
        te_root / "blastdb" / "hg38_rmsk_6_frame_singleLine_done.npy",
    )
    te_class_npy = _pick_first_existing(te_root / "TE_class_dic.npy")
    if te_protein_npy is None:
        raise FileNotFoundError(f"Missing TE protein npy under: {te_root}")
    if te_class_npy is None:
        raise FileNotFoundError(f"Missing TE class npy under: {te_root}")

    _safe_symlink(te_protein_npy, sample_path / "resources" / "te" / te_protein_npy.name)
    _safe_symlink(te_class_npy, sample_path / "resources" / "te" / te_class_npy.name)

    te_blastdb_dir = te_root / "blastdb"
    if not te_blastdb_dir.exists():
        raise FileNotFoundError(f"Missing TE blastdb directory under: {te_root}")
    _safe_symlink(te_blastdb_dir, sample_path / "resources" / "te" / "blastdb")

    # Search engine params
    comet_params = search_params_root / "comet.params"
    msfragger_params = search_params_root / "msfragger.params"
    msgfplus_params = search_params_root / "msgfplus.params"
    for p in [comet_params, msfragger_params, msgfplus_params]:
        if not p.exists():
            raise FileNotFoundError(f"Missing search params: {p}")

    _safe_symlink(comet_params, sample_path / "resources" / "search_params" / "comet.params")
    _safe_symlink(msfragger_params, sample_path / "resources" / "search_params" / "msfragger.params")
    _safe_symlink(msgfplus_params, sample_path / "resources" / "search_params" / "msgfplus.params")

    # Guess TE blast db prefix
    te_blast_prefix = _pick_first_existing(
        te_blastdb_dir / "hg38_rmsk_6_frame_singleLine_done.fasta",
        te_blastdb_dir / "hg38_rmsk_6_frame_singleLine_done",
    )
    te_blast_prefix_str = "resources/te/blastdb/hg38_rmsk_6_frame_singleLine_done.fasta"
    if te_blast_prefix is not None:
        te_blast_prefix_str = f"resources/te/blastdb/{te_blast_prefix.name}"

    # Guess canonical blastdb prefix
    human_blast_prefix_str = ""
    if human_blastdb_dir.exists():
        known = human_blastdb_dir / "human_proteome_UP000005640_noERV_LINE"
        if (known.with_suffix(".pin")).exists() or (known.with_suffix(".psq")).exists() or (known.with_suffix(".phr")).exists():
            human_blast_prefix_str = "resources/db/human_blastpdb/human_proteome_UP000005640_noERV_LINE"

    return {
        "human_fasta": f"resources/db/{human_fasta.name}",
        "contaminants_fasta": f"resources/db/{contaminants_fasta.name}",
        "human_blastdb_prefix": human_blast_prefix_str,
        "te_protein_npy": f"resources/te/{te_protein_npy.name}",
        "te_class_npy": f"resources/te/{te_class_npy.name}",
        "te_blast_db_prefix": te_blast_prefix_str,
        "comet_params": "resources/search_params/comet.params",
        "msfragger_params": "resources/search_params/msfragger.params",
        "msgfplus_params": "resources/search_params/msgfplus.params",
    }


@app.command()
def init(
    sample_path: Path = typer.Argument(..., help="Create TIPs folder layout in sample_path."),
    resources: Path | None = typer.Option(
        None,
        help="Optional resources folder. If provided, symlink resources into sample_path/resources and auto-generate config/denovo_te.yaml.",
    ),
    overwrite_config: bool = typer.Option(False, help="Overwrite config/denovo_te.yaml if it already exists."),
):
    sample_path = sample_path.resolve()

    # Base TIPs layout
    for d in ["raw", "mzML", "mgf", "logs", "tmp", "denovo", "db", "search", "integrate", "config", "resources"]:
        (sample_path / d).mkdir(parents=True, exist_ok=True)

    # denovo-te layout
    (sample_path / "denovo" / "tags").mkdir(parents=True, exist_ok=True)
    (sample_path / "denovo" / "te").mkdir(parents=True, exist_ok=True)

    # resources layout
    (sample_path / "resources" / "te").mkdir(parents=True, exist_ok=True)
    (sample_path / "resources" / "db").mkdir(parents=True, exist_ok=True)
    (sample_path / "resources" / "search_params").mkdir(parents=True, exist_ok=True)

    resolved = {}
    if resources is not None:
        resolved = _materialize_resources(sample_path=sample_path, resources_path=resources)

    cfg_path = sample_path / "config" / "denovo_te.yaml"
    if (not cfg_path.exists()) or overwrite_config:
        cfg_text = f"""te:
  protein_npy: {resolved.get('te_protein_npy', 'resources/te/hg38_rmsk_6_frame_singleLine_done.npy')}
  class_npy: {resolved.get('te_class_npy', 'resources/te/TE_class_dic.npy')}
  # BLAST DB prefix (prefix of .pin/.psq/.phr). Do NOT include the index suffixes.
  # The prefix may include ".fasta" if that is part of the -out name used by makeblastdb.
  blast_db_prefix: {resolved.get('te_blast_db_prefix', 'resources/te/blastdb/hg38_rmsk_6_frame_singleLine_done.fasta')}

db:
  # Canonical fasta paths under sample_path.
  human_fasta: {resolved.get('human_fasta', 'resources/db/UP000005640_9606_processed.fasta')!r}
  contaminants_fasta: {resolved.get('contaminants_fasta', 'resources/db/Contaminants.fasta')!r}
  # If empty, search-te may build BLAST DB for canonical removal (see search-te behavior).
  # integrate-te requires an existing BLAST DB prefix with .pin/.psq/.phr present.
  human_blastdb_prefix: {resolved.get('human_blastdb_prefix', 'resources/db/human_blastpdb/human_proteome_UP000005640_noERV_LINE')!r}

search_params:
  comet: {resolved.get('comet_params', 'resources/search_params/comet.params')!r}
  msfragger: {resolved.get('msfragger_params', 'resources/search_params/msfragger.params')!r}
  msgfplus: {resolved.get('msgfplus_params', 'resources/search_params/msgfplus.params')!r}

blastp:
  bin: blastp
  task: blastp-short
  evalue: 20000
  min_pident: 75
  outfmt: "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"

tag_extract:
  min_score: 0.75
  min_length: 8
  strip_non_aa: true

tools:
  instanovo:
    glob: "*.result.csv"
  casanovo:
    glob: "*_result.mztab"
  pepnet:
    glob: "*"
"""
        cfg_path.write_text(cfg_text, encoding="utf-8")

    typer.echo(f"Initialized: {sample_path}")
    if resources is not None:
        typer.echo(f"Resources linked from: {resources.resolve()}")
        typer.echo(f"Config: {cfg_path}")


@app.command()
def prepare(
    sample_path: Path = typer.Argument(..., help="TIPs sample folder (sample_path)."),
    image: str = typer.Option("tips-core:0.2", help="Docker image for core tools."),
    formats: str = typer.Option("mzml,mgf", help="Comma-separated: mzml,mgf"),
    stage: str = typer.Option("move", help="Stage RAW into raw/: move|copy|none"),
    jobs: int = typer.Option(1, help="Parallel jobs for RAW conversion."),
    dry_run: bool = typer.Option(False, help="Write commands but do not execute."),
    resume: bool = typer.Option(False, help="Skip outputs that already exist."),
):
    sample_path = sample_path.resolve()

    fmt_set = {x.strip().lower() for x in formats.split(",") if x.strip()}
    do_mzml = "mzml" in fmt_set
    do_mgf = "mgf" in fmt_set
    if not (do_mzml or do_mgf):
        raise typer.BadParameter("formats must include at least one of: mzml,mgf")

    outputs = PrepareOutputs(
        sample_path=sample_path,
        mzml_dir=sample_path / "mzML",
        mgf_dir=sample_path / "mgf",
        raw_dir=sample_path / "raw",
        log_dir=sample_path / "logs",
        commands_sh=sample_path / "logs" / "prepare_commands.sh",
        log_file=sample_path / "logs" / "prepare.log",
    )

    prepare_step(
        outputs=outputs,
        image=image,
        do_mzml=do_mzml,
        do_mgf=do_mgf,
        stage_mode=stage.lower(),
        jobs=jobs,
        dry_run=dry_run,
        resume=resume,
    )

    typer.echo(f"Done. Commands: {outputs.commands_sh}")
    typer.echo(f"Log: {outputs.log_file}")


@app.command()
def denovo(
    sample_path: Path = typer.Argument(..., help="TIPs sample folder (sample_path)."),
    tool: str = typer.Option("instanovo", help="De novo tool: instanovo|casanovo|pepnet|all"),
    param_file: Path | None = typer.Option(None, help="Optional parameter file (currently used by Casanovo as --config)."),
    script: Path | None = typer.Option(None, help="Override the built-in denovo batch script path."),
    sample: str | None = typer.Option(
        None,
        help="Optional sample prefix for output filenames. If not set, defaults to sample_path folder name.",
    ),
    gpu: str = typer.Option("auto", help="GPU selection: auto | -1 | 0,1,2,..."),
    dry_run: bool = typer.Option(False, help="Write commands but do not execute."),
    resume: bool = typer.Option(False, help="Skip if expected outputs already exist."),
):
    from tips_cli.steps.denovo import default_denovo_script_path

    script_path = script.resolve() if script is not None else default_denovo_script_path()
    sample_path = sample_path.resolve()
    mgf_dir = sample_path / "mgf"

    # The batch script requires --sample; default to the sample folder name if not provided.

    sample_id = sample.strip() if (sample is not None and sample.strip()) else sample_path.name

    outputs = DenovoOutputs(
        sample_path=sample_path,
        mgf_dir=mgf_dir,
        out_dir=sample_path / "denovo",
        logs_dir=sample_path / "logs",
        commands_sh=(sample_path / "logs" / f"denovo_{tool}_commands.sh"),
        log_file=(sample_path / "logs" / f"denovo_{tool}.log"),
    )

    denovo_step_batch_script(
        outputs=outputs,
        script_path=script_path,
        tool=tool,
        sample=sample_id,
        gpu=gpu,
        param_file=param_file,
        dry_run=dry_run,
        resume=resume,
    )

    typer.echo(f"Done. Commands: {outputs.commands_sh}")
    typer.echo(f"Log: {outputs.log_file}")


@app.command()
def denovo_te(
    sample_path: Path = typer.Argument(..., help="TIPs sample folder (sample_path)."),
    image: str = typer.Option("tips-core:0.2", help="Docker image for core tools."),
    config: Path | None = typer.Option(None, help="denovo-te config yaml. Default: sample_path/config/denovo_te.yaml"),
    tools: str = typer.Option("instanovo,casanovo,pepnet", help="Comma-separated: instanovo,casanovo,pepnet"),
    blastp_workers: int = typer.Option(15, help="Parallel blastp tasks (split fasta count)."),
    threads_per_task: int = typer.Option(4, help="blastp -num_threads per task."),
    dry_run: bool = typer.Option(False, help="Write commands but do not execute."),
    resume: bool = typer.Option(False, help="Skip if outputs already exist."),
    overwrite: bool = typer.Option(False, help="Overwrite outputs if they exist."),
    keep_intermediate: bool = typer.Option(False, help="Keep intermediate blastp folder."),
):
    sample_path = sample_path.resolve()
    if config is None:
        config = sample_path / "config" / "denovo_te.yaml"
    config = config.resolve()

    if not config.exists():
        raise typer.BadParameter(f"Config not found: {config}. Run `tips init {sample_path}` to create a template.")

    tool_list = [x.strip() for x in tools.split(",") if x.strip()]
    if not tool_list:
        raise typer.BadParameter("tools must include at least one of: instanovo,casanovo,pepnet")

    outputs = DenovoTeOutputs(
        sample_path=sample_path,
        denovo_dir=sample_path / "denovo",
        tags_dir=sample_path / "denovo" / "tags",
        te_dir=sample_path / "denovo" / "te",
        log_file=sample_path / "logs" / "denovo_te.log",
        commands_sh=sample_path / "logs" / "denovo_te_commands.sh",
    )

    denovo_te_step(
        outputs=outputs,
        image=image,
        config_path=config,
        tools=tool_list,
        blastp_workers=blastp_workers,
        threads_per_task=threads_per_task,
        dry_run=dry_run,
        resume=resume,
        overwrite=overwrite,
        keep_intermediate=keep_intermediate,
    )

    typer.echo(f"Done. Commands: {outputs.commands_sh}")
    typer.echo(f"Log: {outputs.log_file}")
    typer.echo(f"TE fasta: {outputs.te_dir / 'Denovo_TE_SoftMerged_Merged.fasta'}")


@app.command(name="search-te")
def search_te(
    sample_path: Path = typer.Argument(..., help="TIPs sample folder (sample_path)."),
    image: str = typer.Option("tips-core:0.2", help="Docker image for core tools."),
    config: Path | None = typer.Option(None, help="Config yaml. Default: sample_path/config/denovo_te.yaml"),
    engines: str = typer.Option("comet,msfragger,msgfplus", help="Comma-separated: comet,msfragger,msgfplus"),
    mzml_glob: str = typer.Option("*.mzML", help="Glob under sample_path/mzML/ to select inputs."),
    fdr: float = typer.Option(0.03, help="TE-only PSM FDR threshold."),
    blastp_threads: int = typer.Option(20, help="blastp -num_threads for canonical removal."),
    dry_run: bool = typer.Option(False, help="Write commands but do not execute."),
    resume: bool = typer.Option(False, help="Skip if outputs already exist."),
    overwrite: bool = typer.Option(False, help="Overwrite outputs if they exist."),
):
    sample_path = sample_path.resolve()
    if config is None:
        config = sample_path / "config" / "denovo_te.yaml"
    config = config.resolve()

    if not config.exists():
        raise typer.BadParameter(f"Config not found: {config}. Run `tips init {sample_path}` to create a template.")

    engine_list = [x.strip().lower() for x in engines.split(",") if x.strip()]
    if not engine_list:
        raise typer.BadParameter("engines must include at least one of: comet,msfragger,msgfplus")

    outputs = SearchTeOutputs(
        sample_path=sample_path,
        search_dir=sample_path / "search",
        mzml_dir=sample_path / "mzML",
        log_file=sample_path / "logs" / "search_te.log",
        commands_sh=sample_path / "logs" / "search_te_commands.sh",
    )

    search_te_step(
        outputs=outputs,
        image=image,
        config_path=config,
        engines=engine_list,
        mzml_glob=mzml_glob,
        fdr=fdr,
        blastp_threads=blastp_threads,
        dry_run=dry_run,
        resume=resume,
        overwrite=overwrite,
    )

    typer.echo(f"Done. Commands: {outputs.commands_sh}")
    typer.echo(f"Log: {outputs.log_file}")
    typer.echo(f"Search dir: {outputs.search_dir}")


@app.command(name="integrate-te")
def integrate_te(
    sample_path: Path = typer.Argument(..., help="TIPs sample folder (sample_path)."),
    image: str = typer.Option("tips-core:0.2", help="Docker image for core tools."),
    fdr: float = typer.Option(0.03, help="Philosopher filter threshold for both --pep and --psm."),
    threads: int = typer.Option(20, help="Threads for philosopher iprophet."),
    blastp_db_prefix: str = typer.Option(
        "resources/db/human_blastpdb/human_proteome_UP000005640_noERV_LINE",
        help="BLASTP database prefix under sample_path (no extension).",
    ),
    blastp_threads: int = typer.Option(20, help="blastp -num_threads for canonical removal."),
    hla_alleles: str | None = typer.Option(
        None,
        help="Optional MixMHCpred alleles string, e.g. A0201,A6801,B1302",
    ),
    dry_run: bool = typer.Option(False, help="Write commands but do not execute."),
    resume: bool = typer.Option(False, help="Skip if outputs already exist."),
    overwrite: bool = typer.Option(False, help="Overwrite outputs if they exist."),
):
    sample_path = sample_path.resolve()

    outputs = IntegrateTeOutputs(
        sample_path=sample_path,
        integrate_dir=sample_path / "integrate",
        log_file=sample_path / "logs" / "integrate_te.log",
        commands_sh=sample_path / "logs" / "integrate_te_commands.sh",
    )

    # Default inputs from search-te outputs (TE-only)
    pepxml_comet = sample_path / "search" / "comet" / "te" / "comet_merged_TE.pep.xml"
    pepxml_msfragger = sample_path / "search" / "msfragger" / "te" / "msfragger_merged_TE.pep.xml"
    pepxml_msgfplus = sample_path / "search" / "msgfplus" / "te" / "msgfplus_merged_TE.pep.xml"
    combined_fasta = sample_path / "search" / "db" / "combined_target_decoy.fasta"

    for p in [pepxml_comet, pepxml_msfragger, pepxml_msgfplus, combined_fasta]:
        if not p.exists():
            raise typer.BadParameter(f"Missing required input: {p}")

    integrate_te_step(
        outputs=outputs,
        image=image,
        combined_target_decoy_fasta=combined_fasta,
        pepxml_comet=pepxml_comet,
        pepxml_msfragger=pepxml_msfragger,
        pepxml_msgfplus=pepxml_msgfplus,
        fdr=fdr,
        threads=threads,
        blastp_db_prefix=blastp_db_prefix,
        blastp_threads=blastp_threads,
        hla_alleles=hla_alleles,
        dry_run=dry_run,
        resume=resume,
        overwrite=overwrite,
    )

    typer.echo(f"Done. Commands: {outputs.commands_sh}")
    typer.echo(f"Log: {outputs.log_file}")
    typer.echo(f"Final: {sample_path / 'integrate' / 'final_TE_results.tsv'}")

@app.command(name="run-all")
def run_all(
    sample_path: Path = typer.Argument(
        Path("."),
        help="TIPs sample folder. Default: current working directory.",
    ),
    resources: Path = typer.Option(
        ...,
        help="Resources folder to symlink into sample_path/resources (same as `tips init --resources`).",
    ),
    image: str = typer.Option("tips-core:0.2", help="Docker image for core tools."),
    overwrite_config: bool = typer.Option(False, help="Overwrite config/denovo_te.yaml if it already exists."),
    # prepare
    prepare_formats: str = typer.Option("mzml,mgf", help="Prepare formats. Comma-separated: mzml,mgf"),
    prepare_stage: str = typer.Option("move", help="Stage RAW into raw/: move|copy|none"),
    prepare_jobs: int = typer.Option(1, help="Parallel jobs for RAW conversion."),
    # denovo
    denovo_tool: str = typer.Option("all", help="De novo tool: instanovo|casanovo|pepnet|all"),
    denovo_param_file: Path | None = typer.Option(None, help="Optional parameter file (currently used by Casanovo as --config)."),
    denovo_script: Path | None = typer.Option(None, help="Override the built-in denovo batch script path."),
    denovo_sample: str | None = typer.Option(
        None,
        help="Optional sample prefix for de novo outputs. If not set, defaults to sample_path folder name.",
    ),
    denovo_gpu: str = typer.Option("auto", help="GPU selection: auto | -1 | 0,1,2,..."),
    # denovo-te
    denovo_te_tools: str | None = typer.Option(
        None,
        help="Comma-separated tools for denovo-te. Default: derived from --denovo-tool.",
    ),
    denovo_te_blastp_workers: int = typer.Option(15, help="Parallel blastp tasks (split fasta count)."),
    denovo_te_threads_per_task: int = typer.Option(4, help="blastp -num_threads per task."),
    denovo_te_keep_intermediate: bool = typer.Option(False, help="Keep intermediate blastp folder."),
    # search-te
    search_engines: str = typer.Option("comet,msfragger,msgfplus", help="Comma-separated: comet,msfragger,msgfplus"),
    search_mzml_glob: str = typer.Option("*.mzML", help="Glob under sample_path/mzML/ to select inputs."),
    search_fdr: float = typer.Option(0.03, help="TE-only PSM FDR threshold."),
    search_blastp_threads: int = typer.Option(20, help="blastp -num_threads for canonical removal."),
    # integrate-te
    integrate_fdr: float = typer.Option(0.03, help="Philosopher filter threshold for both --pep and --psm."),
    integrate_threads: int = typer.Option(20, help="Threads for philosopher iprophet."),
    integrate_blastp_db_prefix: str = typer.Option(
        "resources/db/human_blastpdb/human_proteome_UP000005640_noERV_LINE",
        help="BLASTP database prefix under sample_path.",
    ),
    integrate_blastp_threads: int = typer.Option(20, help="blastp -num_threads for canonical removal."),
    integrate_hla_alleles: str | None = typer.Option(
        None,
        help="Optional MixMHCpred alleles string, e.g. A0201,A6801,B1302",
    ),
    # common flags
    dry_run: bool = typer.Option(False, help="Write commands but do not execute."),
    resume: bool = typer.Option(False, help="Skip if expected outputs already exist."),
    overwrite: bool = typer.Option(False, help="Overwrite outputs if they exist (TE/search/integrate steps)."),
):
    """Run TIPs end-to-end from init -> prepare -> denovo -> denovo-te -> search-te -> integrate-te."""
    sample_path = sample_path.resolve()
    resources = resources.resolve()

    # 1) init (folder layout + config)
    init(sample_path=sample_path, resources=resources, overwrite_config=overwrite_config)

    # 2) prepare (RAW -> mzML/mgf)
    prepare(
        sample_path=sample_path,
        image=image,
        formats=prepare_formats,
        stage=prepare_stage,
        jobs=prepare_jobs,
        dry_run=dry_run,
        resume=resume,
    )

    # 3) denovo
    denovo(
        sample_path=sample_path,
        tool=denovo_tool,
        param_file=denovo_param_file,
        script=denovo_script,
        sample=denovo_sample,
        gpu=denovo_gpu,
        dry_run=dry_run,
        resume=resume,
    )

    # 4) denovo-te
    cfg = sample_path / "config" / "denovo_te.yaml"

    tool_str = denovo_te_tools
    if tool_str is None or not tool_str.strip():
        t = denovo_tool.strip().lower()
        if t == "all":
            tool_str = "instanovo,casanovo,pepnet"
        else:
            tool_str = t

    denovo_te(
        sample_path=sample_path,
        image=image,
        config=cfg,
        tools=tool_str,
        blastp_workers=denovo_te_blastp_workers,
        threads_per_task=denovo_te_threads_per_task,
        dry_run=dry_run,
        resume=resume,
        overwrite=overwrite,
        keep_intermediate=denovo_te_keep_intermediate,
    )

    # 5) search-te
    search_te(
        sample_path=sample_path,
        image=image,
        config=cfg,
        engines=search_engines,
        mzml_glob=search_mzml_glob,
        fdr=search_fdr,
        blastp_threads=search_blastp_threads,
        dry_run=dry_run,
        resume=resume,
        overwrite=overwrite,
    )

    # 6) integrate-te
    integrate_te(
        sample_path=sample_path,
        image=image,
        fdr=integrate_fdr,
        threads=integrate_threads,
        blastp_db_prefix=integrate_blastp_db_prefix,
        blastp_threads=integrate_blastp_threads,
        hla_alleles=integrate_hla_alleles,
        dry_run=dry_run,
        resume=resume,
        overwrite=overwrite,
    )

if __name__ == "__main__":
    app()
