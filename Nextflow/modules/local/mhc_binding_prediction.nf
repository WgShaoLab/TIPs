process MHC_BINDING_PREDICTION {
  tag "${sample_id}"
  container params.container_dbsearch
  publishDir "${params.outdir}/${sample_id}", mode: 'copy'

  input:
  tuple val(sample_id), path(filtered_peptides)

  output:
  tuple val(sample_id), path("final_TE_results.tsv"), emit: final_results
  tuple val(sample_id), path("binding.fasta"), emit: binding_fasta, optional: true
  tuple val(sample_id), path("binding.result"), emit: binding_result, optional: true

  script:
  """
  set -euo pipefail

 
  if [ ! -s "${filtered_peptides}" ]; then
    echo "Empty input file"
    cp "${filtered_peptides}" final_TE_results.tsv
    exit 0
  fi

 
  header=\$(head -n1 "${filtered_peptides}")
  peptide_col=\$(echo "\$header" | tr '\\t' '\\n' | grep -n -i "^peptide\$" | cut -d: -f1 | head -n1)

  if [ -z "\$peptide_col" ]; then
    echo "No Peptide column found"
    cp "${filtered_peptides}" final_TE_results.tsv
    exit 0
  fi

  echo "Peptide column index: \$peptide_col"


  tail -n+2 "${filtered_peptides}" | \\
    awk -F'\\t' -v col="\$peptide_col" '{
      pep = \$col
      gsub(/^[ \\t]+|[ \\t]+\$/, "", pep)  # trim whitespace
      len = length(pep)
      if (len >= 8 && len <= 14) {
        print ">" pep
        print pep
      }
    }' > binding.fasta


  peptide_count=\$(grep -c "^>" binding.fasta || true)
  echo "Valid peptides for MHC prediction: \$peptide_count"

  if [ "\$peptide_count" -eq 0 ]; then
    echo "No valid peptides (length 8-14)"
    cp "${filtered_peptides}" final_TE_results.tsv
    exit 0
  fi


  echo "Running MixMHCpred with alleles: ${params.mhc_alleles}"
  
  if ! command -v MixMHCpred &> /dev/null; then
    echo "ERROR: MixMHCpred not found"
    cp "${filtered_peptides}" final_TE_results.tsv
    exit 1
  fi

  MixMHCpred -i binding.fasta -o binding.result -a ${params.mhc_alleles}

  if [ ! -s binding.result ]; then
    echo "MixMHCpred produced no output"
    cp "${filtered_peptides}" final_TE_results.tsv
    exit 0
  fi



  grep -v "^#" binding.result | tail -n+2 > binding_data.tmp


  result_header=\$(grep -v "^#" binding.result | head -n1)
  
  rank_col=\$(echo "\$result_header" | tr '\\t' '\\n' | grep -n "%Rank_bestAllele" | cut -d: -f1)
  allele_col=\$(echo "\$result_header" | tr '\\t' '\\n' | grep -n "BestAllele" | cut -d: -f1)
  pep_result_col=\$(echo "\$result_header" | tr '\\t' '\\n' | grep -n "^Peptide\$" | cut -d: -f1)

  echo "Result columns - Peptide: \$pep_result_col, Rank: \$rank_col, Allele: \$allele_col"

  if [ -z "\$rank_col" ] || [ -z "\$allele_col" ] || [ -z "\$pep_result_col" ]; then
    echo "Cannot find required columns in MixMHCpred output"
    cp "${filtered_peptides}" final_TE_results.tsv
    exit 0
  fi


  awk -F'\\t' -v pcol="\$pep_result_col" -v rcol="\$rank_col" -v acol="\$allele_col" \\
    'NR>0 { print \$pcol "\\t" \$rcol "\\t" \$acol }' binding_data.tmp > lookup.tmp



  echo -e "\${header}\\t%Rank_bestAllele\\tBestAllele" > final_TE_results.tsv


  tail -n+2 "${filtered_peptides}" | while IFS=\$'\\t' read -r line; do

    pep=\$(echo "\$line" | cut -f"\$peptide_col" | tr -d ' ')
    

    match=\$(grep -P "^\${pep}\\t" lookup.tmp || true)
    
    if [ -n "\$match" ]; then
      rank=\$(echo "\$match" | cut -f2)
      allele=\$(echo "\$match" | cut -f3)
    else
      rank="NA"
      allele="NA"
    fi
    
    echo -e "\${line}\\t\${rank}\\t\${allele}"
  done >> final_TE_results.tsv


  total=\$(tail -n+2 final_TE_results.tsv | wc -l)
  with_pred=\$(tail -n+2 final_TE_results.tsv | awk -F'\\t' '\$NF != "NA" && \$(NF-1) != "NA"' | wc -l)
  echo "Final results: \$with_pred/\$total peptides with MHC predictions"


  rm -f binding_data.tmp lookup.tmp

  echo "MHC binding prediction completed successfully"
  """
}