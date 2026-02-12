nextflow.enable.dsl=2

process BLASTP_DENOVO_PEPTIDES {
  tag "${sample_id}"
  
  publishDir "${params.outdir}/blastp_denovo", mode: 'copy'

  input:
  tuple val(sample_id), path(peptides_fasta)
  path te_fasta

  output:
  tuple val(sample_id), path("*_Child_fasta"), emit: blastp_dir

  script:
  def num_task = params.blastp_workers ?: 15
  def basename = peptides_fasta.baseName
  """
  #!/bin/bash
  set -e

  CHILD_DIR="${basename}_Child_fasta"
  mkdir -p "\${CHILD_DIR}"

  if [ -f "${te_fasta}.phr" ] || [ -f "${te_fasta}.pin" ]; then
    echo "Found existing BLAST database index files"
    BLAST_DB="${te_fasta}"
  else
    echo "Building new database..."
    makeblastdb -in ${te_fasta} -dbtype prot -out te_db
    BLAST_DB="te_db"
  fi

  echo "Using BLAST database: \${BLAST_DB}"

  cat > run_blastp.py << 'PYEOF'
import os
import sys
import subprocess
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed

def read_fasta(fasta_file):
    sequences = []
    header = None
    seq = []
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    sequences.append((header, ''.join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header is not None:
            sequences.append((header, ''.join(seq)))
    return sequences

def split_fasta(sequences, output_dir, num_files):
    total = len(sequences)
    chunk_size = max(1, total // num_files)
    file_list = []
    
    for i in range(num_files):
        start = i * chunk_size
        if i == num_files - 1:
            end = total
        else:
            end = start + chunk_size
        
        if start >= total:
            break
            
        chunk = sequences[start:end]
        if not chunk:
            continue
            
        output_file = os.path.join(output_dir, "sub_fasta_{}.fasta".format(i))
        with open(output_file, 'w') as f:
            for header, seq in chunk:
                f.write(">" + header + chr(10) + seq + chr(10))
        file_list.append(output_file)
    
    return file_list

def run_blastp(fasta_file, db_path):
    output_file = fasta_file.replace('.fasta', '.txt')
    cmd = [
        'blastp',
        '-task', 'blastp-short',
        '-query', fasta_file,
        '-db', db_path,
        '-out', output_file,
        '-outfmt', '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq',
        '-evalue', '20000',
        '-num_threads', '4'
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return output_file
    except subprocess.CalledProcessError as e:
        print("BLASTP error for {}: {}".format(fasta_file, e.stderr))
        return None

def filter_and_process_results(input_file, output_file, final_output):
    column_names = ['qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 
                    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qseq', 'sseq']
    
    if os.path.getsize(input_file) == 0:
        open(output_file, 'w').close()
        with open(final_output, 'w') as f:
            writer = csv.writer(f, delimiter=chr(9))
            writer.writerow(column_names + ['q_length', 'qseq_I2L', 'sseq_I2L', 'I/L_identity'])
        return
    
    filtered_rows = []
    final_rows = []
    
    with open(input_file, 'r') as f:
        reader = csv.reader(f, delimiter=chr(9))
        for row in reader:
            if len(row) < 14:
                continue
            try:
                qaccver = row[0]
                pident = float(row[2])
                qstart = int(row[6])
                gapopen = int(row[5])
                length = int(row[3])
                qseq = row[12]
                sseq = row[13]
                
                if qstart == 1 and gapopen == 0 and pident >= 75:
                    filtered_rows.append(row)
                    
                    q_length = len(qaccver)
                    if length == q_length and len(qseq) == len(sseq) and gapopen == 0:
                        qseq_i2l = qseq.replace('I', 'L')
                        sseq_i2l = sseq.replace('I', 'L')
                        il_identity = str(qseq_i2l == sseq_i2l)
                        final_rows.append(row + [str(q_length), qseq_i2l, sseq_i2l, il_identity])
            except (ValueError, IndexError) as e:
                continue
    
    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter=chr(9))
        for row in filtered_rows:
            writer.writerow(row)
    
    with open(final_output, 'w') as f:
        writer = csv.writer(f, delimiter=chr(9))
        writer.writerow(column_names + ['q_length', 'qseq_I2L', 'sseq_I2L', 'I/L_identity'])
        for row in final_rows:
            writer.writerow(row)
    
    print("Filtered results: {} rows".format(len(filtered_rows)))
    print("Final filtered results: {} rows".format(len(final_rows)))

def main(peptides_fasta, child_dir, num_task, blast_db):
    print("Reading fasta: {}".format(peptides_fasta))
    sequences = read_fasta(peptides_fasta)
    print("Total sequences: {}".format(len(sequences)))

    if len(sequences) == 0:
        print("Warning: No sequences found!")
        open(os.path.join(child_dir, "merged_blastp.txt"), 'w').close()
        open(os.path.join(child_dir, "merged_blastp_filter.txt"), 'w').close()
        with open(os.path.join(child_dir, "merged_blastp_TE.csv"), 'w') as f:
            f.write(chr(9).join(['qaccver', 'saccver', 'pident', 'length', 'mismatch', 
                                 'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                 'evalue', 'bitscore', 'qseq', 'sseq', 'q_length',
                                 'qseq_I2L', 'sseq_I2L', 'I/L_identity']) + chr(10))
        return

    print("Splitting fasta into {} files...".format(num_task))
    fasta_files = split_fasta(sequences, child_dir, num_task)
    print("Created {} sub-fasta files".format(len(fasta_files)))

    print("Running BLASTP...")
    result_files = []
    with ProcessPoolExecutor(max_workers=min(num_task, len(fasta_files))) as executor:
        futures = {executor.submit(run_blastp, f, blast_db): f for f in fasta_files}
        for future in as_completed(futures):
            result = future.result()
            if result:
                result_files.append(result)
                print("Completed: {}".format(result))

    print("Merging results...")
    merged_file = os.path.join(child_dir, "merged_blastp.txt")
    with open(merged_file, 'w') as outfile:
        for rf in result_files:
            if os.path.exists(rf):
                with open(rf, 'r') as infile:
                    outfile.write(infile.read())

    filtered_file = os.path.join(child_dir, "merged_blastp_filter.txt")
    final_file = os.path.join(child_dir, "merged_blastp_TE.csv")
    
    print("Filtering and processing...")
    filter_and_process_results(merged_file, filtered_file, final_file)

    print("BLASTP analysis completed!")

if __name__ == "__main__":
    peptides_fasta = sys.argv[1]
    child_dir = sys.argv[2]
    num_task = int(sys.argv[3])
    blast_db = sys.argv[4]
    main(peptides_fasta, child_dir, num_task, blast_db)
PYEOF

  python3 run_blastp.py "${peptides_fasta}" "\${CHILD_DIR}" "${num_task}" "\${BLAST_DB}"
  """
}