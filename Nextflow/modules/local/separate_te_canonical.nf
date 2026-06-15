nextflow.enable.dsl=2

process SEPARATE_TE_CANONICAL {
  tag "${sample_id}_${engine}"

  input:
  tuple val(sample_id), val(engine), path(pepxml_files), path(te_fasta)

  output:
  tuple val(sample_id), val(engine), path("*_TE.pep.xml"), emit: te_pepxml
  tuple val(sample_id), val(engine), path("*_Canonical.pep.xml"), emit: canonical_pepxml

  script:
  '''
  python3 << 'PYTHON_EOF'
import os
import subprocess
from lxml import etree

def load_te_headers(te_fasta_path):
    """Load TE FASTA headers for exact matching (mirrors CLI search_te logic)."""
    headers = set()
    if not te_fasta_path or not os.path.exists(te_fasta_path):
        return headers
    with open(te_fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                h = line[1:].strip()
                if not h:
                    continue
                headers.add(h)
                headers.add(h.split()[0])
    return headers

def is_te_hit(protein, te_headers):
    """Determine if a search_hit protein belongs to TE (mirrors CLI search_te logic)."""
    if not protein:
        return False

    raw = str(protein)
    p0 = raw.split("#")[0].strip()
    p1 = p0.split()[0] if p0 else ""

    # First try exact matching against TE FASTA headers
    if te_headers and ((p0 in te_headers) or (p1 in te_headers)):
        return True

    # Fallback: check for 'Denovo' token
    return ("Denovo" in raw)


# Load TE headers
te_fasta_path = "!{te_fasta}"
te_headers = load_te_headers(te_fasta_path)
print(f"Loaded {len(te_headers)} TE FASTA headers")

# Normalize pepXML extensions (.pepXML -> .pep.xml)
for fname in os.listdir('.'):
    if fname.endswith('.pepXML'):
        new_name = fname.replace('.pepXML', '.pep.xml')
        if fname != new_name:
            os.rename(fname, new_name)

# Process each pep.xml file
engine = "!{engine}"

file_list = [x for x in os.listdir('.')
             if x.endswith('.pep.xml')
             and '_TE.pep.xml' not in x
             and '_Canonical.pep.xml' not in x]

if not file_list:
    raise RuntimeError("No pepxml files to process!")

type_list = ['TE', 'Canonical']

for file in file_list:
    file_path = os.path.join('.', file)

    for type_ in type_list:
        tree_copy = etree.parse(file_path)
        root_copy = tree_copy.getroot()

        for child in root_copy:
            spectrum_queries_to_remove = []
            for spectrum_query in child:
                if 'spectrum_query' not in spectrum_query.tag:
                    continue

                if len(spectrum_query) == 0:
                    continue

                search_result = spectrum_query[0]
                search_hits_to_remove = []

                for search_hit in search_result:
                    protein = search_hit.get("protein", "")
                    is_te = is_te_hit(protein, te_headers)

                    if type_ == 'TE' and not is_te:
                        search_hits_to_remove.append(search_hit)
                    elif type_ == 'Canonical' and is_te:
                        search_hits_to_remove.append(search_hit)

                for hit in search_hits_to_remove:
                    search_result.remove(hit)

                # Remove spectrum_query if no search_hits remain after filtering
                if len(search_result) == 0:
                    spectrum_queries_to_remove.append(spectrum_query)

            for sq in spectrum_queries_to_remove:
                child.remove(sq)

        output = file_path.replace('.pep.xml', f'_{engine}_{type_}.pep.xml')
        tree_copy.write(output, pretty_print=True, xml_declaration=True, encoding="UTF-8")
        print(f"Wrote: {output}")

print("TE/Canonical separation complete.")
PYTHON_EOF
  '''
}