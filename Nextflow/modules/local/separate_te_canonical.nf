nextflow.enable.dsl=2

process SEPARATE_TE_CANONICAL {
  tag "${sample_id}_${engine}"

  input:
  tuple val(sample_id), val(engine), path(pepxml_files)

  output:
  tuple val(sample_id), val(engine), path("*_TE.pep.xml"), emit: te_pepxml
  tuple val(sample_id), val(engine), path("*_Canonical.pep.xml"), emit: canonical_pepxml

  script:
  """
  python << 'PYTHON_EOF'
import sys
sys.path.append('${params.code_base}')
import os
import subprocess
from lxml import etree

def get_TE_canonical_pepxml_fixed(path, engine='${engine}'):

    if engine in ['MSfragger', 'MS_GF']:
        file_list = [os.path.join(path, x) for x in os.listdir(path)
                     if '.pepXML' in x and '_TE' not in x and '_Canonical' not in x]
        for i in file_list:
            new_name = i.replace('.pepXML', '.pep.xml')
            if i != new_name and os.path.exists(i):
                subprocess.run(['mv', i, new_name], check=True)

    file_list = [x for x in os.listdir(path)
                 if x.endswith('.pep.xml')
                 and '_TE.pep.xml' not in x
                 and '_Canonical.pep.xml' not in x]

    if not file_list:
        raise RuntimeError("No pepxml files to process!")

    type_list = ['TE', 'Canonical']

    for file in file_list:
        file_path = os.path.join(path, file)

        for type_ in type_list:
            tree_copy = etree.parse(file_path)
            root_copy = tree_copy.getroot()

            for child in root_copy:
                spectrum_queries_to_remove = []

                for spectrum_query in child:
                    if 'spectrum_query' not in spectrum_query.tag:
                        continue

                    if len(spectrum_query) == 0:
                        spectrum_queries_to_remove.append(spectrum_query)
                        continue

                    search_result = spectrum_query[0]
                    search_hits_to_remove = []

                    for search_hit in search_result:
                        protein = search_hit.get("protein", "")
                        is_te = 'Denovo' in protein

                        if type_ == 'TE' and not is_te:
                            search_hits_to_remove.append(search_hit)
                        elif type_ == 'Canonical' and is_te:
                            search_hits_to_remove.append(search_hit)

                    for hit in search_hits_to_remove:
                        search_result.remove(hit)

                    if len(search_result) == 0:
                        spectrum_queries_to_remove.append(spectrum_query)

                for query in spectrum_queries_to_remove:
                    child.remove(query)

            output = file_path.replace('.pep.xml', f'_{type_}.pep.xml')
            tree_copy.write(output, pretty_print=True, xml_declaration=True, encoding="UTF-8")

get_TE_canonical_pepxml_fixed('.')
PYTHON_EOF
  """
}