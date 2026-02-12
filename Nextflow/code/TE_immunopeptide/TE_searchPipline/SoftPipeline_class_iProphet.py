from random import sample
import pandas as pd
import os
import numpy as np
import subprocess
from io import StringIO
from Bio import SeqIO
import time,re,os
import multiprocessing
from lxml import etree
import sys
import ast
import functools

from TE_immunopeptide.TE_searchPipline.Class.fasta_class import Fasta

def generate_decoy(fasta_file):
    # 在原fasta的基础上生成decoy序列
    f = open(fasta_file, 'r').read()
    record_list = f.split('>')[1:]
    record_dic = {}
    for record in record_list:
        record_dic[record.split('\n')[0]] = record.split('\n')[1]
    decoy_fasta = fasta_file.split('.')[0]+'_decoy.fasta'
    with open(decoy_fasta, 'a') as out_file:
        for header,seq in record_dic.items():
            out_file.write(f'>rev_{header}\n{seq[::-1]}\n')

def extract_subsequences_sub(sequence, scores, min_score=0.75, min_length=8):
    result = []
    n = len(sequence)
    i = 0

    while i < n:
        if scores[i] > min_score:
            start = i
            while i < n and scores[i] > min_score:
                i += 1
            end = i

            # 判断子序列长度是否满足要求
            if end - start >= min_length:
                # 提取初始子序列
                subseq = sequence[start:end]
                subseq_scores = scores[start:end]
                result.append(subseq)

                # 对子序列进一步细分以找到最高平均分数的子序列
                best_avg_score = -1
                best_subseq = ""

                for j in range(len(subseq) - min_length + 1):
                    for k in range(j + min_length, len(subseq) + 1):
                        temp_subseq = subseq[j:k]
                        temp_scores = subseq_scores[j:k]
                        avg_score = sum(temp_scores) / len(temp_scores)

                        if avg_score > best_avg_score:
                            best_avg_score = avg_score
                            best_subseq = temp_subseq

                result.append(best_subseq)
        else:
            i += 1

    return result

def timing_decorator(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # args[0] 是类实例本身（self）
        self = args[0]
        class_name = self.__class__.__name__
        if hasattr(self, 'Soft_name'):  # 检查是否有 soft_name 属性
            print(f"{class_name}: Executing {func.__name__} for soft_name: {self.Soft_name}")

        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()

        print(f"\n{class_name}: {func.__name__} executed in {end_time - start_time:.4f} seconds\n")
        return result

    return wrapper

def get_TE_canonical_pepxmlComet(path):
    if 'MSfragger' in path or 'MS_GF' in path:
        file_list = [os.path.join(path,x) for x in os.listdir(path) if '.pepXML' in x]
        for i in file_list:
            subprocess.run(f'mv {i} {i.replace(".pepXML",".pep.xml")}', shell=True)
    type_list = ['TE','Canonical']
    file_list = [x for x in os.listdir(path) if '.pep.xml' in x]
    for type in type_list:
        for file in file_list:
            file_path = os.path.join(path,file)
            output = file_path.replace('.pep.xml',f'_{type}.pep.xml')
            tree = etree.parse(file_path)
            root = tree.getroot()
            # 遍历所有子元素
            for child in root:
                # 遍历子元素的子元素
                for spectrum_query in child:
                    if 'spectrum_query' in spectrum_query.tag:
                        # print(spectrum_query.tag, spectrum_query.attrib, spectrum_query.text)  # 打印子标签名、属性和文本内容
                        search_result = spectrum_query[0]
                        for search_hit in search_result:
                            protein = search_hit.get("protein")
                            if type=='TE' and 'Denovo' not in protein:
                                # 删除不符合条件的search_hit节点
                                parent = search_hit.getparent()
                                parent.remove(search_hit)
                            elif type=='Canonical' and 'Denovo' in protein:
                                # 删除不符合条件的search_hit节点
                                parent = search_hit.getparent()
                                parent.remove(search_hit)
            # TE 将修改后的 XML 保存到新文件中
            tree.write(output, pretty_print=True, xml_declaration=True, encoding="UTF-8")

def filter_by_fdr(peptideProphet_resultPath, fdr_threshold):
    """
    根据设定的 FDR 阈值筛选质谱数据。
    """
    # peptideProphet_resultPath = f'/data/48/wuqian/Noncanonical/mutation/NCI_3784Mel/mzML/test/merged.pep.csv'
    # peptideProphet_resultPath = f'/data/48/wuqian/TE_DBsearch/0D5P_all/DB_search_iProphet/Comet/Comet_merged_Canonical.pep.csv'
    outdf_path = peptideProphet_resultPath.replace('.pep.csv','_fdr.csv')
    df_merged = pd.DataFrame()
    try:
        df_all = pd.read_csv(peptideProphet_resultPath, sep='\t')
    except:
        print('there are no records in the csv')
        df_merged.to_csv(outdf_path, index=False, sep='\t')
        return

    if df_all.__len__() < 20:
        df_merged.to_csv(outdf_path, index=False, sep='\t')
        return

    df_all['protein_primary'] = df_all['protein'].apply(
        lambda x: x.split('#')[0]
    )
    df_all['peptideProphet'] = df_all['analysis_result'].apply(ast.literal_eval)
    df_all['peptideProphet_probability'] = df_all['peptideProphet'].apply(
        lambda x: float(x[0]['peptideprophet_result']['probability'])
    )
    df_all['decoy'] = df_all['protein_primary'].apply(
        lambda x: 1 if 'rev_' in x else 0
    )
    df_all = df_all.sort_values(by=['peptideProphet_probability'], ascending=False).reset_index(drop=True)
    for charge, df in df_all.groupby('assumed_charge'):
        df = df_all[df_all['assumed_charge'] == charge]
        # 累计统计 target 和 decoy 数量
        df['target_cumsum'] = (df['decoy'] == 0).cumsum()
        df['decoy_cumsum'] = (df['decoy'] == 1).cumsum()

        # 计算 FDR = Decoy / Target
        df['fdr'] = df['decoy_cumsum'] / df['target_cumsum']
        # 筛选出满足 FDR 阈值的行
        df.reset_index(drop=True, inplace=True)
        reversed_df = df[::-1]
        if reversed_df[reversed_df['fdr'] < fdr_threshold].__len__() == 0:
            print(f'Charge {charge} has no results with FDR < {fdr_threshold}.')
            continue
        first_index = reversed_df[reversed_df['fdr'] < fdr_threshold].index[0]
        filtered_df = df.iloc[:first_index + 1]
        filtered_df = filtered_df[filtered_df['decoy'] == 0]
        print(f'Charge {charge} has {filtered_df.__len__()} results with FDR < {fdr_threshold}.')
        if df_merged.empty:
            df_merged = filtered_df
        else:
            df_merged = pd.concat([df_merged, filtered_df])

    # 只返回满足 FDR 阈值的数据
    if df_merged.empty:
        df_merged.to_csv(outdf_path, index=False, sep='\t')
    else:
        df_merged = df_merged.sort_values(by=['peptideProphet_probability'], ascending=False)
        df_merged.to_csv(outdf_path, index=False,sep='\t')

def Blastp_remove(csv_Path):
    # remove canonical peptide for iProphet or PeptideProphet
    df = pd.read_csv(csv_Path, sep='\t')
    assert 'peptide' in df.columns or 'Peptide' in df.columns, \
        "Error : the input df must have 'peptide' or 'Peptide' column"
    if 'peptide' in df.columns:
        peptide_column = 'peptide'
    else:
        peptide_column = 'Peptide'

    blastp_fasta_path = csv_Path.split('.')[0] + '_blastp.fasta'
    with open(blastp_fasta_path, 'a') as f:
        peptide_list = list(set(df[peptide_column].tolist()))
        for seq  in peptide_list:
            f.write(f'>{seq}\n{seq}\n')

    # run blastp
    blastp_result_path = csv_Path.split('.')[0] + '_blastpResult.txt'
    blastp_order = f'/home/qwu/anaconda3/envs/blast/bin/blastp -task blastp-short \
        -query {blastp_fasta_path} \
        -db /data/48/wuqian/genome/human/hg38/blastp_human_UP000005640_noERV_LINE/human_proteome_UP000005640_noERV_LINE \
        -out {blastp_result_path} \
        -outfmt 6 -evalue 20000 -num_threads 20 &'
    subprocess.run(blastp_order, shell=True)

    try:
        while int(subprocess.check_output("ps aux | grep -v grep | grep -c 'blastp'", shell=True)) > 0:
            time.sleep(5)
    except:
        print('Blastp Done')

    # process the result of blastp
    column_names = [
        'Identifier',
        'Protein',
        'Score',
        'Length',
        'Start',
        'End',
        'MatchStart',
        'MatchEnd',
        'DatabaseStart',
        'DatabaseEnd',
        'Value1',
        'Value2'
    ]
    blastp_df = pd.read_csv(f'{blastp_result_path}', sep='\t',
                            names=column_names)
    blastp_df['seq_len'] = blastp_df['Identifier'].apply(
        lambda x: len(x)
    )
    blastp_df = blastp_df[(blastp_df['MatchStart'] == 1) & (blastp_df['MatchEnd'] == blastp_df['seq_len'])]
    blastp_df = blastp_df[blastp_df['Score'] == 100]
    remove_peptide_list = blastp_df['Identifier'].drop_duplicates().to_list()

    # filter the df and output
    df_filter = df[~df[peptide_column].isin(remove_peptide_list)]
    df_filter_path = csv_Path.split('.')[0] + '_blastp.filter.txt'

    df_filter.to_csv(df_filter_path, index=False, sep='\t')
    print('Blastp Canonical peptide Filter Done')

def predict_bindingAffinity_iprophet(sample_path,geneType):
    # sample_path = '/data/48/wuqian/fast/Denovo_test'
    # geneType = 'A0201,A6801,B1302,B4001,C0304,C0602'
    peptide_file = os.path.join(sample_path,'DB_search_iProphet','peptide_blastp.filter.txt')
    if not os.path.exists(peptide_file):
        print(f'{sample_path} is not completed')
        return 1

    df = pd.read_csv(f'{sample_path}/DB_search/soft_merged.txt',sep='\t')
    binding_fasta = os.path.join(sample_path,'DB_search_iProphet','binding.fasta')
    with open(binding_fasta,'a') as f:
        for idx,row in df.iterrows():
            seq = row['Peptide']
            if seq.__len__() < 15:
                f.write(f'>{seq}\n{seq}\n')

    # predict using MixMHCpred
    binding_result = binding_fasta.replace("fasta","result")
    predict_sh = f'/home/qwu/soft/MixMHCpred-master/MixMHCpred \
    -i {binding_fasta} \
    -o {binding_result} \
    -a {geneType}'
    subprocess.call(predict_sh, shell=True)

    # merged the result df
    result_df = pd.read_csv(binding_result,sep='\t',comment='#')
    result_df = result_df[['Peptide','%Rank_bestAllele','BestAllele']]
    rank_dic  = dict(zip(result_df['Peptide'],result_df['%Rank_bestAllele']))
    Alle_dic = dict(zip(result_df['Peptide'],result_df['BestAllele']))
    df['%Rank_bestAllele'] = df['Peptide'].apply(
        lambda x: rank_dic[x] if x in rank_dic.keys() else np.nan
    )
    df['BestAllele'] = df['Peptide'].apply(
        lambda x: Alle_dic[x] if x in Alle_dic.keys() else ''
    )

    df.to_csv(f'{sample_path}/DB_search_iProphet/final_TE_results.csv',sep='\t',index=False)

class DenovoSoftPipeline:
    def __init__(self, sample_path, soft_name,GPU_id = '0'):
        # create the Casanovo work directory
        self.sample_path = sample_path
        self.Soft_name = soft_name
        if not os.path.exists(os.path.join(self.sample_path,'Denovo')):
            subprocess.run(f"mkdir {os.path.join(self.sample_path,'Denovo')}", shell=True)
        self.Soft_workpath = os.path.join(self.sample_path, 'Denovo', self.Soft_name)
        subprocess.run(f'mkdir {self.Soft_workpath}', shell=True)

        self.mgf_list = os.listdir(f'{sample_path}/mgf')
        self.output_fasta_path = os.path.join(self.Soft_workpath, 'merged_raw_peptide.fasta')
        self.GPU_id = GPU_id
        os.environ["CUDA_VISIBLE_DEVICES"] = self.GPU_id  # change the GPU using by casanovo

    def runOrder(self):
        pass

    def Result_process(self, min_score=0.75, min_length=8):
        return pd.DataFrame()

    def WirteInFasta(self, min_score=0.75, min_length=8):
        # 处理结果，输出成fasta，方便后续多个软件的合并
        filter_df = self.Result_process(min_score, min_length)
        with open(self.output_fasta_path, 'a') as f:
            peptide_list = []
            for idx, row in filter_df.iterrows():
                seq_list = row['Filtered_Subsequences']
                for seq in seq_list:
                    peptide_list.append(seq)
            peptide_list = list(set(peptide_list))
            for seq in peptide_list:
                f.write(f'>{seq}\n{seq}\n')

        return 0

    def run(self):
        self.runOrder()
        self.WirteInFasta()

class CasaNovo_searchPipline(DenovoSoftPipeline):
    def __init__(self, sample_path, soft_name,GPU_id = '0'):
        # create the Casanovo workdirectory
        super().__init__(sample_path, soft_name,GPU_id)

    @timing_decorator
    def runOrder(self):
        for mgf in self.mgf_list:
            mgf_file_path = os.path.join(self.sample_path, 'mgf', mgf)
            casanovo_order = f'/home/qwu/anaconda3/envs/casanovo_env/bin/casanovo sequence \
            {mgf_file_path} -m /home/qwu/soft/casanovo/casanovo_v4_2_0.ckpt \
            -o {self.Soft_workpath}/{mgf.split(".")[0]}_result.txt \
            --config /home/qwu/soft/casanovo/casanovo.yaml'
            subprocess.run(casanovo_order, shell=True)

    @timing_decorator
    def Result_process(self, min_score=0.75, min_length=8):
        df = pd.DataFrame()
        path = self.Soft_workpath
        file_list = [x for x in os.listdir(path) if 'mztab' in x]
        for file in file_list:
            with open(os.path.join(path, file), 'r') as file:
                lines = file.readlines()
            # 过滤掉以 "MTD" 开头的行
            filtered_lines = [line for line in lines if not line.startswith('MTD')]
            # 将过滤后的数据转换为 DataFrame
            # 使用 StringIO 将过滤后的文本作为文件来处理
            filtered_data = StringIO(''.join(filtered_lines))
            # 读取为 pandas DataFrame
            df_single = pd.read_csv(filtered_data, sep='\t')
            if df.empty:
                df = df_single
            else:
                df = pd.concat([df, df_single])

        filter_df = df[['sequence', 'opt_ms_run[1]_aa_scores']]
        filter_df['sequence'] = filter_df['sequence'].apply(lambda x: re.sub('[^A-Z]', '', x))
        filter_df = filter_df.rename(columns={'opt_ms_run[1]_aa_scores': 'Score',
                                              'sequence': 'Sequence', })

        filter_df["Score"] = filter_df["Score"].apply(lambda x: list(map(float, x.split(','))))

        # 对每一行数据应用函数
        filter_df["Filtered_Subsequences"] = filter_df.apply(
            lambda row: extract_subsequences_sub(row["Sequence"], row["Score"], min_score=min_score,
                                                 min_length=min_length), axis=1)
        filter_df = filter_df[filter_df["Filtered_Subsequences"].str.len() > 0]

        return filter_df

class Pepent_searchPipline(DenovoSoftPipeline):
    def __init__(self, sample_path, soft_name,GPU_id = '0'):
        # create the Casanovo workdirectory
        super().__init__(sample_path, soft_name,GPU_id)

    @timing_decorator
    def runOrder(self):
        # 例如，设置 CUDA 和 cuDNN 的环境变量
        env = os.environ.copy()
        env['CUDA_HOME'] = '/usr/local/cuda'
        env['LD_LIBRARY_PATH'] = '/usr/local/cuda/lib64:' + env.get('LD_LIBRARY_PATH', '')
        env['PATH'] = '/usr/local/cuda/bin:' + env['PATH']
        for mgf in self.mgf_list:
            mzML_file_path = os.path.join(self.sample_path, 'mgf', mgf)
            Pepnet_order = f'/home/qwu/anaconda3/envs/pepnet/bin/python \
            /home/qwu/soft/PepNet/PepNet-master/denovo.py \
            --input {mzML_file_path} \
            --model /home/qwu/soft/PepNet/PepNet-master/model.h5 \
            --output {self.Soft_workpath}/{mgf.split(".")[0]}.result'
            with open(f'{self.Soft_workpath}/Pepnet_order.sh', 'w') as f:
                f.write(Pepnet_order)
            order = f'sh {self.Soft_workpath}/Pepnet_order.sh'
            subprocess.run(order, shell=True, executable='/bin/bash')

    @timing_decorator
    def Result_process(self, min_score=0.75, min_length=8):
        # 处理结果，输出成fasta
        file_list = [x for x in os.listdir(self.Soft_workpath) if 'result' in x]
        df = pd.DataFrame()
        for file in file_list:
            if df.empty:
                df = pd.read_csv(os.path.join(self.Soft_workpath, file), sep='\t')
            else:
                df = pd.concat([df, pd.read_csv(os.path.join(self.Soft_workpath, file), sep='\t')])
        filter_df = df[['DENOVO', 'Positional Score']]
        filter_df = filter_df.rename(columns={'Positional Score': 'Score', 'DENOVO': 'Sequence'})
        filter_df['Score'] = filter_df['Score'].apply(lambda x: re.sub('[\[\]\s]', '', x))
        filter_df["Score"] = filter_df["Score"].apply(lambda x: list(map(float, [x for x in x.split(',') if x != ''])))
        filter_df = filter_df.dropna(axis=0, how='any')
        filter_df['Sequence'] = filter_df['Sequence'].apply(lambda x: re.sub('[^A-Z]', '', x))

        # 对每一行数据应用函数
        filter_df["Filtered_Subsequences"] = filter_df.apply(
            lambda row: extract_subsequences_sub(row["Sequence"], row["Score"], min_score=min_score,
                                                 min_length=min_length), axis=1)
        filter_df = filter_df[filter_df["Filtered_Subsequences"].str.len() > 0]

        return filter_df

class Instanovo_searchPipline(DenovoSoftPipeline):
    def __init__(self, sample_path, soft_name ,GPU_id = '0'):
        # create the Casanovo workdirectory
        super().__init__(sample_path, soft_name,GPU_id)


    @timing_decorator
    def runOrder(self):
        Instanovo_order = f'export PYTHONPATH=/home/qwu/soft/InstaNovo/InstaNovo:$PYTHONPATH && /home/qwu/anaconda3/envs/instanovo2_5090/bin/python \
        -m instanovo.transformer.predict data_path={os.path.join(self.sample_path,"mgf")}/*.mgf \
        model_path=/home/qwu/soft/InstaNovo/InstaNovo/instanovo_extended.ckpt\
        output_path={self.Soft_workpath}/Instanovo.result denovo=True'
        subprocess.run(Instanovo_order, shell=True)

    @timing_decorator
    def Result_process(self, min_score=0.75, min_length=8):
        # 处理结果，输出成fasta
        df = pd.read_csv(os.path.join(self.Soft_workpath, 'Instanovo.result'))
        df = df.rename(columns={'preds': 'DENOVO', 'token_log_probs': 'Positional Score_raw'})
        filter_df = df[['DENOVO', 'Positional Score_raw']]
        filter_df['Positional Score'] = filter_df['Positional Score_raw'].apply(
            lambda x: [float(np.exp(float(val))) for val in x[1:-1].split(',')])
        filter_df = filter_df[['DENOVO', 'Positional Score']]
        filter_df = filter_df.rename(columns={'Positional Score': 'Score', 'DENOVO': 'Sequence'})
        filter_df['Sequence'] = filter_df['Sequence'].apply(lambda x: re.sub('[^A-Z]', '', x))

        # 对每一行数据应用函数
        filter_df["Filtered_Subsequences"] = filter_df.apply(
            lambda row: extract_subsequences_sub(row["Sequence"], row["Score"], min_score=min_score,
                                                 min_length=min_length), axis=1)
        filter_df = filter_df[filter_df["Filtered_Subsequences"].str.len() > 0]

        return filter_df

class DB_SoftPipline:
    def __init__(self, sample_path : str, soft_name : str):
        self.sample_path = sample_path
        self.Soft_name = soft_name
        if not os.path.exists(os.path.join(self.sample_path,'DB_search_iProphet')):
            subprocess.run(f"mkdir {os.path.join(self.sample_path,'DB_search_iProphet')}", shell=True)
        self.Soft_workpath = os.path.join(self.sample_path, 'DB_search_iProphet', self.Soft_name)
        if not os.path.exists(self.Soft_workpath):
            subprocess.run(f'mkdir {self.Soft_workpath}', shell=True)
            self.mzML_path = f'{self.sample_path}/mzML'
            self.mzML_list = os.listdir(f'{self.sample_path}/mzML')
            subprocess.run(f'ln -s {self.sample_path}/mzML/*.mzML {self.Soft_workpath}', shell=True)

        self.fasta_file = os.path.join(self.Soft_workpath,f'{self.Soft_name}_DBsearch.fasta')
        self.TE_Canonical_list = ['TE','Canonical']

    def is_running(self,soft_path) -> bool:
        try:
            # 如果返回的进程数大于0，则说明进程还在运行
            return int(subprocess.check_output(
                f"ps aux | grep -v grep | grep -c '{soft_path}'",
                shell=True)) > 0
        except Exception:
            # 出现异常时认为进程已结束
            return False

    @timing_decorator
    def DB_fasta(self):
        # generate fasta file
        if not os.path.exists(f'{self.sample_path}/Denovo/Denovo_result_merged_tag_1.fasta'):
            sample_Denovo_fasta = f'{self.sample_path}/Denovo/Denovo_TE_SoftMerged_InstaNovo.fasta'
            if not os.path.exists(sample_Denovo_fasta):
                raise FileNotFoundError(f"The file '{sample_Denovo_fasta}' does not exist.")
            fasta = Fasta(sample_Denovo_fasta, uniform_PE_tag=1)
            fasta.output(f'{self.sample_path}/Denovo/Denovo_result_merged_tag_1.fasta')

        human_fasta = '/data/48/wuqian/TE_DBsearch/public_fasta/UP000005640_9606_processed_tag_0.fasta'
        sample_Denovo_fasta = f'{self.sample_path}/Denovo/Denovo_result_merged_tag_1.fasta'
        contaminants_fasta = f'/data/48/wuqian/TE_DBsearch/public_fasta/Contaminants_tag_0.fasta'

        Target_fasta = f'cat {human_fasta} {sample_Denovo_fasta} {contaminants_fasta}\
               > {self.Soft_workpath}/{self.Soft_name}.fasta'
        subprocess.run(Target_fasta, shell=True)

        # 生成MSGF使用的database
        generate_decoy(f'{self.Soft_workpath}/{self.Soft_name}.fasta')

        merged_fasta_order = (f'cat {self.Soft_workpath}/{self.Soft_name}.fasta {self.Soft_workpath}/{self.Soft_name}_decoy.fasta'
                                f' > {self.fasta_file}')
        subprocess.run(merged_fasta_order, shell=True)

    @timing_decorator
    def runSoft(self):
        pass

    @timing_decorator
    def PeptideProphet(self):
        # sperate the TE PSM and the Canonical PSM
        get_TE_canonical_pepxmlComet(self.Soft_workpath)

        # process the PSM
        fasta_DB = os.path.join(self.Soft_workpath,f'{self.Soft_name}_DBsearch.fasta')
        for type in self.TE_Canonical_list:
            order = f'/home/qwu/soft/TPP_6.0.0/tpp/bin/xinteract -OAPNd -PPM -eN -p0 -THREADS=16 -drev_ \
            -D{fasta_DB} \
            -N{self.Soft_workpath}/{self.Soft_name}_merged_{type}.pep.xml \
            {self.Soft_workpath}/*_{type}.pep.xml'
            subprocess.run(order, shell=True)

    @timing_decorator
    def pepXML_fdr(self,fdr=0.03):
        for type in self.TE_Canonical_list:
            type_pepXML = f'{self.Soft_workpath}/{self.Soft_name}_merged_{type}.pep.xml'
            assert os.path.exists(type_pepXML), f"The file '{type_pepXML}' does not exist."

            subprocess.run(f'python /home/qwu/soft/TPP_6.0.0/tpp/bin/pepxml2csv.py {type_pepXML}', shell=True)
            print(f'{type} pepXML has transform to csv')

            # filter using fdr
            type_csv = type_pepXML.replace('.xml', '.csv')
            filter_by_fdr(type_csv,fdr)

    @timing_decorator
    def blastp_filter(self):
        TE_csv = f'{self.Soft_workpath}/{self.Soft_name}_TE_fdr.csv'
        Blastp_remove(TE_csv)

    @timing_decorator
    def runPipeline(self):
        pass

class Comet_Pipeline(DB_SoftPipline):
    def __init__(self,sample_path):
        super().__init__(sample_path,'Comet')
        # generate fasta file
        if not os.path.exists(self.fasta_file):
            self.DB_fasta()

    @timing_decorator
    def runSoft(self):
        comet_order = f'/home/qwu/soft/TPP_6.0.0/tpp/bin/comet {self.Soft_workpath}/*.mzML \
         -P/data/48/wuqian/TE_DBsearch/config/comet_iprophet.params -D{self.fasta_file}'
        subprocess.run(comet_order, shell=True)

        while self.is_running('/home/qwu/soft/TPP_6.0.0/tpp/bin/comet'):
            time.sleep(5)

        print('raw Comet Done')

    @timing_decorator
    def runPipeline(self):

        self.runSoft()
        self.PeptideProphet()
        self.pepXML_fdr()

class MSfragger_Pipeline(DB_SoftPipline):
    def __init__(self,sample_path):
        super().__init__(sample_path,'MSfragger')
        self.MSfragger_config = f'{self.Soft_workpath}/MSfragger_Params_iprophet.params'
        if not os.path.exists(self.fasta_file):
            self.DB_fasta()

    @timing_decorator
    def generateConfig(self):
        # generate config file
        temp_config_file = open(r'/data/48/wuqian/TE_DBsearch/config/closed_fragger.params', 'r').read()
        temp_config_file = temp_config_file.replace('database_name =', f'database_name = {self.fasta_file}')
        with open(f'{self.MSfragger_config}', 'w') as f:
            f.write(temp_config_file)

    @timing_decorator
    def runSoft(self):
        MSfragger_order = f'java -Xmx150g -jar /home/qwu/soft/fragpipe/tools/MSFragger-4.1/MSFragger-4.1.jar \
        {self.MSfragger_config} {self.Soft_workpath}/*.mzML '
        subprocess.run(MSfragger_order, shell=True)

        print('raw Comet Done')

    @timing_decorator
    def runPipeline(self):
        # generate fasta file
        self.generateConfig()

        self.runSoft()
        self.PeptideProphet()
        self.pepXML_fdr()

class MS_GF_Pipeline(DB_SoftPipline):
    def __init__(self,sample_path):
        super().__init__(sample_path,'MS_GF')
        self.MS_GF_config = '/data/48/wuqian/TE_DBsearch/config/MSGFPlus_Params_iprophet.txt'
        if not os.path.exists(self.fasta_file):
            self.DB_fasta()

    @timing_decorator
    def runSoft(self,mzML_path):
        MSfragger_order = f'java -jar -Xmx150g /home/qwu/soft/MS-GF+/MSGFPlus.jar -s {mzML_path} \
        -conf "{self.MS_GF_config}" -d "{self.fasta_file}" -o "{mzML_path.split(".mzML")[0]}.mzid" '
        subprocess.run(MSfragger_order, shell=True)

    @timing_decorator
    def mzidToPepXML(self):
        # modify the mzid file using script
        mzid_list = [os.path.join(self.Soft_workpath, x) for x in os.listdir(self.Soft_workpath) if x.endswith('.mzid')]
        for mzid_file in mzid_list:
            modify_order = f'/home/qwu/soft/sh_script/mzid_conver.sh {mzid_file} {mzid_file.split(".")[0]}'
            subprocess.run(modify_order, shell=True)

            convert_oder = f'/home/qwu/soft/pwiz/idconvert {mzid_file.split(".")[0]}_modify.mzid --pepXML -o {self.Soft_workpath}'
            subprocess.run(convert_oder, shell=True)

    @timing_decorator
    def runPipeline(self,num_task=10):
        # generate fasta file
        mzML_list = [os.path.join(self.Soft_workpath, x) for x in os.listdir(self.Soft_workpath) if x.endswith('.mzML')]
        pool = multiprocessing.Pool(num_task)  # 3个进程执行
        for mzML in mzML_list:
            if not os.path.exists(f'{self.Soft_workpath}/MS_GF_DBsearch.cnlcp'):
                pool.apply_async(func=self.runSoft, args=(mzML,))
                time.sleep(30)
            else:
                pool.apply_async(func=self.runSoft, args=(mzML,))
        pool.close()
        pool.join()

        self.mzidToPepXML()
        self.PeptideProphet()
        self.pepXML_fdr()

def parall_MS_GF(path,sample_list,task_num=15):
    # generate the run order of MS_GF+
    parall_order_list = []
    for sample in sample_list:
        sample_path = os.path.join(path,sample)
        MS_GF_Pipeline(sample_path)
        SoftWork_dir = os.path.join(sample_path,'DB_search_iProphet','MS_GF')
        mzML_list = [os.path.join(SoftWork_dir,x) for x in os.listdir(SoftWork_dir) if x.endswith('.mzML')]
        MS_GF_config = '/data/48/wuqian/TE_DBsearch/config/MSGFPlus_Params_iprophet.txt'
        fasta_file = os.path.join(SoftWork_dir,'MS_GF_DBsearch.fasta')
        for mzML in mzML_list:
            single_order = f'java -jar -Xmx150g /home/qwu/soft/MS-GF+/MSGFPlus.jar -s {mzML} \
            -conf "{MS_GF_config}" -d "{fasta_file}" -o "{mzML.split(".mzML")[0]}.mzid" '
            parall_order_list.append(single_order)

    # parall running the MS-GF+
    pool = multiprocessing.Pool(task_num)  # 3个进程执行
    for order in parall_order_list:
        SoftWork_dir = os.path.dirname(order.split('-o ')[1][1:-2])
        if not os.path.exists(f'{SoftWork_dir}/MS_GF_DBsearch.cnlcp'):
            pool.apply_async(func=subprocess.run, args=(order,), kwds={'shell': True})
            time.sleep(30)
        else:
            pool.apply_async(func=subprocess.run, args=(order,), kwds={'shell': True})
    pool.close()
    pool.join()

    # processing the runnning result
    for sample in sample_list:
        sample_path2 = os.path.join(path,sample)
        print(sample_path2)
        MS_GF_Pipeline(sample_path2).mzidToPepXML()
        MS_GF_Pipeline(sample_path2).PeptideProphet()
        MS_GF_Pipeline(sample_path2).pepXML_fdr()

def Sample_iprophet_blastpFilter(sample_path,type='TE',fdr=0.03):
    # intergted the MS-GF+, Comet and MSfragger results using iprophet
    # sample_path = '/data/48/wuqian/fast/Denovo_test'
    soft_list = ['Comet','MSfragger','MS_GF']
    # type = 'TE'
    pepXML_list = [os.path.join(sample_path,'DB_search_iProphet',soft,f'{soft}_merged_{type}.pep.xml') for soft in soft_list]

    # init
    iprophet_path = os.path.join(sample_path, 'DB_search_iProphet')
    os.chdir(iprophet_path)
    # init / prepare the DB / run in single soft / filter and report
    subprocess.run('/home/qwu/soft/philosopher-master/philosopher workspace --init', shell=True)
    subprocess.run(f'/home/qwu/soft/philosopher-master/philosopher database --annotate {iprophet_path}/MS_GF/MS_GF_DBsearch.fasta',
                   shell=True)

    # iprophet
    subprocess.run(f'/home/qwu/soft/TPP_6.0.0/tpp/bin/RefreshParser {pepXML_list[0]} \
        {iprophet_path}/MS_GF/MS_GF_DBsearch.fasta',shell=True)
    subprocess.run(f'/home/qwu/soft/TPP_6.0.0/tpp/bin/RefreshParser {pepXML_list[1]} \
        {iprophet_path}/MS_GF/MS_GF_DBsearch.fasta',shell=True)
    subprocess.run(f'/home/qwu/soft/philosopher-master/philosopher iprophet --threads 20 \
                        {pepXML_list[0]} {pepXML_list[1]} {pepXML_list[2]}', shell=True)

    # proteinprophet
    subprocess.run(f'/home/qwu/soft/philosopher-master/philosopher proteinprophet --iprophet interact.iproph.pep.xml', shell=True)

    # filter
    subprocess.run(f'/home/qwu/soft/philosopher-master/philosopher filter --pepxml interact.iproph.pep.xml \
     --pep {fdr} --psm {fdr}',shell=True)
    subprocess.run('/home/qwu/soft/philosopher-master/philosopher report',shell=True)

    #blastp
    Blastp_remove(os.path.join(iprophet_path,'peptide.tsv'))

def Sample_iprophet_blastpFilter_Canonical(sample_path,type='Canonical',fdr=0.03):
    # intergted the MS-GF+, Comet and MSfragger results using iprophet
    # sample_path = '/data/48/wuqian/fast/Denovo_test'
    soft_list = ['Comet','MSfragger','MS_GF']
    # type = 'TE'
    pepXML_list = [os.path.join(sample_path,'DB_search_iProphet',soft,f'{soft}_merged_{type}.pep.xml') for soft in soft_list]

    # init
    iprophet_path = os.path.join(sample_path, 'DB_search_iProphet','Canonical')
    subprocess.run(f'mkdir {iprophet_path}',shell=True)
    os.chdir(iprophet_path)
    # init / prepare the DB / run in single soft / filter and report
    subprocess.run('/home/qwu/soft/philosopher-master/philosopher workspace --init', shell=True)
    subprocess.run(f'/home/qwu/soft/philosopher-master/philosopher database --annotate {iprophet_path}/../MS_GF/MS_GF_DBsearch.fasta',
                   shell=True)

    # iprophet
    subprocess.run(f'/home/qwu/soft/TPP_6.0.0/tpp/bin/RefreshParser {pepXML_list[0]} \
        {iprophet_path}/MS_GF/MS_GF_DBsearch.fasta',shell=True)
    subprocess.run(f'/home/qwu/soft/TPP_6.0.0/tpp/bin/RefreshParser {pepXML_list[1]} \
        {iprophet_path}/MS_GF/MS_GF_DBsearch.fasta',shell=True)
    subprocess.run(f'/home/qwu/soft/philosopher-master/philosopher iprophet --threads 20 \
                        {pepXML_list[0]} {pepXML_list[1]} {pepXML_list[2]}', shell=True)

    # filter
    subprocess.run(f'/home/qwu/soft/philosopher-master/philosopher filter --pepxml interact.iproph.pep.xml \
     --pep {fdr}',shell=True)
    subprocess.run('/home/qwu/soft/philosopher-master/philosopher report',shell=True)

def predict_bindingAffinity_iprophet(sample_path,geneType):
    peptide_file = os.path.join(sample_path,'DB_search_iProphet','peptide_blastp.filter.txt')
    if not os.path.exists(peptide_file):
        print(f'{sample_path} is not completed')
        return 1

    df = pd.read_csv(peptide_file,sep='\t')
    if df.__len__() > 0:
        binding_fasta = os.path.join(sample_path,'DB_search_iProphet','binding.fasta')
        if os.path.exists(binding_fasta):
            subprocess.run(f'rm {binding_fasta}', shell=True)
        with open(binding_fasta,'a') as f:
            for idx,row in df.iterrows():
                seq = row['Peptide']
                if  7< seq.__len__() < 15:
                    f.write(f'>{seq}\n{seq}\n')

        # predict using MixMHCpred
        binding_result = binding_fasta.replace("fasta","result")
        predict_sh = f'/home/qwu/soft/MixMHCpred-master/MixMHCpred \
        -i {binding_fasta} \
        -o {binding_result} \
        -a {geneType}'
        print(predict_sh)
        subprocess.call(predict_sh, shell=True)

        # merged the result df
        result_df = pd.read_csv(binding_result,sep='\t',comment='#')
        result_df = result_df[['Peptide','%Rank_bestAllele','BestAllele']]
        rank_dic  = dict(zip(result_df['Peptide'],result_df['%Rank_bestAllele']))
        Alle_dic = dict(zip(result_df['Peptide'],result_df['BestAllele']))
        df['%Rank_bestAllele'] = df['Peptide'].apply(
            lambda x: rank_dic[x] if x in rank_dic.keys() else np.nan
        )
        df['BestAllele'] = df['Peptide'].apply(
            lambda x: Alle_dic[x] if x in Alle_dic.keys() else ''
        )

        df.to_csv(f'{sample_path}/DB_search_iProphet/final_TE_results.csv',sep='\t',index=False)



if __name__ == '__main__':
    pass
