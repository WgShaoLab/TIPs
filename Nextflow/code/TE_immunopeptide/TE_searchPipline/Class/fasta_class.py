import re,os
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
import time

class Fasta:
    def __init__(self, raw_fasta_path,fasta_source='',uniform_PE_tag='',other_tag='', whether_subfasta=False):
        self.raw_fasta_path = raw_fasta_path
        self.fasta_source = fasta_source
        self.whether_subfasta = whether_subfasta
        ###
        self.fastaFile_str = open(self.raw_fasta_path, 'r').read()
        ###
        self.record_num = self.fastaFile_str.split('>')[1:].__len__()
        # 是否是多行的fasta，如果是多行的fasta，自动转为单行的fasta
        if self.fastaFile_str.split('\n').__len__() > self.record_num*2 + 4:
            self.tranMultiFastaTosingle()
        # 将PE标签统一为tag值
        if uniform_PE_tag:
            self.Uniform_PE_tag(uniform_PE_tag)
        # 添加其他的tag
        if other_tag:
            self.add_other_tag(other_tag)
        # 统计
        PE_list = re.findall(r'PE=(.)', self.fastaFile_str)
        self.PE_list = list(set(PE_list))

    def tranMultiFastaTosingle(self):
        fasta_list = self.fastaFile_str.split('>')[1:]
        tmp_fastaFile_str = ''
        for record in fasta_list:
            line_list = record.split('\n')
            tmp_fastaFile_str += '>' + line_list[0] +'\n' + ''.join(line_list[1:]) + '\n'
        self.fastaFile_str = tmp_fastaFile_str
        print('transform the multi fasta file to single fasta')

    def Uniform_PE_tag(self,tag):
        fasta_list = self.fastaFile_str.split('>')[1:]
        tmp_fastaFile_str = ''
        for record in fasta_list:
            if 'PE=' in record:
                record = re.sub(' PE=(\d) ',f' PE={tag} ',record)
                tmp_fastaFile_str += '>' + record
            else:
                header = record.split('\n')[0] + f' PE={tag}'
                seq = record.split('\n')[1]
                tmp_fastaFile_str += '>' + header + '\n' + seq + '\n'
        self.fastaFile_str = tmp_fastaFile_str
        print('Uniform the fasta PE tage')

    def add_other_tag(self,other_tag):
        fasta_list = self.fastaFile_str.split('>')[1:]
        tmp_fastaFile_str = ''
        for record in fasta_list:
            header = record.split('\n')[0] + f' {other_tag}'
            seq = record.split('\n')[1]
            tmp_fastaFile_str += '>' + header + '\n' + seq + '\n'
        self.fastaFile_str = tmp_fastaFile_str
        print(f'add {other_tag} tag in the fasta')

    def split_fasta(self, output_dir, num_files=15):
        '''
        将单个的fasta拆分，以利于fasta的并行blastp
        '''
        # 读取FASTA文件中的所有记录
        records = self.fastaFile_str.split('>')[1:]

        # 计算每个子文件应该有的记录数量
        records_per_file = len(records) // num_files
        remainder = len(records) % num_files  # 处理记录数不能均分的情况

        start = 0
        for i in range(num_files):
            # 计算每个文件应该包含的记录数量
            end = start + records_per_file + (1 if i < remainder else 0)
            sub_records = records[start:end]
            start = end

            # 输出子文件
            output_fasta = os.path.join(output_dir, f"sub_fasta_{i + 1}.fasta")
            with open(output_fasta, "w") as output_handle:
                sub_fasta = ''
                for i in sub_records:
                    sub_fasta += '>' + i
                output_handle.write(sub_fasta)
        print(f'split fasta done')

    def Denovo_rawPeptide_blastp_IL(self,blastp_path='/home/qwu/anaconda3/envs/blast/bin/blastp',
        blastp_db_path='/home/qwu/data/6_frame/hg38_rmsk_6_frame_singleLine_done.fasta',num_task=15):
        # /data/48/wuqian/genome/human/hg38/repeat_anno/6_frame/hg38_rmsk_6_frame_singleLine_done.fasta
        '''
        :param fasta_file: 其中head和seq均为peptide序列
        :return:
        '''
        def I_L_identity(q_seq, t_seq):
            if q_seq == t_seq:
                return True
            diff_list = []
            for i in range(len(q_seq)):
                if q_seq[i] != t_seq[i] and set([q_seq[i], t_seq[i]]) == set(['I', 'L']):
                    diff_list.append(t_seq[i])
            diff_list = list(set(diff_list))

            if not diff_list:
                return False
            else:
                return ''.join(sorted(diff_list)) in 'IL'

        child_fasta_dir = self.raw_fasta_path.split('.')[0]+"_Child_fasta"
        subprocess.run(f'mkdir {child_fasta_dir}', shell=True)
        self.split_fasta(child_fasta_dir, num_files=num_task)

        # 对于每一个subfasta进行blastp

        for i in [x for x in os.listdir(child_fasta_dir) if 'sub_fasta_' in x]:
            blastp_order = f"nohup {blastp_path} -task blastp-short \
            -query {os.path.join(child_fasta_dir, i).split('.')[0]}.fasta \
            -db {blastp_db_path} \
            -out {os.path.join(child_fasta_dir, i).split('.')[0]}.txt \
            -outfmt ' 6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'\
            -evalue 20000 -num_threads 4 &"
            subprocess.run(blastp_order, shell=True)

        # 检查blastp是否运行完成
        try:
            while int(subprocess.check_output("ps aux | grep -v grep | grep -c 'blastp'", shell=True)) > 0:
                time.sleep(5)
        except:
            print('Blastp Done')

        # 合并blastp的结果文件并使用awk进行筛选
        merged_order = f'cat {child_fasta_dir}/*.txt > {child_fasta_dir}/merged_blastp.txt'
        subprocess.run(merged_order, shell=True)
        awk_filter_order = f"awk '$7 == 1 && $6 == 0 && $3 >= 75' {child_fasta_dir}/merged_blastp.txt > \
                            {child_fasta_dir}/merged_blastp_filter.txt"
        subprocess.run(awk_filter_order, shell=True)

        # 处理blastp的dataframe
        column_names = ['qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart',
                        'send', 'evalue', 'bitscore', 'qseq', 'sseq']
        df = pd.read_csv(f"{child_fasta_dir}/merged_blastp_filter.txt", names=column_names, sep='\t')
        df['q_length'] = df['qaccver'].str.len()
        df['qseq_I2L'] = df['qseq'].apply(lambda x: x.replace('I','L'))
        df['sseq_I2L'] = df['sseq'].apply(lambda x: x.replace('I', 'L'))
        df_filter = df[df['length'] == df['q_length']]
        df_filter = df_filter[df_filter['qseq'].str.len() == df_filter['sseq'].str.len()]
        df_filter = df_filter[df_filter['gapopen'] == 0]
        df_filter['I/L_identity'] = df_filter.apply(lambda x:x['qseq_I2L']== x['sseq_I2L'], axis=1)

        df_filter.to_csv(f"{child_fasta_dir}/merged_blastp_TE.csv", index=False, sep='\t')
        return 0

    def generateOther_fasta(self,output_path):
        # 在原fasta的基础上生成decoy序列
        with open(output_path, 'a') as out_file:
            record_list = []
            for record in SeqIO.parse(self.raw_fasta_path,'fasta'):
                # 生成反向序列
                decoy_seq = str(record.seq)[::-1]
                record_list.append(f'>rev_{record.id}\n{decoy_seq}\n')
            out_file.writelines(''.join(record_list))

    def output(self,new_fasta_path):
        with open(new_fasta_path, "w") as output_handle:
            output_handle.write(self.fastaFile_str)




if __name__ =='__main__':
    pass


