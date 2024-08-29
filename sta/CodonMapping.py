#!/usr/bin/python
# -*- coding:utf-8 -*-
# Auther: TianYou96
import re
import os
import argparse
import sys

import pandas as pd
from pandas import DataFrame
from Bio import SeqIO, AlignIO
from tqdm import tqdm


def mkdir(dir_name):
    current_dir=os.getcwd()
    New_Dir_Path=os.path.join(current_dir, dir_name)
    if not os.path.exists(New_Dir_Path):
        os.makedirs(New_Dir_Path)
    else:
        print('directory already exists : {}'.format(New_Dir_Path))
        sys.exit(0)
    return New_Dir_Path
def Determine_input_folder_exists(input_folder):
    if not os.path.exists(input_folder):
        print('Input folder not exists : {}'.format(input_folder))
        sys.exit(0)
    else:
        input_abs_path=os.path.abspath(input_folder)
        return input_abs_path

def read_align_to_dict(alignment_folder,id_list):

    matched_alignments_dict={}
    files = os.listdir(alignment_folder)

    for gene_id in id_list:
        pattern = re.compile(r'^{}.*'.format(gene_id))
        for file in files:
            if pattern.match(file):
                input_align=os.path.join(alignment_folder, file)
                alignment_list = list(AlignIO.read(input_align, "fasta"))
                matched_alignments_dict[gene_id]=alignment_list
    return matched_alignments_dict


def Map_Codon_module(site_df,alignment_folder,output_dir):
    output_dir = mkdir(output_dir)
    alignment_folder_path = Determine_input_folder_exists(alignment_folder)

    id_list=list(set(list(site_df['Gene name'])))
    matched_alignments_dict=read_align_to_dict(alignment_folder_path,id_list)


    for gene_id in tqdm(id_list, desc="    Output Codon Tables", unit="item", ncols=80):
        gene_site_codon_tab=site_df[site_df['Gene name'] == gene_id]
        alignment_list=matched_alignments_dict[gene_id]
        seq_id_list=[align.id for align in alignment_list]
        columns_list=['Gene name','Score', 'Codon Position'] + seq_id_list
        gene_site_codon_output_table = DataFrame(columns=columns_list)
        for index, row in gene_site_codon_tab.iterrows():
            Codon_Start = int(row['Codon Start Loc.'])-1
            Codon_End = row['Codon End Loc.']
            align_info_dict={}
            for align_info in alignment_list:
                align_id=align_info.id
                codon = align_info.seq[Codon_Start:Codon_End]
                align_info_dict[align_id]=codon
            Codon_Position = f"{(row['Codon Start Loc.'])}-{Codon_End}"
            codon_list = [align_info_dict[seq_id] for seq_id in seq_id_list]
            gene_site_codon_output_table.loc[index] = [gene_id,row['Score'], Codon_Position] + codon_list
        gene_site_codon_output_table.to_csv(f"{output_dir}/{gene_id}_codon.csv")


if __name__ == "__main__":
    # 设置输入使用argparse
    parser = argparse.ArgumentParser(description='Codon mapping module in SiteKIT')
    parser.add_argument('--codon_dir',
                        type=str,required=False,
                        help='Input folder include MSA of codon')

    parser.add_argument('--codon_site_table',
                        type=str,required=False,
                        help='gene site table of [Disjoint_Site_codon_map.csv], it needs to be in CSV format.')
    parser.add_argument('-o','--output_prefix',
                        type=str,default='4-Disjoint_site_Codon',
                        help='Output folder name. [Default= 4-Disjoint_site_Codon ]')

    args = parser.parse_args()

    alignment_folder = os.path.abspath(args.codon_dir)
    output_dir=args.output_prefix
    input_csv=args.Codon_site_table
    site_df = pd.read_csv(input_csv)
    Map_Codon_module(site_df, alignment_folder, output_dir)






