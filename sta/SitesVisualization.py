#!/usr/bin/python
# -*- coding:utf-8 -*-
# Auther: TianYou96
import argparse
import os
import sys

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

def Determine_input_folder_exists(input_folder):
    if not os.path.exists(input_folder):
        print('Input folder not exists : {}'.format(input_folder))
        sys.exit(0)
    else:
        input_abs_path=os.path.abspath(input_folder)
        return input_abs_path

def list_files_in_folder(folder_path):
    file_list = []
    absolute_path=os.path.abspath(folder_path)
    for file in os.listdir(folder_path):
        file_list.append(os.path.join(absolute_path, file))
    return file_list

def mkdir(dir_name):
    current_dir=os.getcwd()
    New_Dir_Path=os.path.join(current_dir, dir_name)
    if not os.path.exists(New_Dir_Path):
        os.makedirs(New_Dir_Path)
    else:
        print('directory already exists : {}'.format(New_Dir_Path))
        sys.exit(0)
    return New_Dir_Path

def output_itol(input_fasta,site,output_dir):

    gene=str(input_fasta).split('/')[-1].split('.')[0]
    gene_site_label = f"{gene}_{site}"

    sequence_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
    seq_id_list = list(sequence_dict.keys())
    align_length=len(sequence_dict[seq_id_list[0]].seq)

    if site < 5:
        start_position = 0
        end_position = 10
    elif site > int(align_length-5):
        start_position = align_length-10
        end_position = align_length
    else:
        start_position = site-5
        end_position = site+5

    dataset_alignment='''DATASET_ALIGNMENT
    SEPARATOR COMMA
    CUSTOM_COLOR_SCHEME,MY_SCHEME_1,A=#d2d0c9,M=#d2d0c9,I=#d2d0c9,L=#d2d0c9,V=#d2d0c9,P=#746f69,G=#746f69,C=#746f69,F=#d0ad16,Y=#d0ad16,W=#d0ad16,S=#34acfb,T=#34acfb,N=#34acfb,Q=#34acfb,R=#34fb54,K=#34fb54,H=#34fb54,D=#fb4034,E=#fb4034 \t
    CUSTOM_COLOR_SCHEME,MY_SCHEME_2,A=#30a040,R=#2015a5,N=#10ffa0,D=#6048c0,C=#608a80,Q=#601f00,E=#5048c0,G=#702048,H=#a5a4a4,I=#1a30f0,L=#11a0f0,K=#003505,M=#00a0f0,F=#10a0f0,P=#0ff300,S=#93f300,T=#33ff00,W=#0a30f0,Y=#25a4a4,V=#90a3f0 \t
    SIZE_FACTOR,1
    COLOR_SCHEME,clustal
    DISPLAY_CONSENSUS,0
    DISPLAY_CONSERVATION,0
    '''
    DATASET_LABEL=f'DATASET_LABEL,{gene_site_label} \nCOLOR,#ff0000'
    site_START=f"START_POSITION,{start_position}"
    site_END=f"END_POSITION,{end_position}"

    if align_length < 4000:
        gene_site_iTOL = open(f"{output_dir}/{gene_site_label}_iTOL.txt", 'w+')
        print(f"{dataset_alignment}\n{DATASET_LABEL}\n{site_START}\n{site_END}\nDATA\n", file=gene_site_iTOL)
        seq_id_list = list(sequence_dict.keys())
        for seq_id in seq_id_list:
            sequence = str(sequence_dict[seq_id].seq)
            seq_id = str(seq_id)
            idline = '>' + seq_id
            print(idline, file=gene_site_iTOL)
            print(str(sequence), file=gene_site_iTOL)
        gene_site_iTOL.close()



def Tree_Site_Plot_module(input_dir,gene_site_df,output_name):
    input_fasta_path = Determine_input_folder_exists(input_dir)
    fasta_list = list_files_in_folder(input_fasta_path)
    output_dir = mkdir(f'{output_name}')
    for index, row in tqdm(gene_site_df.iterrows(), total=len(gene_site_df), desc="  Output iTOL txt", ncols=80):
        gene=row['Gene name']
        site=int(row['Location'])

        input_fasta=next((item for item in fasta_list if gene in item), None)
        if input_fasta == None:

            continue
        else:

            output_itol(input_fasta, site, output_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Sites visualization module in SiteKIT')
    parser.add_argument('--aa_dir',
                        type=str,required=True,
                        help='Input folder include MSA of amino acid ')
    parser.add_argument('--gene_site_table',
                        type=str,default='',required=False,
                        help='gene site table of [Disjoint_Site_summary.csv], it needs to be in CSV format. ')
    parser.add_argument('-o','--output_prefix',
                        type=str,default='3-Disjoint_site_iTOL',
                        help='Output folder name. [Default= 3-Disjoint_site_iTOL ]')
    args = parser.parse_args()
    input_aa_dir = os.path.abspath(args.aa_dir)
    input_csv=args.gene_site_table
    output_name=args.output_prefix

    gene_site_df = pd.read_csv(input_csv, index_col=0)
    Tree_Site_Plot_module(input_aa_dir, gene_site_df,output_name)