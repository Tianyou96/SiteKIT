#!/usr/bin/python
# -*- coding:utf-8 -*-
# Auther: TianYou96
import os
import re
import argparse
import sys
from tqdm import tqdm
from collections import Counter


def mkdir(dir_name):
    current_dir=os.getcwd()
    New_Dir_Path=os.path.join(current_dir, dir_name)
    if not os.path.exists(New_Dir_Path):
        os.makedirs(New_Dir_Path)
    else:
        print('Directory already exists : {}'.format(New_Dir_Path))
        sys.exit(0)
    return New_Dir_Path
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
def get_file_lines(infile_name):
    infile = open(infile_name, "r")
    inlines_list = infile.readlines()
    infile.close()
    cleaned_lines = [line.strip() for line in inlines_list]
    return cleaned_lines

def get_fasta_dict_and_check(input_fasta):
    lines = get_file_lines(input_fasta)
    current_id = None
    sequence_dict = {}
    all_id_list= []
    Important_error = None

    output_err = open("{0}.log".format('err'), 'a+')

    for line in lines:
        if line.startswith('>'):
            #Create a ID
            current_id = line[1:]
            sequence_dict[current_id] = ""
            if current_id in all_id_list:
                print(f'ERR Duplicate ID in {input_fasta} : {current_id}', file=output_err)
                Important_error = True
            else:
                all_id_list.append(current_id)
        else:
            if ">" in str(line):
                #If the sequence and ID are mixed
                seq_line=line.split('>')[0]
                sequence_dict[current_id] += seq_line
                id_line=line.split('>')[1]
                # Create a ID
                current_id=id_line
                sequence_dict[current_id] = ""
                if current_id in all_id_list:
                    print(f'ERR Duplicate ID in {input_fasta} : {current_id}', file=output_err)
                    Important_error = True
                else:
                    all_id_list.append(current_id)
            else:
                #Append sequence
                sequence_dict[current_id] += line
    #Check for empty values or inconsistent lengths
    sequence_length_list=[]
    for seq_id, sequence in sequence_dict.items():
        if sequence == "":
            print(f'ERR No sequence in {input_fasta} : {seq_id}', file=output_err)
            Important_error = True
        sequence_length = len(sequence)
        sequence_length_list.append(sequence_length)

    count_sequence_lengths = Counter(sequence_length_list)
    Most_sequence_lengths= count_sequence_lengths.most_common(1)[0][0]

    for seq_id, sequence in sequence_dict.items():
        sequence_length = len(sequence)

        if sequence_length == Most_sequence_lengths:

            continue
        else:
            #print(f'{sequence_length}/{Most_sequence_lengths}')
            print(f'ERR Inconsistent sequence length in {input_fasta} : {seq_id}', file=output_err)
            Important_error = True

    output_err.close()
    if Important_error == True:
        return 'Important_error'
    else:
        return sequence_dict


def keep_id(input_fasta, output_name):
    sequence_dict = get_fasta_dict_and_check(input_fasta)
    if sequence_dict == 'Important_error':
        return False
    else:
        #print(f"\t {input_fasta}")
        output_fas = open("{0}.fasta".format(output_name), 'a+')
        for seq_id, sequence in sequence_dict.items():
            old_seq_id = str(seq_id)
            idline = '>' + old_seq_id
            print(idline, file=output_fas)
            print(str(sequence), file=output_fas)
        output_fas.close()
def relapce_id(input_fasta, output_name,Separator):
    sequence_dict = get_fasta_dict_and_check(input_fasta)
    if sequence_dict == 'Important_error':
        return False
    else:
    #print(f"\t {input_fasta}")
        output_fas = open("{0}.fasta".format(output_name), 'a+')
        for seq_id, sequence in sequence_dict.items():
            if Separator == 'no':
                old_seq_id=str(seq_id)
            else:
                old_seq_id = str(seq_id).split(Separator)[0].split('\t')[0]
            new_seq_id = re.sub(r'[^a-zA-Z0-9]', '_', old_seq_id)
            idline = '>' + new_seq_id
            print(idline, file=output_fas)
            print(str(sequence), file=output_fas)
        output_fas.close()

def Fasta_Tool_module(DIR_fasta,output_name,mode,Separator):
    DIR_fasta_path=Determine_input_folder_exists(DIR_fasta)
    file_list=list_files_in_folder(DIR_fasta_path)
    output_path=mkdir(output_name)
    if mode == 'mode1':
        for input_fasta in tqdm(file_list, desc="  Check fasta", unit="item", ncols=80):
            file_name=str(input_fasta).split('/')[-1].split('.')[0]
            output_path_name=os.path.join(output_path, file_name)
            keep_id(input_fasta, output_path_name)
    elif mode == 'mode2':
        for input_fasta in tqdm(file_list, desc="  Check fasta", unit="item", ncols=80):
            file_name=str(input_fasta).split('/')[-1].split('.')[0]
            output_path_name=os.path.join(output_path, file_name)
            relapce_id(input_fasta, output_path_name,Separator)
    return output_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Fasta tool in SiteKIT')
    parser.add_argument('--aa_dir',
                        type=str,required=True,
                        help='Input folder include MSA of amino acid ')
    parser.add_argument('-o','--output_prefix',
                        type=str,default='1-Standard_MSA',
                        help='Output folder name. [Default= 1-Standard_MSA ]')
    parser.add_argument('-m','--mode',
                        type=str,default='mode1',
                        help="mode1: Preserves the sequence ID exactly as it is. [Default]; "
                        "mode2: Replaces any characters in the sequence ID that are not alphanumeric with '_' (underscore). "
                        "Meanwhile, only retain the content before the first space."
                        )
    parser.add_argument('--separator',
                        type=str,default='no',
                        help="It is not enabled by default, and it will be enabled after specifying the separator in mode 2. "
                        )
    args = parser.parse_args()

    DIR_fasta=args.aa_dir
    output_name = args.output_prefix
    mode = args.mode
    Separator = args.separator

    Fasta_Tool_module(DIR_fasta,output_name,mode,Separator)


