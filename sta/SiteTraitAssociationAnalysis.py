#!/usr/bin/python
# -*- coding:utf-8 -*-
# Auther: TianYou96
import os
import re
import argparse
import multiprocessing
import sys
import pandas as pd
from tqdm import tqdm
from Bio import AlignIO
from pandas import DataFrame
from collections import Counter
from scipy.stats import chi2_contingency


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
def get_unique_sites_coverage(Trait_Group_seq,Control_Group_seq):
    #count Unique Sites Percentage
    unique_sites_list = []
    for site in Trait_Group_seq:
        if site not in Control_Group_seq:
            unique_sites_list.append(site)
    # Statistics on the percentage of unique amino acids in a specified group to their total number of loci
    if len(unique_sites_list) != 0:
        unique_sites_coverage =(len(unique_sites_list)/len(Trait_Group_seq))
    else:
        unique_sites_coverage = 0
    return unique_sites_coverage
def Generate_Column_Matrix(alignment):
    align_length = len(alignment[0].seq)
    seq_col_list = []
    for site in range(align_length):
        site_list = []
        for align in alignment:
            site_info = align.seq[site]
            site_list.append(site_info)
        seq_col_list.append(site_list)
    return seq_col_list

def significantly_element(Trait_set, Control_set, significance_level=0.05):
    intersection_keys = set(Trait_set.keys()) & set(Control_set.keys())
    significance_element_Trait_set={}
    for element in intersection_keys:
        Count_element_Trait = Trait_set.get(element, 0)
        Count_element_Control = Control_set.get(element, 0)
        total_Trait = int(sum(Trait_set.values()))
        total_Control = int(sum(Control_set.values()))

        element_Trait_per =( Count_element_Trait / total_Trait )/ total_Trait
        element_Control_per =( Count_element_Control / total_Control )/ total_Control

        observed = [Count_element_Control,total_Control]
        expected = [Count_element_Trait,total_Trait]

        chi2, chi2_p_value, dof, expected = chi2_contingency([observed, expected])

        if  chi2_p_value < significance_level and element_Trait_per > element_Control_per:
            #print(f"{element} : {chi2_p_value}")
            significance_element_Trait_set[element]=f"{round((element_Trait_per)*100,2)}%"

    return significance_element_Trait_set

def gte_site_sum_tab(alignment,Trait_Group_id_list,seq_id_list):
    # Generate a matrix by grouping columns
    seq_col_list = Generate_Column_Matrix(alignment)
    site_info_tab = DataFrame(
        columns=['Location', 'Unique Sites Percentage','Potential trait-associated elements', 'Trait Set', 'Control Set', 'Raw list'])
    # Get the row number of the specified class group
    Trait_Group_id_position = []
    for id in Trait_Group_id_list:
        Trait_Group_id_position.append(seq_id_list.index(id))
    # alignment starting from the first column
    col_position = 0
    for col_list in seq_col_list:
        col_position = col_position + 1
        # Obtain a list of site infor for the Trait group
        Trait_Group_id_position_site_list = []
        for i in Trait_Group_id_position:
            Trait_Group_id_position_site_list.append(col_list[i])
        # Obtain a list of site infor for the Control group
        Control_Group_site_position = [i for i in range(len(col_list)) if i not in Trait_Group_id_position]
        Control_Group_site_list = [col_list[i] for i in Control_Group_site_position]
        # other_site_list & Trait_id_position_site_list
        unique_sites_coverage = round((get_unique_sites_coverage(Trait_Group_id_position_site_list, Control_Group_site_list))*100,2)
        Trait_id_site_dict = Counter(Trait_Group_id_position_site_list)
        other_dict = Counter(Control_Group_site_list)
        #most_common_element, most_common_count=Trait_id_site_dict.most_common(1)[0]

        significance_element_Trait_set=significantly_element(Trait_id_site_dict, other_dict)
        if significance_element_Trait_set == {}:
            Potential_unique_elements=''
        else:
            Potential_unique_elements=significance_element_Trait_set
        output_var_info = [col_position, unique_sites_coverage,Potential_unique_elements,
                           dict(Trait_id_site_dict), dict(other_dict), col_list]
        site_info_tab.loc[col_position] = output_var_info
    return site_info_tab

def getdf_site_fst_mpicore(alignment_info_list):
    #alignment_info_list=[alignment,Trait_Group_prefix_id_list,output_dir]

    input_align_path = alignment_info_list[0]
    output_name = input_align_path.split('/')[-1].split('.')[0]
    Trait_Group_prefix_id_list = alignment_info_list[1]
    output_dir = alignment_info_list[2]
    try:
        alignment = AlignIO.read(input_align_path, "fasta")
        #print('Read species position in sequences ...')
        all_seq_id_list=[]
        for align in alignment:
            seq_id = align.id
            all_seq_id_list.append(seq_id)
        #Obtain the ID of the Trait group for this sequence

        pattern = re.compile('|'.join(['^{}.*'.format(prefix) for prefix in Trait_Group_prefix_id_list]))
        Trait_Group_full_id_list = [elem for elem in all_seq_id_list if pattern.match(elem)]

        if len(Trait_Group_full_id_list) > 0:
            site_info_tab=gte_site_sum_tab(alignment,Trait_Group_full_id_list,all_seq_id_list)
            site_info_tab.to_csv('{0}/{1}.site_info.csv'.format(output_dir, output_name))
        else:
            return str("WAR! No specified group : {}".format(output_name))
    except:
        return str("ERR in read MSA: "+str(input_align_path))

def Output_gene_site_tab(input_dir,Number_processes,Trait_Group_prefix_id_list,output_name):
    '''Output of site statistical table for each aligned sequence position'''
    output_dir=mkdir(f'{output_name}')
    alignments_list=list_files_in_folder(input_dir)
    Processing_info_list=[]
    for alignment  in alignments_list:
        alignment_info_list=[alignment,Trait_Group_prefix_id_list,output_dir]
        Processing_info_list.append(alignment_info_list)
    pool = multiprocessing.Pool(processes=int(Number_processes), maxtasksperchild=10000)
    res_info_iter = list(tqdm(pool.imap(func=getdf_site_fst_mpicore, iterable=Processing_info_list), total=len(Processing_info_list),
                                  desc='  Output Site Table', ncols=80))
    pool.close()

    for return_info in res_info_iter:
        if isinstance(return_info,str):
            output_err = open(f"{output_name}.log", 'a+')
            index_number = res_info_iter.index(return_info)
            print(str(str(index_number)+':'+return_info),file=output_err)
            output_err.close()
    return output_dir

def Calculate_Homogenity(Group_dict):
    Number_total_Group = sum(Group_dict.values())
    most_common_element = max(Group_dict, key=Group_dict.get)
    Homogeneity = int(Group_dict[most_common_element]) / Number_total_Group
    return Homogeneity

def Site_evolution_tab(input_site_tab_dir,Trait_Group_prefix_id_list,mode):
    '''Read all site tables and output filtered summary tables
    '''
    site_tab_list = list_files_in_folder(input_site_tab_dir)
    sum_div_site_tab = DataFrame(columns=['Gene name','Location', 'USP','OEA','DES',
                                         'Trait-Homo.','Control-Homo.','Score',
                                          'Dominant element','Potential trait-associated elements',
                                         'Statistics of Trait','Statistics of Control'])

    codon_map_table = DataFrame(columns=['Gene name','Location','Score', 'Codon Start Loc.','Codon End Loc.'])
    Expectations_Trait_Group=len(Trait_Group_prefix_id_list)
    index_number=0
    for site_tab_path in tqdm(site_tab_list, desc="  Summary Site Table", unit="item", ncols=80):
        input_site_tab_path=site_tab_path
        df_site_info = pd.read_csv(input_site_tab_path,index_col=0)
        #
        info_site_df = df_site_info[
            ((df_site_info['Unique Sites Percentage']) == 100) | (df_site_info['Potential trait-associated elements'].notnull())]
        gene_name=str(site_tab_path).split('/')[-1].split('.')[0]
        for index, row in info_site_df.iterrows():
            try:
                Trait_Group_dict = dict(row['Trait Set'])
                Control_Group_dict = dict(row['Control Set'])

            except:
                Trait_Group_dict = eval(row['Trait Set'])
                Control_Group_dict = eval(row['Control Set'])

            Observations_Trait_Group = sum(Trait_Group_dict.values())
            Sample_Coverage = Observations_Trait_Group/Expectations_Trait_Group
            Trait_Homogeneity = round((Calculate_Homogenity(Trait_Group_dict))*100,2)
            Control_Homogeneity = round((Calculate_Homogenity(Control_Group_dict))*100,2)
            USP=row['Unique Sites Percentage']

            #Dominant_Element
            most_common_element = max(Trait_Group_dict, key=Trait_Group_dict.get)
            if int(Trait_Group_dict[most_common_element])/Observations_Trait_Group > 0.5:
                Dominant_Element = most_common_element
            else:
                Dominant_Element =''

            # dominant element state (DES)
            if USP == 100:
                Predict_element_consistency = 100
            elif 100 > USP:
                Predict_element_consistency = 0
            #mk Potential_unique_elements
            if row['Potential trait-associated elements'] == '' or pd.isnull(row['Potential trait-associated elements']):
                Potential_unique_elements=''
            else:
                try:
                    Potential_unique_elements_dict = dict(row['Potential trait-associated elements'])
                except:
                    Potential_unique_elements_dict = eval(str(row['Potential trait-associated elements']))
                Consistency_predicted_elements_list = list(Potential_unique_elements_dict.keys())

                if Dominant_Element in Consistency_predicted_elements_list:
                    Predict_element_consistency = 50
                else:
                    Predict_element_consistency = 0

                if len(Consistency_predicted_elements_list) >1 :
                    Potential_unique_elements = '|'.join(list(Potential_unique_elements_dict.keys()))
                else:
                    Potential_unique_elements = ''.join(list(Potential_unique_elements_dict.keys()))

            Score=round(( (Sample_Coverage*0.2) + (Trait_Homogeneity*0.2) + (Control_Homogeneity*0.2) +
                       (row['Unique Sites Percentage'])*0.2+ (Predict_element_consistency)*0.2),2)
            sum_div_site_tab.loc[index_number] = [gene_name,row['Location'],USP,Sample_Coverage,Predict_element_consistency,
                                                  Trait_Homogeneity,Control_Homogeneity,
                                                  Score,Dominant_Element,Potential_unique_elements,row['Trait Set'],row['Control Set']]
            index_number+=1

            algin_codon_End = int(row['Location']) * 3
            algin_codon_Start = algin_codon_End - 2
            codon_map_table.loc[index_number] = [gene_name, row['Location'],Score, algin_codon_Start, algin_codon_End]

    sum_div_site_tab=sum_div_site_tab.sort_values(by=['Score','OEA','Trait-Homo.'],ascending=False)

    if mode == 'AA':
        return sum_div_site_tab
    elif  mode == 'Codon':
        return sum_div_site_tab, codon_map_table



def Site_Trait_Association_module(input_align_dir,Trait_Group_id_str,Number_processes,output_name,mode):
    input_align_path = Determine_input_folder_exists(input_align_dir)
    Trait_Group_id_list = Trait_Group_id_str.split(',')
    Trait_Group_id_list = [id.strip() for id in Trait_Group_id_list]
    print(f"\tTrait set prefix id list :{Trait_Group_id_list}")
    gene_site_tab_dir=Output_gene_site_tab(input_align_path,int(Number_processes),Trait_Group_id_list,output_name)
    Disjoint_output_name = 'Focal_Sites'
    if mode == 'AA':
        Disjoint_Site_table=Site_evolution_tab(gene_site_tab_dir,Trait_Group_id_list,mode)
        Disjoint_Site_table.to_csv(f"{Disjoint_output_name}.csv",index=False)
        return Disjoint_Site_table
    elif  mode == 'Codon':
        Disjoint_Site_table, codon_map_table = Site_evolution_tab(gene_site_tab_dir, Trait_Group_id_list, mode)
        Disjoint_Site_table.to_csv(f"{Disjoint_output_name}_summary.csv",index=False)
        codon_map_table.to_csv(f"{Disjoint_output_name}_codon_map.csv",index=False)
        return Disjoint_Site_table, codon_map_table



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Site-trait association analysis module in SiteKIT')
    parser.add_argument('--aa_dir',
                        type=str,required=True,
                        help='Input folder include MSA of amino acid ')
    parser.add_argument( '--group_id',
                        type=str,required=True,
                        help="Specified IDs or prefixes for trait group. \n"
                             "Include Trait Group prefix list in single quotes('') separated by commas(,). \n"
                             "Do not include spaces. \n"
                             "For example: 'sample1,sample2,...' ")
    parser.add_argument('-o','--output_prefix',
                        type=str,default='2-AA_site_USP',
                        help='Output folder name. [Default= 2-AA_site_USP ]')
    parser.add_argument( '--multiprocessing',
                        type=str,default=6,
                        help="The number of multiple processes \n")
    parser.add_argument('-m','--mode',
                        type=str,default='AA',
                        help="AA: only folder MSA of amino acid. [Default]; "
                        "Codon: additional folder MSA of codon. "
                        )
    args = parser.parse_args()

    input_align_dir=os.path.abspath(args.aa_dir)
    Trait_Group_id_str=args.group_id
    Number_processes=args.multiprocessing
    output_name=args.output_prefix
    mode=args.mode

    Site_Trait_Association_module(input_align_dir, Trait_Group_id_str, Number_processes, output_name,mode)




