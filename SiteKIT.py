#!/usr/bin/python
# -*- coding:utf-8 -*-
# Auther: TianYou96
import os
import argparse
import sys

from sta import FastaTool
from sta import CodonMapping
from sta import SitesVisualization
from sta import SiteTraitAssociationAnalysis


def mian():

    print(f"\t---------------------------\n"
          f"\tWelcome to SiteKIT (￣∇￣)! \n"
          f"\t---------------------------")

    if args.aa_dir is None or args.group_id is None:
        parser.print_help()
        sys.exit(0)

    input_aa_dir = os.path.abspath(args.aa_dir)
    output_name = args.output_prefix
    mode = args.mode
    Separator = args.separator
    Trait_Group_id_str = args.group_id
    Number_processes = args.multiprocessing
    align_mode='AA'
    if args.codon_dir is not None:
        input_codon_dir = os.path.abspath(args.codon_dir)
        align_mode = 'Codon'
    output_path=FastaTool.mkdir(output_name)
    os.chdir(output_path)


    if align_mode == 'AA':
        print(f"Step 1 : Standardize and Validate Input Sequence")
        Standard_input_AA = FastaTool.Fasta_Tool_module(input_aa_dir, '1-Standard_MSA_AA', mode, Separator)
        print(f"Step 2 : Generate Site Information and Disjoint_Site Summary Tables")
        Disjoint_Site_table = SiteTraitAssociationAnalysis.Site_Trait_Association_module(Standard_input_AA,Trait_Group_id_str,
                                                                                   Number_processes,'2-AA_site_USP',align_mode)

        print(f"Step 3 : Create iTOL Visualization File")
        SitesVisualization.Tree_Site_Plot_module(Standard_input_AA, Disjoint_Site_table,'3-Focal_sites_iTOL')
    elif align_mode == 'Codon':
        print(f"Step 1 : Standardize and Validate Input Sequence")
        Standard_input_AA = FastaTool.Fasta_Tool_module(input_aa_dir, '1-Standard_MSA_AA', mode, Separator)
        Standard_input_codon=FastaTool.Fasta_Tool_module(input_codon_dir, '1-Standard_MSA_Codon', mode, Separator)
        print(f"Step 2 : Generate Site Information and Disjoint Site Summary Tables")
        Disjoint_Site_table, codon_map_table = SiteTraitAssociationAnalysis.Site_Trait_Association_module(Standard_input_AA,Trait_Group_id_str,
                                                                                   Number_processes,'2-AA_site_USP',align_mode)

        print(f"Step 3 : Create iTOL Visualization File")
        SitesVisualization.Tree_Site_Plot_module(Standard_input_AA, Disjoint_Site_table,'3-Focal_sites_iTOL')
        print(f"Step 4 : Output Codon Mapping Table")
        CodonMapping.Map_Codon_module(codon_map_table, Standard_input_codon, '4-Focal_sites_Codon')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='SiteKIT for site-trait association analysis')
    parser.add_argument('--aa_dir',
                        type=str,required=True,
                        help='Input folder include MSA of amino acid ')
    parser.add_argument( '--group_id',
                        type=str,required=True,
                        help="Specified IDs or prefixes for trait group. \n"
                             "Include Trait Group prefix list in single quotes('') separated by commas(,). \n"
                             "Do not include spaces. \n"
                             "For example: 'sample1,sample2,...' ")

    parser.add_argument('--codon_dir',
                        type=str,required=False,
                        help='Input folder include MSA of codon')

    parser.add_argument('-o','--output_prefix',
                        type=str,default='Site_Trait_Association_Analysis',
                        help='Output folder name. [Default= Site_Trait_Association_Analysis]')
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
    parser.add_argument( '--multiprocessing',
                        type=str,default=10,
                        help="The number of multiple processes \n")
    args = parser.parse_args()

    mian()









