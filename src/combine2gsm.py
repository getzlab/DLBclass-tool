import pandas as pd
import numpy as np
import argparse
from datetime import datetime
import time

start = time.time()
TODAY = datetime.now().strftime("%d%b%Y")

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--id',
                    help='Cohort set name.',
                    required=True, type=str)
parser.add_argument('-v', '--sv_gsm',
                    help='SV GSM.',
                    required=True, type=str)
parser.add_argument('-m', '--mutation_gsm',
                    help='MAF GSM.',
                    required=True, type=str)
parser.add_argument('-c', '--copy_number_gsm',
                    help='SEG GSM.',
                    required=True, type=str)
parser.add_argument('-f', '--feature_order_file',
                    help='GSM Feature order .',
                    required=True, type=str)
#parser.add_argument('-s', '--sample_order_file',
#                    help='GSM Sample order.',
#                    required=True, type=str)
parser.add_argument('-o', '--output_dir',
                    help='Output directory.',
                    required=False, type=str,default='./')
args = parser.parse_args()

scna_df = pd.read_csv(args.copy_number_gsm, sep='\t', index_col=0) #'../../data_tables/gsm/DLBCL.699.scnaGSM.Sep_23_2022.tsv'
col = scna_df.index.values.tolist()
# prefix with "X" if not already prefixed
if not (col[0][0]=="X"):
    xcol = ['X'+x for x in col]
    scna_df.index = xcol

mut_df = pd.read_csv(args.mutation_gsm, sep='\t', index_col=0)     #'../../data_tables/gsm/DLBCL.699.mutationsGSM.Sep_23_2022.tsv'
sv_df = pd.read_csv(args.sv_gsm, sep='\t', index_col=0)            #'../../data_tables/gsm/DLBCL.699.svGSM.Sep_23_2022.tsv'

sample_order = list(sv_df.columns)
sample_order.sort()
scna_df = scna_df.reindex(columns=sample_order)
mut_df = mut_df.reindex(columns=sample_order)
sv_df = sv_df.reindex(columns=sample_order)


full_gsm = pd.concat([mut_df, scna_df, sv_df]) #,ignore_index=True)

# Had to add this after analysis because I was using set() which is not deterministically ordered
#sample_order = pd.read_csv(args.sample_order_file, index_col=0, header=None).index #'../../data_tables/gsm/sample_order.tsv'
feature_order = list(pd.read_csv(args.feature_order_file, index_col=0, header=None).index)   #'../../data_tables/gsm/feature_order.tsv'
#feature_order = ['X'+x for x in feature_order]

gsm = full_gsm.loc[feature_order, sample_order]
gsm = gsm.reset_index().copy()
gsm = gsm.rename(columns={"index":"classifier_name"}).copy()

outfile = args.output_dir + args.id + '.' + TODAY + '.GSM.tsv'
print('output :', outfile)
gsm.to_csv(outfile, sep='\t',  index=False)
