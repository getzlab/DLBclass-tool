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
parser.add_argument('-s', '--sample_set',
                    help='Sample file listing all samples\' in cohort, even w/o SV event.',
                    required=True, type=str)
parser.add_argument('-v', '--sv_file',
                    help='The SV event file (BreakPointer format).',
                    required=True, type=str)
parser.add_argument('-x', '--extra_sv_file_in_gsm_format',
                    help='additional SV events (GSM format).',
                    required=False, type=str,default='')
parser.add_argument('-o', '--output_dir',
                    help='Output directory.',
                    required=False, type=str,default='./')
parser.add_argument('-g','--genome_build',
                    help='Genome build: hg19, hg38.',
                    required=False, type=str, default='hg19')
args = parser.parse_args()

#OUTPUT_FN = '../../data_tables/gsm/DLBCL.699.svGSM.Sep_23_2022.tsv'

#sv_file = '../../data_tables/additional_gsm_inputs/DLBCL_Shipp_Staudt.SVs.14-Dec-2021.txt'
#sample_set = list(pd.read_csv('../../data_tables/train_test_sets/TrainingSet_550Subset_May2021.txt', sep='\t', header=None, index_col=0).index)
#sample_set = set(sample_set + list(pd.read_csv('../../data_tables/train_test_sets/TestingSet_149Subset_May2021.txt', sep='\t', header=None, index_col=0).index))

SV = pd.read_csv(args.sv_file, sep='\t') #, index_col=0)
SV['genes'] = SV['gene1'].fillna('') + '>' + SV['gene2'].fillna('')

SV['evt'] = '---'

SV.loc[SV['genes'].str.contains('BCL2'), 'evt'] = 'SV.BCL2'
SV.loc[SV['genes'].str.contains('BCL6'), 'evt'] = 'SV.BCL6'
SV.loc[SV['genes'].str.contains('MYC'), 'evt'] = 'SV.MYC'

SV = SV.loc[SV['evt'] != '---']

S = pd.read_csv(args.sample_set, sep='\t', index_col=0)
sample_set = list(S.index)

sv_df = pd.DataFrame(0, columns=sample_set, index=['SV.BCL2', 'SV.BCL6', 'SV.MYC'])
alt_counts_df = pd.DataFrame(0, columns=sample_set, index=['SV.BCL2', 'SV.BCL6', 'SV.MYC'])

sv_bcl2_samples = set(SV.loc[(SV['gene1'] == 'BCL2')|(SV['gene2'] == 'BCL2'), 'individual']).intersection(set(sample_set))
sv_bcl6_samples = set(SV.loc[(SV['gene1'] == 'BCL6')|(SV['gene2'] == 'BCL6'), 'individual']).intersection(set(sample_set))
sv_myc_samples = set(SV.loc[(SV['gene1'] == 'MYC')|(SV['gene2'] == 'MYC'), 'individual']).intersection(set(sample_set))

#sv_bcl2_samples = set(SV.loc[SV['genes'].str.contains('BCL2'), 'individual']).intersection(set(sample_set))
#sv_myc_samples = set(SV.loc[SV['genes'].str.contains('MYC'), 'individual']).intersection(sample_set)
#sv_bcl6_samples = set(SV.loc[SV['genes'].str.contains('BCL6'), 'individual']).intersection(sample_set)

sv_df.loc['SV.BCL2', sv_bcl2_samples] = 3
sv_df.loc['SV.BCL6', sv_bcl6_samples] = 3
sv_df.loc['SV.MYC', sv_myc_samples] = 3

for s in sv_bcl2_samples:
    vcf_talt = SV.loc[(SV['individual'] == s) &
                      (SV['genes'].str.contains('BCL2')) , 'VCF_TALT']
    if not vcf_talt.shape[0]:
        print(s)
    alt_counts_df.loc['SV.BCL2', s] = vcf_talt.values[0]

for s in sv_bcl6_samples:
    vcf_talt = SV.loc[(SV['individual'] == s) &
                      (SV['genes'].str.contains('BCL6')), 'VCF_TALT']
    if not vcf_talt.shape[0]:
        print(s)
    alt_counts_df.loc['SV.BCL6', s] = vcf_talt.values[0]

for s in sv_myc_samples:
    vcf_talt = SV.loc[(SV['individual'] == s) &
                      (SV['genes'].str.contains('MYC')) , 'VCF_TALT']
    if not vcf_talt.shape[0]:
        print(s)
    alt_counts_df.loc['SV.MYC', s] = vcf_talt.values[0]

# TCGA SV fix............
if len(args.extra_sv_file_in_gsm_format)>0:
    tcga_svs = pd.read_csv(args.extra_sv_file_in_gsm_format,sep='\t', index_col=0) #../../data_tables/gsm/tcga_svs.tsv'
    sv_df[tcga_svs.columns] = tcga_svs

end = time.time()
print("Execution time :", (end - start), "sec")

gsm = sv_df.reset_index().copy()
gsm = gsm.rename(columns={"index":"classifier_name"}).copy()

outfile = args.output_dir + args.id + '.' + TODAY + '.SV.GSM.tsv'
print('output :')
print( outfile)
gsm.to_csv(outfile, sep='\t', index=False)

alt_counts_df = alt_counts_df.reset_index().copy()
alt_counts_df = alt_counts_df.rename(columns={"index":"classifier_name"}).copy()

outfileA = args.output_dir + args.id + '.' + TODAY + '.SV.alt_counts.tsv'
alt_counts_df.to_csv(outfileA, sep='\t', index=False)
