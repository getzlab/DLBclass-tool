import argparse
import pandas as pd
import numpy as np
from datetime import datetime
TODAY = datetime.now().strftime("%d%b%Y")
import time

start = time.time()
# --maf ../../data_tables/maf_files/DLBCL_combined_699.hg38B.noPDE4DIP.noHISTartifacts.maf --mutsig_sig_genes ../../data_tables/mutsig2cv_gistic_qvalues/DLBCL_550_training_noPDE4DIP_noHISTartifacts.sig_genes.txt --mutsig_q_thresh 0.10 --output_fn ../../data_tables/gsm/DLBCL.699.mutationsGSM.Sep_23_2022.tsv --additional_sig_genes ../../data_tables/additional_gsm_inputs/NatMed_104_sig_genes.tableS3.tsv --include_myd88_L265P --include_ccf CCF_hat --blacklist ../../data_tables/additional_gsm_inputs/dlbclass_blacklist.tsv --ploidy ../../data_tables/purities_ploidies/PloidyDataFrame.txt --purity ../../data_tables/purities_ploidies/ALLPurities_fixednames.tsv --coo ../../data_tables/phenotypes/COO_and_genetic_lables.txt


parser = argparse.ArgumentParser()
parser.add_argument('--id',
                    help='Cohort set name.',
                    required=True, type=str)
parser.add_argument('-s', '--sample_set',
                    help='Sample file listing all samples\' in cohort, even w/o SV event.',
                    required=True, type=str)
parser.add_argument('-m','--maf',
                    help='The maf file (MAF format) to call your samples\' mutation events.',
                    required=True, type=str)
parser.add_argument('-g','--gsm_genes',
                    help='GSM gene SNVs and indels  to be included .',
                    required=False, type=str, default='')
parser.add_argument('-o','--output_dir',
                    help='Output directory.',
                    required=False, type=str, default='./')

args = parser.parse_args()

cols = ['Hugo_Symbol', 'Tumor_Sample_Barcode', 'Variant_Classification', 'Protein_Change']

maf = pd.read_csv(args.maf, sep='\t', dtype=str, low_memory=False, usecols=cols)
S = pd.read_csv(args.sample_set, sep='\t', index_col=0)
sample_set = list(S.index)

# determine which samples to call events for
samples = set(maf['Tumor_Sample_Barcode'])

# original code removed Splice_Site outside of exons, we want to keep them in the GSM  
#maf = maf.loc[~maf['Protein_Change'].isna()]
Functional_Variant_Classifications=["Splice_Site","De_novo_Start_InFrame","De_novo_Start_OutOfFrame","Start_Codon_Del"]
maf = maf.loc[(maf['Protein_Change'].notna()) | (maf['Variant_Classification'].isin(Functional_Variant_Classifications))]

GENES = ['ACTB', 'ATP2A2', 'B2M', 'BCL10', 'BCL11A', 'BCL2', 'BCL6', 'BCL7A', 'BRAF', 'BTG1', \
'BTG2', 'CARD11', 'CCDC27', 'CD274', 'CD58', 'CD70', 'CD79B', 'CD83', 'CREBBP', 'CRIP1', \
'CXCR4', 'DTX1', 'DUSP2', 'EBF1', 'EEF1A1', 'EP300', 'ETS1', 'ETV6', 'EZH2', 'FADD', 'FAS', \
'GNA13', 'GNAI2', 'GRHPR', 'HIST1H1B', 'HIST1H1C', 'HIST1H1D', 'HIST1H1E', 'HIST1H2AC', \
'HIST1H2AM', 'HIST1H2BC', 'HLA.A', 'HLA.B', 'HLA.C', 'HVCN1', 'IGLL5', 'IKZF3', 'IRF2BP2', \
'IRF4', 'IRF8', 'KLHL6', 'KMT2D', 'KRAS', 'LTB', 'LYN', 'MAP2K1', 'MEF2B', 'MEF2C', 'METAP1D', \
'MYD88.L265P', 'MYD88.OTHER', 'MYD88', 'NFKBIA', 'NFKBIE', 'NOTCH2', 'OSBPL10', 'PABPC1', \
'PIM1', 'PIM2', 'POU2AF1', 'POU2F2', 'PRDM1', 'PTEN', 'PTPN6', 'RAC2', 'RHOA', 'SESN3', \
'SF3B1', 'SGK1', 'SMG7', 'SOCS1', 'SPEN', 'STAT3', 'TBL1XR1', 'TET2', 'TMEM30A', 'TMSB4X', \
'TNFAIP3', 'TNFRSF14', 'TNIP1', 'TOX', 'TP53', 'TUBGCP5', 'UBE2A', 'YY1', 'ZC3H12A', 'ZEB2', \
'ZFP36L1', 'ZNF423']

GENE_DICT = {'HLA-A':'HLA.A','HLA-B':'HLA.B', 'HLA-C':'HLA.C'}
maf['Hugo_Symbol'].replace(GENE_DICT,inplace=True)

if len(args.gsm_genes)>0:
    GENES = pd.read_csv(args.gsm_genes, sep='\t', index_col=0)


# maf limited to GSM genes
maf = maf.loc[(maf['Hugo_Symbol'].isin(GENES))]


# determine blacklist

# make a samples x genes_to_call GSM
GSM = pd.DataFrame(0, index=GENES, columns=sample_set).reset_index()
GSM = GSM.reset_index().set_index('index', drop=False)
GSM.index = GSM.level_0
GSM.drop(columns=['level_0'], inplace=True)
#GSM.reset_index() #inplace=True,drop=False)
GSM.rename(columns={"index":"classifier_name"},inplace=True)
GSM.index = GSM.classifier_name
# Fill in non-silent events first.

non_sil_maf =  maf.loc[~(maf['Variant_Classification'].isin(['Silent']))]

for _, event in non_sil_maf.iterrows():
    gene = event['Hugo_Symbol']
    s = event['Tumor_Sample_Barcode']

    GSM.at[gene, s] = 2

    if gene == 'MYD88':
        if event['Protein_Change'] == 'p.L265P':
            GSM.loc['MYD88.L265P', s] = 2
        else:
            GSM.loc['MYD88.OTHER', s] = 2

# Next call silents. Don't clobber non silent events if there.

sil_maf = maf.loc[maf['Variant_Classification'] == 'Silent']

for _, event in sil_maf.iterrows():
    gene = event['Hugo_Symbol']
    s = event['Tumor_Sample_Barcode']

    if gene == 'MYD88':
        GSM.loc['MYD88.OTHER', s] = 1

    if GSM.loc[gene, s] == 2:
        continue

    GSM.loc[gene, s] = 1



GSM = GSM.round(4)
GSM = GSM.replace(-1, 0.0)

#GSM.index = GSM.index.str.upper()
outfile = args.output_dir + args.id + '.' + TODAY + '.MAF.GSM.tsv'
print('output :', outfile)
GSM.to_csv(outfile, sep='\t',  index=False)
