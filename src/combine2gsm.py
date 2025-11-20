import pandas as pd
import numpy as np
import argparse
from datetime import datetime
import time

start = time.time()
TODAY = datetime.now().strftime("%d%b%Y")

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--id", help="Cohort set name.", required=True, type=str)
parser.add_argument("-v", "--sv_gsm", help="SV GSM.", required=False, type=str)
parser.add_argument("-m", "--mutation_gsm", help="MAF GSM.", required=True, type=str)
parser.add_argument(
    "-c", "--copy_number_gsm", help="SEG GSM.", required=False, type=str
)
parser.add_argument(
    "-p",
    "--published_gsm",
    help="GSM used to provide average values of missing features .",
    required=True,
    type=str,
    default="gsm/DLBclass_published.GSM.tsv",
)
parser.add_argument(
    "-f", "--feature_order_file", help="GSM Feature order .", required=True, type=str
)
parser.add_argument(
    "-o",
    "--output_dir",
    help="Output directory.",
    required=False,
    type=str,
    default="./",
)
args = parser.parse_args()

# Load and sort mutation data frame
mut_df = pd.read_csv(
    args.mutation_gsm, sep="\t", index_col=0
)  #'../../data_tables/gsm/DLBCL.699.mutationsGSM.Sep_23_2022.tsv'

sample_order = list(mut_df.columns)
sample_order.sort()

mut_df = mut_df.reindex(columns=sample_order)

full_gsm = mut_df.copy()

# Load and sort copy number data frame if provided
if args.copy_number_gsm:
    scna_df = pd.read_csv(
        args.copy_number_gsm, sep="\t", index_col=0
    )  #'../../data_tables/gsm/DLBCL.699.scnaGSM.Sep_23_2022.tsv'
    col = scna_df.index.values.tolist()
    # prefix with "X" if not already prefixed
    if not (col[0][0] == "X"):
        xcol = ["X" + x for x in col]
        scna_df.index = xcol
    scna_df = scna_df.reindex(columns=sample_order)
    full_gsm = pd.concat([full_gsm, scna_df])

# Load and sort structural variant data frame if provided
if args.sv_gsm:
    sv_df = pd.read_csv(
        args.sv_gsm, sep="\t", index_col=0
    )  #'../../data_tables/gsm/DLBCL.699.svGSM.Sep_23_2022.tsv'
    sv_df = sv_df.reindex(columns=sample_order)
    full_gsm = pd.concat([full_gsm, sv_df])

# Load feature order
feature_order = list(
    pd.read_csv(args.feature_order_file, index_col=0, header=None).index
)  #'../../data_tables/gsm/feature_order.tsv'

# Load published GSM matrix from training data
# Calculate average of each feature across samples
# This code can be reimplemented if you want to substitue 
# missing SV values with a 1 (the only features with a non-zero average)
# pub_gsm = pd.read_csv(args.published_gsm, sep='\t', index_col=0)
# pub_gsm = pd.DataFrame(pub_gsm.mean(axis=1).round().astype("Int64"))
# for sample in sample_order:
#     if sample not in pub_gsm.columns:
#         pub_gsm[sample] = pub_gsm[0]
# pub_gsm = pub_gsm.reindex(index = feature_order, columns = sample_order)
# print(pub_gsm)

# Fill missing values in full_gsm with zeros
# The paper says missing features should be replaced with the average 
# from the published gsm matrix, but that value is 0 for all features 
# except for the SVs, which get an average of 1 because they are 
# assigned a value of 3 when present. So, for consistency across all 
# feature types, missing data will be given a 0 value. 
full_gsm = full_gsm.reindex(index=feature_order, columns=sample_order)
full_gsm = full_gsm.fillna(0)

# Ensure the correct feature order and that every value is an integer
gsm = full_gsm.copy()
gsm = gsm.loc[feature_order, sample_order]
gsm = gsm.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
gsm = gsm.reset_index().copy()
gsm = gsm.rename(columns={"index": "classifier_name"}).copy()

outfile = args.output_dir + "/" + args.id + "." + TODAY + ".GSM.tsv"
print("output :", outfile)
gsm.to_csv(outfile, sep="\t", index=False)
