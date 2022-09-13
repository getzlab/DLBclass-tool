# Import required libraries
import pandas as pd
import src.classify_generic as cg
import src.format_data as fd

# Read in GSM

# *** Replace filename with your GSM file to classify your own data ***

gsm_file = '../gsm/DLBCL_testset_gsm.tsv'

gsm = pd.read_csv(gsm_file, sep='\t', index_col=0)
gsm.head()

reduced_gsm = fd.construct_reduced_winning_version(gsm)
reduced_gsm.head()

classified_samples = cg.classify_samples_winning_model(reduced_gsm)


labels = pd.read_csv('/Users/twood/Desktop/DLBCL-Classifier-Public/data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.tsv',
                    sep='\t', index_col=0)
labels = labels[labels.index.isin(classified_samples.index)]
labels = labels.loc[classified_samples.index]

print(sum(labels['cluster'] == classified_samples['PredictedCluster']) / 149)