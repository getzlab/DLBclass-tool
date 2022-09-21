import src.classify_generic as cg
import src.format_data as fd
import pandas as pd
import matplotlib.pyplot as plt

def classify_samples(gsm_file, cohort):
    output_fn = './classifications/' + cohort.replace(' ', '').strip() + '_classified_samples.tsv'
    print('Predictions will be written to:\n\n', output_fn)

    gsm = pd.read_csv(gsm_file, sep='\t', index_col=0)
    reduced_gsm = fd.construct_reduced_winning_version(gsm)

    classified_samples = cg.classify_samples_winning_model(reduced_gsm)
    classified_samples.index.name = 'sample'
    classified_samples['PredictedCluster'] = classified_samples['PredictedCluster'].map({1.0: 'C1', 2.0: 'C2',
                                                                                         3.0: 'C3', 4.0: 'C4',
                                                                                         5.0: 'C5'})
    classified_samples.to_csv(output_fn, sep='\t', index=True)
    return classified_samples


def plot_sample_barplot(sample, classified_samples):
    sample_prediction = classified_samples.loc[sample][['C1', 'C2', 'C3', 'C4', 'C5']]

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.bar(sample_prediction.index, sample_prediction, color=["#803e98", "#00a4d1", "#f29123", "#4d872d", "#dd2929"])
    ax.set_ylabel('Confidence', size=25)
    ax.tick_params('y', labelsize=15)
    ax.set_ylim((0, 1))
    ax.set_title(sample, size=35)
    plt.show()

