import glob
import pandas as pd
import numpy as np
import torch
import src.nn as nn


def classify_samples_winning_model(data, pMax=None, nets=None, modelname=None):
    modelpath = 'saved_models/FINALMODEL_NN_evaluation_seeds1_100_folds5_reducedV3.4_removeN5/*'
    files = glob.glob(modelpath)

    if not nets:
        nets = []
        for file in files:
            loadedNet = torch.load(file)
            NFEATURES = (len(list(loadedNet.items())[1][1][1]))
            net = nn.Net(10, NFEATURES, 5)
            net.load_state_dict(loadedNet)
            net.eval()
            nets.append(net)

    print('loaded neural networks')

    print('Classifying', data.shape[0], 'samples...')

    pMax = 0.93856484

    pred_df = pd.DataFrame()
    i = 0
    for idx, row in data.iterrows():
        net_inputs = torch.tensor(row, dtype=torch.float)
        average_output = None
        for net in nets:
            out = net.forward(net_inputs)
            if average_output is None:
                average_output = out.detach().numpy()
            else:
                average_output = average_output + out.detach().numpy()
        average_output = average_output / len(nets)
        pred_df[idx] = average_output
        i += 1

    pred_df = pred_df.T

    # pMax = pMax
    pred_df_adjusted = pd.DataFrame(columns=['C1', 'C2', 'C3', 'C4', 'C5', 'Confidence', 'PredictedCluster'])
    for idx, row in pred_df.iterrows():
        new_arr = np.float64(np.array(row))
        pWin = np.float64(np.max(row))
        newConf = np.min([np.float64(pWin / pMax), 1])
        shrinkfactor = np.max([np.float64((pMax - pWin) / (pMax * np.float64(np.sum(np.delete(new_arr, new_arr.argmax()))))), 0])
        for idx2 in range(len(new_arr)):
            if idx2 == new_arr.argmax():
                new_arr[idx2] = newConf
            else:
                new_arr[idx2] = new_arr[idx2] * shrinkfactor
        if sum(new_arr) > 1.001 or sum(new_arr) < 0.999:
            print('Sum != 1 (floating point precision), 0.999 < epsilon < 1.001. Renormalizing', new_arr, idx, modelname)
            new_arr = new_arr / sum(new_arr)
        pred_df_adjusted.loc[idx] = np.append(new_arr, [[np.max(new_arr)], [np.argmax(new_arr) + 1]])

    print('Done!')
    return pred_df_adjusted


def classify_samples_generic(data, modelname):
    modelpath = '../saved_models/' + modelname + '/*'
    files = glob.glob(modelpath)

    nets = []
    for file in files:
        loadedNet = torch.load(file)
        NFEATURES = (len(list(loadedNet.items())[1][1][1]))
        net = nn.Net(10, NFEATURES, 5)
        net.load_state_dict(loadedNet)
        net.eval()
        nets.append(net)

    pMax = pd.read_csv('../evaluation_validation_set/' + modelname + '_nfeatures' + str(NFEATURES) + '.tsv',
                       sep='\t', index_col=0).max().max()

    pred_df = pd.DataFrame()
    i = 0
    for idx, row in data.iterrows():
        net_inputs = torch.tensor(row, dtype=torch.float)
        average_output = None
        for net in nets:
            out = net.forward(net_inputs)
            if average_output is None:
                average_output = out.detach().numpy()
            else:
                average_output = average_output + out.detach().numpy()
        average_output = average_output / len(nets)
        pred_df[idx] = average_output
        i += 1

    pred_df = pred_df.T

    # pMax = pMax
    pred_df_adjusted = pd.DataFrame(columns=['C1', 'C2', 'C3', 'C4', 'C5', 'Confidence', 'PredictedCluster'])
    for idx, row in pred_df.iterrows():
        new_arr = np.float64(np.array(row))
        pWin = np.float64(np.max(row))
        newConf = np.min([np.float64(pWin / pMax), 1])
        shrinkfactor = np.max([np.float64((pMax - pWin) / (pMax * np.float64(np.sum(np.delete(new_arr, new_arr.argmax()))))), 0])
        for idx2 in range(len(new_arr)):
            if idx2 == new_arr.argmax():
                new_arr[idx2] = newConf
            else:
                new_arr[idx2] = new_arr[idx2] * shrinkfactor
        if sum(new_arr) > 1.001 or sum(new_arr) < 0.999:
            print('Sum != 1 (floating point precision), 0.999 < epsilon < 1.001. Renormalizing', new_arr, idx, modelname)
            new_arr = new_arr / sum(new_arr)
        pred_df_adjusted.loc[idx] = np.append(new_arr, [[np.max(new_arr)], [np.argmax(new_arr) + 1]])

    return pred_df_adjusted
