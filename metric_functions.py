import numpy as np
from scipy.stats import pearsonr

def get_correlations(x1, x2, rowwise = True):

    if rowwise:
        cors = [pearsonr(x1[idx,], x2[idx,])[0] for idx in range(x1.shape[0])]
    else:
        cors = [pearsonr(x1[:,idx], x2[:,idx])[0] for idx in range(x1.shape[1])]

    return(cors)


def get_mse(A, B, rowwise = True):

    if rowwise:
        ax = 0
        mse = ((A - B)**2).mean(axis=ax)
    else:
        ax = 1
        mse = ((A - B)**2).mean(axis=ax)

    return(mse)


def get_mae(A, B, rowwise = True):

    if rowwise:
        ax = 0
        mae = np.abs(A - B).mean(axis=ax)
    else:
        ax = 1
        mae = np.abs(A - B).mean(axis=ax)

    return(mae)
        

def matrix_cor(x, y):

    return np.corrcoef(np.corrcoef(x, rowvar = False).flatten(), np.corrcoef(y, rowvar = False).flatten())[0,1]


def matrix_mse(x, y):
    
    mses = get_mse(np.corrcoef(x, rowvar = False).flatten(), np.corrcoef(y, rowvar = False).flatten())
    return np.mean(mses)

def matrix_mae(x, y):
    
    maes = get_mae(np.corrcoef(x, rowvar = False).flatten(), np.corrcoef(y, rowvar = False).flatten())
    return np.mean(maes)
    