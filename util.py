import numpy as np
from numpy.linalg import inv, pinv
import nibabel as nb

def zstandarize_ts(X):
    X = X - X.mean(axis = 0, keepdims = True)
    X = X / np.sqrt(np.nansum(X**2, axis=0)/X.shape[0])
    return X

def correlate(X, Y):
    """ Correlate X and Y numpy arrays after standardizing them"""
    X = zstandarize_ts(X)
    Y = zstandarize_ts(Y)
    return Y.T @ X / X.shape[0]
