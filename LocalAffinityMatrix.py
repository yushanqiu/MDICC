import numpy as np
import pandas as pd


def kscale(matrix, k=7, dists=None, minval=0.004):
    """ Returns the local scale based on the k-th nearest neighbour """
    r,c = matrix.shape
    scale = np.zeros((r,c))
    for i in range(1,(k+1)):
        ix = (np.arange(len(matrix)), matrix.argsort(axis=0)[i])
        d = matrix[ix][np.newaxis].T
        dists = (d if dists is None else dists)
        scale1 = dists.dot(dists.T)
        scale = scale1 + scale
    scale = scale/k
    return np.clip(scale, minval, np.inf)


def affinity(matrix, k):
    scale = kscale(matrix, k)
    msq = matrix * matrix
    scaled = -msq /(0.5*scale+0.5*matrix)
    scaled[np.where(np.isnan(scaled))] = 0.0
    a = np.exp(scaled)
    a.flat[::matrix.shape[0]+1] = 0.0  # zero out the diagonal
    return a


def testaff(matrix,k):
    k = int(k)
    a = matrix
    affi = affinity(a,k)
    return affi


# dists = pd.read_csv(r'D:\AS_RBP\数据去噪\a.csv',index_col = 0)
# dists = dists.values
#
# aff = affinity(dists, 2)
# print(aff)

