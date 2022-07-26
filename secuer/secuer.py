#!/usr/secuer_console/env python
# coding: utf-8
import os
from pathlib import Path
import yaml
from scipy.linalg import eigh
from scipy.sparse import issparse, csr_matrix
from scipy import spatial
import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import time
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.cluster import KMeans,DBSCAN,AgglomerativeClustering
import gc
import igraph as ig
import louvain
warnings.filterwarnings("ignore")
import numba
from numba import njit, prange
numba.config.NUMBA_DEFAULT_NUM_THREADS=8

# __all__ = ["secuer","Read", "Tcut_for_bipartite_graph"]

def Read(filename,sheet = None,istranspose=False):
    '''  
    Read 10x-Genomics-formatted hdf5 file.
    '''
    try:
        data = sc.read(filename,sheet=sheet)
        if istranspose:
            data=data.T
    except:
        try:
            data = sc.read_10x_h5(filename)
        except:
            data = sc.read_10x_mtx(filename)
    # if type=='10x_h5':
    #     data = sc.read_10x_h5(filename)
    # elif type == '10x_mtx':
    #     data = sc.read_10x_mtx(filename)
    # else:
    #     data = sc.read(filename)
    return data

# This code is from Scanpy
def get_indices_distance_from_dense_matrix(D, neighSize: int, returnDis=True):
    '''
    :param D: a distance matrix with
    :param neighSize:
    :param returnDis:
    :return:
    '''
    sample_range = np.arange(D.shape[0])[:, None]
    indices = np.argpartition(D, neighSize, axis=1)[:, :neighSize]
    indices = indices[sample_range, np.argsort(D[sample_range, indices])]
    if returnDis:
        distances = D[sample_range, indices]
        return indices, distances
    else:
        return indices
        
def pdist2_fast(A, B, metric='euclidean'):
    '''
    compute distabce between A and B with different method
    A: 2D array
    B: 2D array
    euclidean: ['sqeuclidean','euclidean','L1','cosine']
    return: dist matrix of A and B
    '''
    if metric == 'L1':
        res = spatial.distance.cdist(A, B, metric='minkowski', p=1, V=None,
                                     VI=None, w=None)
    elif metric == 'sqeuclidean':
        res = spatial.distance.cdist(A, B, metric='sqeuclidean', p=None, V=None,
                                     VI=None, w=None)
    elif metric == 'euclidean':
        res = spatial.distance.cdist(A, B, metric='euclidean', p=None, V=None,
                                     VI=None, w=None)
    elif metric == 'cosine':
        res = spatial.distance.cdist(A, B, metric='cosine', p=None, V=None,
                                     VI=None, w=None)
    else:
        print('unknown metric')
    return res

def fast_kmeans_scipy(ds, k, max_iter=20):
    '''
    :param ds: Data used for clustering.
    :param k: The number of clustering.
    :param max_iter: Maximum number of iterations.
    '''
    m, n = ds.shape
    result = np.empty(m, dtype=int)
    np.random.seed(1)
    cores = ds[np.random.choice(m, k, replace=False)]
    maxiter: int = 1
    while True:
        if (maxiter != 1 and len(np.unique(result)) != k):
            np.random.seed(2)
            cores = ds[np.random.choice(m, k, replace=False)]
        distance = pdist2_fast(ds, cores)
        index_min = np.argmin(distance, axis=1)
        if (index_min == result).all():
            return result, cores
        if (maxiter >= max_iter):
            return result, cores
        result[:] = index_min
        for i in range(k):
            items = ds[result == i]
            cores[i] = np.mean(items, axis=0)
        maxiter += 1
        
def getRepresentativesByHybridSelection(fea, pSize, cntTimes=10, seed=1):
    '''
    Select pSize representatives by hybrid selection.
    :param fea: Input data.
    :param pSize: The number of representatives.
    :return:The pSize cluster centers as the final representatives.
    '''
    N = fea.shape[0]
    # randomly select pSize * cntTimes candidate representatives.
    bigPSize = cntTimes * pSize
    if pSize > N:
        pSize = N
    if bigPSize > N:
        bigPSize = N
    np.random.seed(seed)
    selectIdxs = np.random.choice(N, bigPSize, replace=False)
    bigRpFea = fea[selectIdxs, :]
    label, RpFea = fast_kmeans_scipy(bigRpFea, pSize)  # max_iter=20
    return RpFea

def partition(nums, start, end, k):
    """
    During the process, it's guaranteed start <= k <= end
    """
    if start == end:
        return nums[k]

    left, right = start, end
    pivot = nums[(start + end) // 2]
    while left <= right:
        while left <= right and nums[left] < pivot:
            left += 1
        while left <= right and nums[right] > pivot:
            right -= 1
        if left <= right:
            nums[left], nums[right] = nums[right], nums[left]
            left, right = left + 1, right - 1

    # left is not bigger than right
    if k <= right:
        return partition(nums, start, right, k)
    if k >= left:
        return partition(nums, left, end, k)
    return nums[k]
    
def kthLargestElement(k, A):
    '''Find the kth largest number'''
    if not A or k < 1 or k > len(A):
        return None
    return partition(A, 0, len(A) - 1, len(A) - k)

def Estimatekbyeigen(vector, gapth=4):
    '''Estimate the number of k by eigenvectors of bipartite graph'''
    vecbin = pd.cut(np.abs(vector) ** 0.0001, 1000).value_counts()
    diff = np.argwhere(vecbin.values != 0).squeeze()[1:] - np.argwhere(vecbin.values != 0).squeeze()[0:-1]
    thre = kthLargestElement(gapth, diff.tolist())
    n = np.argwhere(diff >= thre)[-1] + 1
    Ks = sum(np.abs(vector) ** 0.0001 < vecbin.index[np.argwhere(vecbin.values != 0).squeeze()[n]][0].left)
    return Ks

def Tcut_for_bipartite_graph(B, Nseg,eskMethod, gapth=4,maxKmIters=100,
                             cntReps=3,clusterMethod='Kmeans',
                             eps=0.1, min_samples=100):
    '''
    Bipartite graph partition by transfer-cuts 
    '''
    Nx, Ny = B.shape
    # if Ny < Nseg:
    #     print('Need more columns!')
    dx = np.sum(B, axis=1)
    dx[dx == 0] = 1e-10  # Just to make 1./dx feasible.
    dx = dx.A1
    Dx = csr_matrix((1 / dx, np.arange(Nx), np.arange(Nx + 1)), shape=(Nx, Nx))
    Dx.eliminate_zeros()
    Wy = B.T @ Dx @ B
    d = np.sum(Wy, axis=1).A1
    d[d <= 0] = 1e-10
    D = csr_matrix((1 / np.sqrt(d), np.arange(Ny), np.arange(Ny + 1)), shape=(Ny, Ny))
    D.eliminate_zeros()
    nWy = D @ Wy @ D
    nWy = (nWy + nWy.T) / 2
    # computer eigenvectors
    eval_, evec = eigh(nWy.toarray())
    eval_, evec = np.real(eval_), np.real(evec)
    if eskMethod != 'subGraph' and not Nseg:
        Nseg = Estimatekbyeigen(-1+2*np.sqrt(eval_),gapth)
    idx = np.argsort(-eval_)
    Ncut_evec = D @ evec[:, idx[0:Nseg]]
    # compute the Ncut eigenvectors on the entire bipartite graph (transfer!)
    evec = Dx @ B @ Ncut_evec
    # normalize each row to unit norm
    evec = evec / (np.sqrt(np.sum(evec ** 2, axis=1)) + 1e-10).reshape([-1, 1])
    # k-means
    if clusterMethod == 'Kmeans':
        lab = KMeans(n_clusters=Nseg, max_iter=maxKmIters, n_init=cntReps,
                     init='k-means++').fit(evec)
    elif clusterMethod == 'DBSCAN':
        lab = DBSCAN(eps=eps, min_samples=min_samples).fit(evec)
    else:
        lab = AgglomerativeClustering(n_clusters=Nseg).fit(evec)
    return lab.labels_,Nseg

def EstimatekbysubGraph(RpfeaDist, RpFeaKnnIdx, resolution, addweights,Knn=5):
    '''
    Estimated the clustering by graph built by anchors
    '''
    # RpFeaKnnIdx5 = RpFeaKnnIdx[:,0:Knn]
    index = np.array([[i] * Knn for i in range(RpFeaKnnIdx.shape[0])]).flatten()
    colunm = RpFeaKnnIdx[:, 0:Knn].flatten()
    weights = RpfeaDist[index, colunm]
    G = ig.Graph()
    G.add_vertices(RpFeaKnnIdx.shape[0])  # this adds adjacency.shape[0] vertices
    G.add_edges(list(zip(index, colunm)))
    G.es['weight'] = weights
    if addweights == True:
        part_new = louvain.find_partition(G, louvain.RBConfigurationVertexPartition,
                                          seed=0, resolution_parameter=resolution, weights=weights)
    else:
        part_new = louvain.find_partition(G, louvain.RBConfigurationVertexPartition,
                                          seed=0, resolution_parameter=resolution)
    return np.unique(part_new.membership).shape[0]

@njit(cache=True, parallel=True)
def NearestRepIndex(fea,RpFea,N,p,
         repClsLabel,centerDist,
         distance):
    cntRepCls = int(p ** 0.5)
    minCenterIdxs = np.argmin(centerDist, axis=1)
    nearestRepInRpFeaIdx = np.empty(N, dtype='int64')
    for i in prange(cntRepCls):
        tmp = np.where(repClsLabel == i)[0]
        nearestRepInRpFeaIdx[minCenterIdxs == i] = np.argmin(pdist2_fast(fea[minCenterIdxs == i, :],
                                                                         RpFea[tmp, :], metric=distance), axis=1)
        nearestRepInRpFeaIdx[minCenterIdxs == i] = tmp[nearestRepInRpFeaIdx[minCenterIdxs == i]]
    return nearestRepInRpFeaIdx

# adj_param=False can remove
def secuer(fea, Ks=None,
          distance='euclidean',
          p=1000,
          Knn=5,
          mode='secuer',
          eskMethod = 'subGraph',
          eskResolution=0.8,
          addweights=False,
          seed=1,
          gapth=4,
          clusterMethod='Kmeans',
          maxTcutKmIters=100,
          cntTcutKmReps=3):
    N = fea.shape[0]  # n*F
    if p > N:
        p = N
    # print(p)
    # Get $p$ representatives by hybrid selection
    RpFea = getRepresentativesByHybridSelection(fea, p,seed=seed)
    # Approx.KNN
    # 1.partition RpFea into $cntRepCls$ rep - clusters
    cntRepCls = int(p ** 0.5)
    # 2. find the center of each rep-cluster
    if distance == 'euclidean':
        repClsLabel, repClsCenters = fast_kmeans_scipy(RpFea, cntRepCls)  # max_iter=20, n_init=1
    else:
        repClsLabel, repClsCenters = fast_kmeans_scipy(RpFea, cntRepCls)
    #  3. Pre-compute the distance between N objects and the $cntRepCls$
    centerDist = pdist2_fast(fea, repClsCenters, metric=distance)
    del repClsCenters
    gc.collect()
    minCenterIdxs = np.argmin(centerDist, axis=1)
    nearestRepInRpFeaIdx = np.empty(N, dtype='int64')
    for i in range(cntRepCls):
        tmp = np.where(repClsLabel == i)[0]
        nearestRepInRpFeaIdx[minCenterIdxs == i] = np.argmin(pdist2_fast(fea[minCenterIdxs == i, :],
                                                                         RpFea[tmp, :], metric=distance), axis=1)
        nearestRepInRpFeaIdx[minCenterIdxs == i] = tmp[nearestRepInRpFeaIdx[minCenterIdxs == i]]
    del repClsLabel, minCenterIdxs, tmp
    gc.collect()
    # For each object, compute its distance to the candidate neighborhood of its nearest representative (in RpFea)
    neighSize = 10 * Knn
    RpfeaDist = pdist2_fast(RpFea, RpFea, metric=distance)
    if neighSize > RpfeaDist.shape[0]:
        neighSize = RpfeaDist.shape[0]-2
    RpFeaKnnIdx = get_indices_distance_from_dense_matrix(RpfeaDist,
                                                         neighSize + 1,
                                                         returnDis=False)  # p * neighSize+1
    # estimate k
    if not Ks and eskMethod=='subGraph':
        # start=time.time()
        Ks = EstimatekbysubGraph(RpfeaDist, RpFeaKnnIdx, eskResolution,
                                 addweights,Knn=Knn)
    del RpfeaDist
    gc.collect()
    RpFeaKnnDist = np.zeros([N, np.shape(RpFeaKnnIdx)[1]])
    for i in range(p):
        RpFeaKnnDist[nearestRepInRpFeaIdx == i, :] = pdist2_fast(fea[nearestRepInRpFeaIdx == i, :],
                                                                 RpFea[RpFeaKnnIdx[i, :], :],
                                                                 metric=distance)
    # Get the final KNN according to the candidate neighborhood.
    RpFeaKnnIdxFull = RpFeaKnnIdx[nearestRepInRpFeaIdx, :]  # 得到每个样本对应的最近的10*knn+1个邻居
    knnIdx1, knnDist = get_indices_distance_from_dense_matrix(RpFeaKnnDist, Knn)  #
    knnIdx = RpFeaKnnIdxFull[np.arange(N)[:, None], knnIdx1]
    del knnIdx1, RpFeaKnnIdxFull, RpFeaKnnDist
    gc.collect()
    # Compute the cross-affinity matrix B for the bipartite graph
    if distance == 'cosine':
        Gsdx = 1 - knnDist
    else:
        knnMeanDiff = np.mean(knnDist, 1)  # use the mean distance of anchors as the kernel parameter $\sigma_i$ for each cell 
        knnMeanDiff[knnMeanDiff==0] = np.mean(knnDist)
        Gsdx = np.exp(-(knnDist ** 2) / (2 * knnMeanDiff[:,None] ** 2))
    Gsdx[Gsdx == 0] = 1e-16
    indptr = np.arange(0, N * Knn + 1, Knn)
    B = csr_matrix((Gsdx.copy().ravel(),  # copy the data, otherwise strange behavior here
                    knnIdx.copy().ravel(), indptr), shape=(N, p))
    try:
        labels, Ks = Tcut_for_bipartite_graph(B, Ks, eskMethod,gapth,
                                              maxTcutKmIters,
                                              cntTcutKmReps,
                                              clusterMethod=clusterMethod)
    except:
        labels, Ks = Tcut_for_bipartite_graph(B, Ks[0], eskMethod,gapth,
                                              maxTcutKmIters, cntTcutKmReps,
                                              clusterMethod=clusterMethod)
    if mode=='secuer':
        return labels
    else:
        return labels, Ks

if __name__ == '__main__':
    file = 'D://My_data//Allproject//Secuer//Clustering0804//gold_label_data/'
    files = os.listdir(file)
    fileh5sd_gold = [i for i in files if i.endswith('.h5ad')]
    nmi_secuer=[]
    ari_secuer=[]
    t_secuer=[]
    for i in range(len(fileh5sd_gold)):
        print(fileh5sd_gold[i])
        data = sc.read(file + fileh5sd_gold[i])
        fea = data.obsm['X_pca']
        # print(f'{mtx[j]}: secuer')
        start = time.time()
        # sc.pp.neighbors(data)
        # sc.tl.louvain(data)
        # res = data.obs['louvain']
        res = secuer(fea, p=1000, Knn=7, maxTcutKmIters=100,
                    cntTcutKmReps=3,seed=1)
        end = time.time() - start
        print(np.unique(res).shape[0],np.unique(data.obs['celltype']).shape[0])

        name_up = 'secuer_lable_' + str(i)
        data.obs[name_up] = pd.Categorical(res)
        nmi_secuer.append(normalized_mutual_info_score(data.obs['celltype'], res))
        ari_secuer.append(adjusted_rand_score(data.obs['celltype'], res))
        t_secuer.append(end)
        print()
    res_secuer = pd.DataFrame([t_secuer, nmi_secuer, ari_secuer],
                             index=['t_secuer', 'nmi_secuer', 'ari_secuer'],
                             columns = fileh5sd_gold)
    print(res_secuer.values)
if __name__ == '__main__':
    file = 'D://My_data//Allproject//Secuer//Clustering0804//gold_label_data/'
    files = os.listdir(file)
    fileh5sd_gold = [i for i in files if i.endswith('.h5ad')]
    nmi_Secuer=[]
    ari_Secuer=[]
    t_Secuer=[]
    for i in range(len(fileh5sd_gold)):
        print(fileh5sd_gold[i])
        data = sc.read(file + fileh5sd_gold[i])
        fea = data.obsm['X_pca']
        # print(f'{mtx[j]}: Secuer')
        start = time.time()
        # sc.pp.neighbors(data)
        # sc.tl.louvain(data)
        # res = data.obs['louvain']
        res = secuer(fea, p=1000, Knn=7, maxTcutKmIters=100,
                    cntTcutKmReps=3,seed=1)
        end = time.time() - start
        print(np.unique(res).shape[0],np.unique(data.obs['celltype']).shape[0])

        name_up = 'Secuer_lable_' + str(i)
        data.obs[name_up] = pd.Categorical(res)
        nmi_Secuer.append(normalized_mutual_info_score(data.obs['celltype'], res))
        ari_Secuer.append(adjusted_rand_score(data.obs['celltype'], res))
        t_Secuer.append(end)
        print()
    res_Secuer = pd.DataFrame([t_Secuer, nmi_Secuer, ari_Secuer],
                             index=['t_Secuer', 'nmi_Secuer', 'ari_Secuer'],
                             columns = fileh5sd_gold)
    print(res_Secuer.values)
