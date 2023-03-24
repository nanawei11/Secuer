#!/usr/secuer_console/env python
# coding: utf-8
# from pathlib import Path
# import yaml
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
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
import gc
import igraph as ig
import louvain

warnings.filterwarnings("ignore")
import numba
from multiprocessing.pool import ThreadPool as Pool
from functools import partial
import logging


# numba.config.NUMBA_DEFAULT_NUM_THREADS=8

# __all__ = ["secuer","Read", "Tcut_for_bipartite_graph"]

class Logger:
    def __init__(self):
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)
        self.logger.handlers.clear()
        fmt = logging.Formatter('[%(asctime)s] [%(levelname)s] %(message)s', '%Y-%m-%d %H:%M:%S')
        logg = logging.StreamHandler()
        logg.setFormatter(fmt)
        self.logger.addHandler(logg)

    def info(self, message):
        self.logger.info(message)

    def war(self, message):
        self.logger.warn(message)

    def error(self, message):
        self.logger.error(message)


def Read(filename, sheet=None, istranspose=False):
    '''  
    Read 10x-Genomics-formatted hdf5 file.
    '''
    try:
        data = sc.read(filename, sheet=sheet)
        if istranspose:
            data = data.T
    except:
        try:
            data = sc.read_10x_h5(filename)
        except:
            data = sc.read_10x_mtx(filename)
    return data


# This code is from Scanpy
def get_indices_distance_from_dense_matrix(D, neighSize: int, returnDis=True):
    '''
    Get the k-nearest neighbors from a distance matrix.
    :param D: A distance matrix.
    :param neighSize: The number of neighbors.
    :param returnDis: whether to return a distance matrix of k neighbors.
    :return: The neighbors' distances or neighbors' indices.
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
    Compute distance between A and B with different methods.
    :param A: 2D array
    :param B: 2D array
    :param euclidean: ['sqeuclidean','euclidean','L1','cosine']
    return: The distance between the corresponding rows of A and B.
    '''
    if metric == 'L1':
        res = spatial.distance.cdist(A, B, metric='minkowski', p=1)
    elif metric == 'sqeuclidean':
        res = spatial.distance.cdist(A, B, metric='sqeuclidean')
    elif metric == 'euclidean':
        res = spatial.distance.cdist(A, B, metric='euclidean')
    elif metric == 'cosine':
        res = spatial.distance.cdist(A, B, metric='cosine')
    else:
        print('unknown metric')
    return res


def fast_kmeans_scipy(ds, k, max_iter=20):
    '''
    Fast implementation of kmeans.
    :param ds: Data used for clustering.
    :param k: The number of clustering.
    :param max_iter: Maximum number of iterations.
    :return The label of each data points.
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
    :return: The pSize cluster centers as the final representatives.
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
    '''
    Find the kth largest number.
    '''
    if not A or k < 1 or k > len(A):
        return None
    return partition(A, 0, len(A) - 1, len(A) - k)


def Estimatekbyeigen(vector, gapth=4):
    '''
    Estimate the number of k by eigenvectors of bipartite graph
    '''
    vecbin = pd.cut(np.abs(vector) ** 0.0001, 1000).value_counts()
    diff = np.argwhere(vecbin.values != 0).squeeze()[1:] - np.argwhere(vecbin.values != 0).squeeze()[0:-1]
    thre = kthLargestElement(gapth, diff.tolist())
    n = np.argwhere(diff >= thre)[-1] + 1
    Ks = sum(np.abs(vector) ** 0.0001 < vecbin.index[np.argwhere(vecbin.values != 0).squeeze()[n]][0].left)
    return Ks


def Tcut_for_bipartite_graph(B, Nseg, eskMethod, gapth=4, maxKmIters=100,
                             cntReps=3, clusterMethod='Kmeans',
                             eps=0.1, min_samples=100):
    '''
    Bipartite graph partition by transfer-cuts.
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
        Nseg = Estimatekbyeigen(-1 + 2 * np.sqrt(eval_), gapth)
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
    elif clusterMethod == 'AgglomerativeClustering':
        lab = AgglomerativeClustering(n_clusters=Nseg).fit(evec)
    else:
        print('Unknown clustering method. please choose one of [kmeans, DBSCAN, AgglomerativeClustering].')
    return lab.labels_, Nseg


def EstimatekbysubGraph(RpfeaDist, RpFeaKnnIdx, resolution, addweights, Knn=5):
    '''
    Estimated the clustering by graph built by anchors.
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


def Multi_process_cntRepCls(fea, repClsLabel_m, nearestRepInRpFeaIdx_m,
                            minCenterIdxs_m, RpFea_m, distance, i):
    # qq = i
    tmp = np.where(repClsLabel_m == i)[0]
    nearestRepInRpFeaIdx_m[minCenterIdxs_m == i] = np.argmin(pdist2_fast(fea[minCenterIdxs_m == i, :],
                                                                         RpFea_m[tmp, :], metric=distance), axis=1)
    pp = tmp[nearestRepInRpFeaIdx_m[minCenterIdxs_m == i]]
    return pp


def Multi_RpFeaKnnDist(fea, nearestRepInRpFeaIdx, RpFea, RpFeaKnnIdx, distance, i):
    pp = pdist2_fast(fea[nearestRepInRpFeaIdx == i, :],
                     RpFea[RpFeaKnnIdx[i, :], :],
                     metric=distance)
    return pp


def secuer(fea,
           Ks=None,
           distance='euclidean',
           p=1000,
           Knn=5,
           mode='secuer',
           eskMethod='subGraph',
           eskResolution=0.8,
           addweights=False,
           seed=1,
           gapth=4,
           clusterMethod='Kmeans',
           Gaussiankernel='localscaled',
           maxTcutKmIters=100,
           cntTcutKmReps=3,
           multiProcessState=False,
           num_multiProcesses=4
           ):
    '''
    Cluster cell using Secuer.
    :param fea: A expression matrix with cells by genes.
    :param distance: str, default=euclidean. The metrics of distance. Can be sqeuclidean, euclidean, L1, cosine.
    :param p: int, default=1000. The number of anchors.
    :param Knn: int, default=7. The number of neighbors of anchors for each cell.
    :param mode: str, default=secuer. secuer for performing secuer or secuer-C when using secuerconsensus function.
    :param eskMethod: str, default=subGraph. The methods to estimate the number of clusters: subGraph or bipart-eigen.
    :param eskResolution: float, default=0.8. The resolution for estimating the clusters when using subGraph.The higher the resolution, the more clusters.
    :param gapth: int, default=4. The gap values of estimating the clusters when using bipart-eigen.The higher the gapth, the more clusters.
    :param clusterMethod: str, default=Kmeans. The final clustering methods in secuer. Can be Kmeans, DBSCAN or AgglomerativeClustering.
    :param Gaussiankernel: str, default=localscaled. The gaussian kernel used to build bipartite graph. Can be localscaled or traditionalscaled.
    :param maxTcutKmIters: int, default=100. The maximum number of iterations for k-means.
    :param cntTcutKmReps: int, default=3. Number of times with different centroid seeds in the k-means.
    :param multiProcessState: bool, default=False. Whether to use parallel. Recommend to default by False.
    :param num_multiProcess: int, default=4. The number of parallel processes.
    :return: 1D-array. The labels for each cell.
    '''
    N = fea.shape[0]  # n*F
    if p > N:
        p = N

    logg = Logger()

    # print(p)
    # Get $p$ representatives by hybrid selection
    if mode == 'secuer':
        logg.info('Selecting representatives...')
    RpFea = getRepresentativesByHybridSelection(fea, p, seed=seed)
    # Approx.KNN
    # 1.partition RpFea into $cntRepCls$ rep - clusters
    if mode == 'secuer':
        logg.info('Approximate KNN...')
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

    #### parallel
    if multiProcessState == False:
        for i in range(cntRepCls):
            tmp = np.where(repClsLabel == i)[0]
            nearestRepInRpFeaIdx[minCenterIdxs == i] = np.argmin(pdist2_fast(fea[minCenterIdxs == i, :],
                                                                             RpFea[tmp, :], metric=distance), axis=1)
            nearestRepInRpFeaIdx[minCenterIdxs == i] = tmp[nearestRepInRpFeaIdx[minCenterIdxs == i]]
        del repClsLabel, minCenterIdxs, tmp

    else:
        pool = Pool(num_multiProcesses)
        func = partial(Multi_process_cntRepCls, fea, repClsLabel, nearestRepInRpFeaIdx, minCenterIdxs, RpFea, distance)
        outputs = pool.map(func, np.arange(0, cntRepCls, 1))
        for i in range(cntRepCls):
            nearestRepInRpFeaIdx[minCenterIdxs == i] = outputs[i]
        del outputs
    gc.collect()
    # For each object, compute its distance to the candidate neighborhood of its nearest representative (in RpFea)
    neighSize = 10 * Knn
    RpfeaDist = pdist2_fast(RpFea, RpFea, metric=distance)
    if neighSize > RpfeaDist.shape[0]:
        neighSize = RpfeaDist.shape[0] - 2
    RpFeaKnnIdx = get_indices_distance_from_dense_matrix(RpfeaDist,
                                                         neighSize + 1,
                                                         returnDis=False)  # p * neighSize+1

    # estimate k
    if not Ks and eskMethod == 'subGraph':
        if mode == 'secuer':
            logg.info('Estimating the number of clustering...')
        # start=time.time()
        Ks = EstimatekbysubGraph(RpfeaDist, RpFeaKnnIdx, eskResolution,
                                 addweights, Knn=Knn)
    del RpfeaDist
    gc.collect()
    RpFeaKnnDist = np.zeros([N, np.shape(RpFeaKnnIdx)[1]])

    # new add
    if multiProcessState == False:
        for i in range(p):
            RpFeaKnnDist[nearestRepInRpFeaIdx == i, :] = pdist2_fast(fea[nearestRepInRpFeaIdx == i, :],
                                                                     RpFea[RpFeaKnnIdx[i, :], :],
                                                                     metric=distance)
    else:
        pool = Pool(num_multiProcesses)
        RpFeaKnnDist_func = partial(Multi_RpFeaKnnDist,
                                    fea, nearestRepInRpFeaIdx, RpFea, RpFeaKnnIdx, distance)
        outputs = pool.map(RpFeaKnnDist_func, np.arange(0, p, 1))
        for i in range(p):
            #### Calculate distance between neighbors and raw data
            RpFeaKnnDist[nearestRepInRpFeaIdx == i, :] = outputs[i]
    # Get the final KNN according to the candidate neighborhood.
    # Get the nearest 10*knn+1 neighbors corresponding to each sample
    RpFeaKnnIdxFull = RpFeaKnnIdx[nearestRepInRpFeaIdx, :]
    knnIdx1, knnDist = get_indices_distance_from_dense_matrix(RpFeaKnnDist, Knn)  #
    knnIdx = RpFeaKnnIdxFull[np.arange(N)[:, None], knnIdx1]
    del knnIdx1, RpFeaKnnIdxFull, RpFeaKnnDist
    gc.collect()
    if mode == 'secuer':
        logg.info('Bipartite graph partitioning...')
    # Compute the cross-affinity matrix B for the bipartite graph
    if distance == 'cosine':
        Gsdx = 1 - knnDist
    else:
        if Gaussiankernel == 'localscaled':
            knnMeanDiff = np.mean(knnDist,
                                  1)  # use the mean distance of anchors as the kernel parameter $\sigma_i$ for each cell
            knnMeanDiff[knnMeanDiff == 0] = np.mean(knnDist)
            Gsdx = np.exp(-(knnDist ** 2) / (2 * knnMeanDiff[:, None] ** 2))
        elif Gaussiankernel == 'traditionalscaled':
            knnMeanDiff = np.mean(knnDist)  # use the mean distance as the kernel parameter $\sigma$
            Gsdx = np.exp(-(knnDist ** 2) / (2 * knnMeanDiff ** 2))

        else:
            print(f'Unknown Gaussiankernel. please choose one of [localscaled, traditionalscaled].')
    Gsdx[Gsdx == 0] = 1e-16
    indptr = np.arange(0, N * Knn + 1, Knn)
    B = csr_matrix((Gsdx.copy().ravel(),  # copy the data, otherwise strange behavior here
                    knnIdx.copy().ravel(), indptr), shape=(N, p))
    try:
        labels, Ks = Tcut_for_bipartite_graph(B, Ks, eskMethod, gapth,
                                              maxTcutKmIters,
                                              cntTcutKmReps,
                                              clusterMethod=clusterMethod)
    except:
        labels, Ks = Tcut_for_bipartite_graph(B, Ks[0], eskMethod, gapth,
                                              maxTcutKmIters, cntTcutKmReps,
                                              clusterMethod=clusterMethod)
    if mode == 'secuer':
        return labels
    else:
        return labels, Ks


if __name__ == '__main__':
    file = 'gold_label_data/'
    files = os.listdir(file)
    fileh5sd_gold = [i for i in files if i.endswith('.h5ad')]
    nmi_secuer = []
    ari_secuer = []
    t_secuer = []
    for i in range(len(fileh5sd_gold)):
        print(fileh5sd_gold[i])
        data = sc.read(file + fileh5sd_gold[i])
        fea = data.obsm['X_pca']
        start = time.time()
        res = secuer(fea, p=1000, Knn=7,
                     maxTcutKmIters=100,
                     cntTcutKmReps=3, seed=1,
                     multiProcessState=True,
                     num_multiProcesses=4
                     )
        end = time.time() - start
        print(np.unique(res).shape[0], np.unique(data.obs['celltype']).shape[0])
        nmi_secuer.append(normalized_mutual_info_score(data.obs['celltype'], res))
        ari_secuer.append(adjusted_rand_score(data.obs['celltype'], res))
        t_secuer.append(end)
        print()
    res_secuer = pd.DataFrame([t_secuer, nmi_secuer, ari_secuer],
                              index=['t_secuer', 'nmi_secuer', 'ari_secuer'],
                              columns=fileh5sd_gold)
    print(res_secuer.values)
