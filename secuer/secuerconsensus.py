import gc
import os
import numpy as np
from collections import Counter
import scanpy as sc
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
import time
import pandas as pd
import sys
from scipy.sparse import issparse, csr_matrix
from secuer.secuer import (secuer, Tcut_for_bipartite_graph)
# random_rep=True,adj_param=False can remove
# __all__ = ['secuerconsensus']
def secuerconsensus(fea,
          k=None,
          run_secuer=True,
          M=20,
          p=1000,
          Knn=5):

    baseCls,ks = secuerC_EnsembleGeneration(run_secuer = run_secuer,
                                       fea = fea, M=M, p= p,
                                       Knn=Knn
                                       )
    print('Performing the consensus function...')
    if not k:
        k=Counter(ks).most_common(1)[0][0]
    L = secuerC_ConsensusFunction(baseCls, k)
    return L,k

# random_rep=True,adj_param=False can remove
def secuerC_EnsembleGeneration(run_secuer,
                             fea, M, p=1000,
                             Knn=5):
    '''
    Generate M base cluserings.
    The number of clusters in each base clustering is randomly selected in
    the range of [lowK, upK].
    '''
    N = fea.shape[0]
    if p > N:
        p = N

    tcutKmIters = 5
    tcutKmRps = 1
    if N<1000:
        resolution = [float(i / 10) for i in range(1, M + 1)]
    else:
        if M < 11:
            resolution = [float(i / 10) for i in range(2, M*2+1, 2)]
        elif  M>=11 and M<21:
            resolution = [float(i / 10) for i in range(2, M+2)]
        else:
            np.random.seed(1)
            resolution = np.random.choice([float(i / 10) for i in range(1, M+1)],M)
    # resolution = [float(i / 10) for i in range(1, M + 1)]
    np.random.seed(1) # set random seet
    distance1 = np.random.choice(['euclidean', 'cosine'], M)
    if run_secuer:
        members = []
        k = []
        for j in range(M):
            print(f'Running secuer {j+1}')
            res,ks = secuer(fea=fea,eskMethod='subGraph',
                        eskResolution
                            =resolution[j],mode='Consensus',
                        distance=distance1[j],p=p,Knn=Knn,
                        maxTcutKmIters=tcutKmIters,cntTcutKmReps=tcutKmRps,
                        seed=1)
            members += [res.tolist()]
            k += [ks]
        members = np.array(members).T
    else:
        members = fea

    return members,k  # N by M cluster

def secuerC_ConsensusFunction(baseCls, k,
                            maxTcutKmIters=100, cntTcutKmReps=3):
    # Combine the M base clusterings in baseCls to obtain the final clustering
    # result (with k clusters).
    N, M = baseCls.shape
    maxCls = np.max(baseCls, axis=0) + 1
    # print(f'k:{k},maxCls:{maxCls}')
    maxCls = np.cumsum(maxCls)
    baseCls = baseCls + np.concatenate([[0], maxCls[0:-1]])
    cntCls = maxCls[-1]
    # Build the bipartite graph.
    indptr = np.arange(0, N * M + 1, M)
    B = csr_matrix(([1] * N * M,  # copy the data, otherwise strange behavior here
                     baseCls.copy().ravel(), indptr.copy().ravel()), shape=(N, cntCls))
    # print(B.shape)
    del baseCls
    gc.collect()
    colB = np.sum(B, axis=0)
    B = B[:, (np.array(colB).flatten()) != 0]
    # Cut the bipartite graph.
    labels,ks = Tcut_for_bipartite_graph(B, k,eskMethod='subGraph',
                                      maxKmIters=maxTcutKmIters,
                                      cntReps = cntTcutKmReps)
    return labels

if __name__ == '__main__':
    file = 'D://My_data//Allproject//Secuer//Clustering0804//gold_label_data/'
    files = os.listdir('D://My_data//Allproject//Secuer//Clustering0804//gold_label_data/')
    fileh5sd_gold = [i for i in files if i.endswith('.h5ad')]
    fileh5sd_gold
    nmi_secuerC = []
    ari_secuerC = []
    t_secuerC = []
    for i in range(len(fileh5sd_gold)):
        print(fileh5sd_gold[i])
        data = sc.read(file + fileh5sd_gold[i])
        fea = data.obsm['X_pca']
        start = time.time()
        try:
            sc.pp.neighbors(data,n_pcs=48)
        except:
            sc.pp.neighbors(data, n_pcs=50)
        res,k = secuerconsensus(run_secuer=True,fea= fea,Knn=5,M=6)
        print(np.unique(res).shape[0])
        end = time.time() - start
        print(np.unique(res).shape[0],np.unique(data.obs['celltype']).shape[0])
        name_up = 'secuerC_lable_' + str(i)
        data.obs[name_up] = pd.Categorical(res)
        nmi_secuerC.append(normalized_mutual_info_score(data.obs['celltype'], res))
        ari_secuerC.append(adjusted_rand_score(data.obs['celltype'], res))
        t_secuerC.append(end)
        print()
    res_secuerC = pd.DataFrame([t_secuerC, nmi_secuerC, ari_secuerC],
                             index=['t_secuerC', 'nmi_secuerC', 'ari_secuerC'],
                             columns=fileh5sd_gold)
    print(res_secuerC.values)
