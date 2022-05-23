import gc
import sys
import os
import psutil
from BigC import *
# random_rep=True,adj_param=False can remove
def BigCC(run_BigC,
          fea,
          k=None,
          M=20,
          p=1000,
          Knn=5,
          eskMethod='subGraph',
          gapth=4,
          eskReso=0.8,
          addweights=False,
          random_rep=True,
          adj_param=False):

    baseCls,ks = BigC_EnsembleGeneration(run_BigC = run_BigC,
                                       fea = fea, M=M, p= p,
                                       Knn=Knn,random_rep=random_rep,
                                       adj_param=adj_param,
                                       eskMethod=eskMethod,
                                       gapth=gapth,
                                       eskReso=eskReso,
                                       addweights=addweights
                                       )
    print('Performing the consensus function...')
    if not k:
        k=Counter(ks).most_common(1)[0][0]
    # tic1 = time.time()
    L = BigC_ConsensusFunction(baseCls, k)
    # print(time.time() - tic1)
    return L,k

# random_rep=True,adj_param=False can remove
def BigC_EnsembleGeneration(run_BigC,
                             fea, M, p=1000,
                             eskMethod='subGraph',
                             gapth=4,
                             eskReso=0.8,
                             addweights=False,
                             Knn=5,random_rep=True,
                             adj_param=False):
    '''
    Generate M base cluserings.
    The number of clusters in each base clustering is randomly selected in
    the range of [lowK, upK].
    '''
    N = fea.shape[0]
    if p > N:
        p = N
    # rand('state',sum(100*clock)*rand(1)) % Reset the clock before generating random numbers
    # Ks = np.random.choice(upK - lowK + 1, M) + lowK - 1
    # In ensemble generation, the iteration number in the kmeans discretization of each base cluserer can be        set to small values, so as to improve
    # diversity of base clusterings and reduce the iteration time costs.
    tcutKmIters = 5
    tcutKmRps = 1
    if N<1000:
        resolution = [float(i / 10) for i in range(1, M + 1)]
    else:
        if M < 11:
            resolution = [float(i / 10) for i in range(2, M*2+1, 2)]
        elif M<21 and M>=11:
            resolution = [float(i / 10) for i in range(2, M+2)]
        else:
            np.random.seed(1)
            resolution = np.random.choice([float(i / 10) for i in range(1, M+1)],M)
    # resolution = [float(i / 10) for i in range(1, M + 1)]
    np.random.seed(1) # set random seet
    distance1 = np.random.choice(['euclidean', 'cosine'], M)
    if run_BigC:
        print(f'random select representatives')
        members = []
        k = []
        for j in range(M):
            # print(j)
            res,ks = BigC(fea=fea,eskMethod=eskMethod,gapth=gapth,
                        eskReso=resolution[j],mode='Consensus',
                        distance=distance1[j],p=p,Knn=Knn,
                        maxTcutKmIters=tcutKmIters,cntTcutKmReps=tcutKmRps,
                        seed=j,addweights=addweights,
                        adj_param=adj_param)
            members += [res.tolist()]
            k += [ks]
        members = np.array(members).T
    else:
        members = fea

    return members,k  # N by M cluster

def BigC_ConsensusFunction(baseCls, k, gapth=4,eskMethod = 'subGraph',
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
    labels,ks = Tcut_for_bipartite_graph(B, k,gapth=gapth,eskMethod=eskMethod,
                                      maxKmIters=maxTcutKmIters,
                                      cntReps = cntTcutKmReps)
    return labels

if __name__ == '__main__':
    file = 'D://My_data//Allproject//单细胞聚类//Clustering0804//gold_label_data/'
    files = os.listdir('D://My_data//Allproject//单细胞聚类//Clustering0804//gold_label_data/')
    fileh5sd_gold = [i for i in files if i.endswith('.h5ad')]
    fileh5sd_gold
    nmi_BigCC = []
    ari_BigCC = []
    t_BigCC = []
    for i in range(len(fileh5sd_gold)):
        print(fileh5sd_gold[i])
        data = sc.read(file + fileh5sd_gold[i])
        fea = data.obsm['X_pca']
        start = time.time()
        try:
            sc.pp.neighbors(data,n_pcs=48)
        except:
            sc.pp.neighbors(data, n_pcs=50)
        res,k = BigCC(run_BigC=True,fea= fea,Knn=5,M=6)
        print(np.unique(res).shape[0])
        end = time.time() - start
        print(np.unique(res).shape[0],np.unique(data.obs['celltype']).shape[0])
        name_up = 'BigCC_lable_' + str(i)
        data.obs[name_up] = pd.Categorical(res)
        nmi_BigCC.append(normalized_mutual_info_score(data.obs['celltype'], res))
        ari_BigCC.append(adjusted_rand_score(data.obs['celltype'], res))
        t_BigCC.append(end)
        print()
    res_BigCC = pd.DataFrame([t_BigCC, nmi_BigCC, ari_BigCC],
                             index=['t_BigCC', 'nmi_BigCC', 'ari_BigCC'],
                             columns=fileh5sd_gold)
    print(res_BigCC.values)
