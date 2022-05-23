import argparse
import sys
import os
sys.path.append(os.path.abspath(r'/public/home/zhengxq/weinn/Clustering_scRNA/scRNA-seq/USPECcode'))
from USPECNEW0430 import *

def USPEC_ESK(fea,distance='euclidean', p=1000, Knn=5,seed=1,resolution = 0.8):
    N = fea.shape[0]  # n*F
    if p > N:
        p = N
    RpFea = getRepresentativesByHybridSelection(fea, p, seed=seed)
    # neighSize = 10 * Knn
    neighSize = 15
    RpfeaDist = pdist2_fast(RpFea, RpFea, distance)
    RpFeaKnnIdx = get_indices_distance_from_dense_matrix(RpfeaDist, neighSize+1,returnDis=False)  # p * neighSize+1
    ks = EstimatekbysubGraph(RpfeaDist=RpfeaDist,RpFeaKnnIdx=RpFeaKnnIdx,addweights=False,resolution=resolution,Knn=Knn)
    return ks

def TestAllMethod(data,label,file,method,sigma='cons',k=None,p=1000,knn=5,seed=1,
                  npc=50,adj_param = False):
    if method=='leiden':
        start=time.time()
        sc.pp.neighbors(data,n_pcs=npc)
        sc.tl.leiden(data)
        end=time.time()-start
        try:
            nmi= normalized_mutual_info_score(data.obs[label],data.obs['leiden'])
            ari= adjusted_rand_score(data.obs[label],data.obs['leiden'])
            label_leiden = data.obs['leiden'].astype(int).tolist()
        except:
            nmi=-1
            ari=-1
            label_leiden = data.obs['leiden'].astype(int).tolist()
        k = np.unique(label_leiden).shape[0]
        with open(file,'a') as f:
            f.write(f"{method}\t{knn}\t{k}\t{nmi}\t{ari}\t{end}\t{label_leiden}\n")
    elif method=='louvain':
        start=time.time()
        sc.pp.neighbors(data,n_pcs=npc)
        sc.tl.louvain(data)
        end=time.time()-start
        try:
            nmi = normalized_mutual_info_score(data.obs[label],data.obs['louvain'])
            ari = adjusted_rand_score(data.obs[label],data.obs['louvain'])
            label_louvain = data.obs['louvain'].astype(int).tolist()
        except:
            nmi=-1
            ari=-1
            label_louvain = data.obs['louvain'].astype(int).tolist()
        k = np.unique(label_louvain).shape[0]
        with open(file,'a') as f:
            f.write(f"{method}\t{knn}\t{k}\t{nmi}\t{ari}\t{end}\t{label_louvain}\n")
    elif method=='kmeans':
        try:
            fea = data.obsm['X_pca'][:,0:npc]
        except:
            fea = data.X
        if not k:
            k = USPEC_ESK(fea, distance='euclidean', p=p, Knn=knn, seed=seed)
        start = time.time()
        res = KMeans(n_clusters=k).fit(fea)
        end= time.time()-start
        try:
            nmi = normalized_mutual_info_score(data.obs[label],res.labels_)
            ari = adjusted_rand_score(data.obs[label],res.labels_)
            label_kmeans = res.labels_.tolist()
        except:
            ari=-1
            nmi=-1
            label_kmeans = res.labels_.tolist()
        with open(file,'a') as f:
            f.write(f"{method}\t{knn}\t{k}\t{nmi}\t{ari}\t{end}\t{label_kmeans}\n")
    elif method=='uspec':
        try:
            fea = data.obsm['X_pca'][:,0:npc]
        except:
            fea = data.X
        try:
            gt=data.obs[label].tolist()
        except:
            gt=1
        start = time.time()
        res = USPEC(fea,Ks=k,sigma=sigma,distance='euclidean',p=p,Knn=knn,seed=seed,adj_param=adj_param)
        end= time.time()-start
        try:
            nmi = normalized_mutual_info_score(gt,res)
            ari=adjusted_rand_score(gt,res)
            label_uspec = res.tolist()
        except:
            nmi=-1
            ari=-1
            label_uspec = res.tolist()
        k=(np.unique(label_uspec)).shape[0]
        with open(file,'a') as f:
            f.write(f"{method}\t{knn}\t{k}\t{nmi}\t{ari}\t{end}\t{label_uspec}\n")
def main():
    description='Clustering scRNA-seq using different methods.'
    parser = argparse.ArgumentParser(prog='Clustring',description=description, usage='%(prog)s [options]')
    # parser.add_argument('-h','--help',help="Show this help message and exit.")
    parser.add_argument('-i','--inputfile',dest='inputh5ad',
                        help='Input a .h5ad file.')
    parser.add_argument('-l','--label',dest='label',default='0',
                        help='Specify the label in input file.')
    parser.add_argument('-o','--outfile',dest='outfile',help='Output file.')
    parser.add_argument('-m','--method',dest='method',help='Choose a method from [louvain,leiden,kmeans,uspec].')
    parser.add_argument('-p', action='extend',dest='params', nargs='+',
                       help='Get the neighbors knn and the representatives p.')
    parser.add_argument('-k', dest='k', type=int,
                        help='Get the number of clusters k.',default=None)
    parser.add_argument('-s', dest='seed', type=int,
                        help='Seet the seed for representatives.',default=1)
    parser.add_argument('--sigma', dest='sigma',
                        help='Set the sigma using one of the [cons,noncons].', default='noncons')
    parser.add_argument('--adj_param',action='store_true' ,dest='adj_param',default=False,
                        help='Wheather esk using after adjusted params.[true or False]')
    # parser.add_argument('-D', action='extend',dest='Descan', nargs='+',
    #                    help='Get the DBSCAN params:[clusterMethod,eps,minsamples].')

    options = parser.parse_args()
    if not options.inputh5ad:
        parser.print_help()
        sys.exit(1)
    # if not options.k:
    #     try:
    #         k = np.unique(data.obs[label]).shape[0]
    #     except:
    #         print('Please input the number of clusters using -k.')
    # else:
    try:
        k = int(options.k)
    except:
        k = options.k
    sigma = options.sigma
    input_file = options.inputh5ad
    label = options.label
    parms = options.params
    method = options.method
    output_file = options.outfile
    seed=options.seed
    adj_param = options.adj_param
    # dparms = options.Descan

    # clusterMethod = str(dparms[0])
    # eps=float(dparms[1])
    # minsamples = int(dparms[2])
    try:
        p = int(parms[0])
        knn=int(parms[1])
        npc = int(parms[2])
    except:
        p=1000
        knn=5
        npc=50
    # filepath = '/public/home/zhengxq/weinn/Clustering_scRNA/scRNA-seq/Cluster0719/scDCC-master/data/Small_Datasets/results/'
    print(f"Your input parameters are: input:{input_file}, label:{label},p and knn:{parms},k:{k},method: {method},outfile:{output_file},seed:{seed},adjusted_para:{adj_param}")
    data=sc.read(input_file)
    savefile=output_file
    if method=='louvain' or method=='leiden':
        TestAllMethod(data=data,label=label,file = savefile,method=method,npc=npc)
    else:
        TestAllMethod(data=data, label=label,sigma=sigma, file=savefile, k=k, p=p, knn=knn, method=method,seed=seed,npc=npc,adj_param = adj_param)

if __name__=='__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)



