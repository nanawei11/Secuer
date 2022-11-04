from secuer.secuer import secuer
from secuer.secuerconsensus import secuerconsensus
import os
import scanpy as sc
file = 'D://My_data//Allproject//Secuer//Clustering0804//gold_label_data/'
files = os.listdir('D://My_data//Allproject//Secuer//Clustering0804//gold_label_data/')
files = os.listdir('D://My_data//Allproject//Secuer//Clustering0804//gold_label_data/')
print(files)
fileh5sd_gold = [i for i in files if i.endswith('.h5ad')]
i=0
data = sc.read(file + fileh5sd_gold[i])
print(data)
fea = data.obsm['X_pca']
res = secuer(fea= fea,
    Knn=5,
    multiProcessState=True,
    num_multiProcesses=4)
res = secuerconsensus(run_secuer=True,
                      fea= fea,
                     Knn=5,
                     M=5,
                    multiProcessState=True,
                    num_multiProcesses=4)

