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



Secuer S -i ./example_data/Biase_k3_FPKM_scRNA --yaml ./config.yaml -o ./Biase_result -p 1000 --knn 5 --transpose


Hello, thank you for your feedback. This issue is due to an upgrade to the  function `scipy.spatial.distance.cdist` in Python. We have uploaded the newest version to [https://github.com/nanawei11/Secuer/](https://github.com/nanawei11/Secuer/) .You can also find it at [https://pypi.org/project/secuer/1.0.9/](https://pypi.org/project/secuer/1.0.9/).

Here is an example:
```
from secuer.secuer import secuer
from secuer.secuerconsensus import secuerconsensus
import os
import scanpy as sc
data = sc.read('Biase_pca.h5ad')
print(data)
res = secuer(fea= fea,
               p=1000)
res = secuerconsensus(fea= fea,
                      Knn=5,
                      p=1000
                      M=5)

```

```AnnData object with n_obs × n_vars = 49 × 3922
    obs: 'celltype', 'n_genes_by_counts', 'total_counts', 'total_counts_ercc', 'pct_counts_ercc'
    var: 'dropouts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'mean', 'std', 'ercc', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts'
    uns: 'hvg', 'log1p', 'pca'
    obsm: 'X_pca'
    varm: 'PCs'

[2022-11-04 19:05:15] [INFO] Selecting representatives...
[2022-11-04 19:05:15] [INFO] Approximate KNN...
[2022-11-04 19:05:15] [INFO] Estimating the number of clustering...
[2022-11-04 19:05:15] [INFO] Bipartite graph partitioning...

[2022-11-04 19:05:51] [INFO] Running secuer 1
[2022-11-04 19:05:51] [INFO] Running secuer 2
[2022-11-04 19:05:51] [INFO] Running secuer 4
[2022-11-04 19:05:51] [INFO] Running secuer 3
[2022-11-04 19:05:51] [INFO] Running secuer 5
```