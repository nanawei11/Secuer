# Secuer: ultrafast, scalable and accurate clustering of single-cell RNA-seq data

`Secuer` is a superfast and scalable clustering algorithm for (ultra-)large scRNA-seq data analysis based on spectral clustering.  `Secuer-consensus` is a consensus clustering algorithm with Secuer as a subroutine. In addition, `Secuer` can also be applied to other large-scale omics data with two-dimensional (features by observations). For more details see [secuer](https://doi.org/10.1371/journal.pcbi.1010753).  

The workflow of Secuer:

<img src="https://github.com/nanawei11/Secuer/raw/main/Figures/Figure1.png" style="zoom: 33%;" />

## Installation

`Secuer` is available in [python](https://www.python.org). 

```python
# use anaconda
conda create -n secuer python=3.9
conda activate secuer 
pip install secuer matplotlib pandas scanpy igraph louvain pyyaml

# or 
pip install secuer
```

## Run Seucer (usage)

#### Essential parameters

To run Secuer with default parameters, you only need to specify:

- -i 

  scRNA-seq data (cells by genes) file for clustering. 

* --yaml 

  The parameters of data preprocessing. see [config.yaml](https://github.com/nanawei11/Secuer/blob/main/config.yaml) for more details.

#### options
You can also specify the following options:

- -p         

  The number of anchors, default by 1000.

- -o 

  Output file directory and file name, default by output.

- --knn 

  The number of k nearest neighbors anchors, default by 7.

- --distance

  The metrics measuring the dissimilarity between cells or anchors, default by euclidean.
  
- --transpose
  
  Require it if your data is a .csv, .txt or tsv file with features by observations.
  
- --eskMethod

  Specify the method used for estimated the number of cluster, default by subGraph.

* --eskResolution

  Specify the resolution when `--eskMethod`  is subGraph, default by 0.8.

* --gapth

  Specify the gapth largest value when `--eskMethod`  is not subGraph.


Example for run Secuer with custom parameters:

```sh
$ Secuer S -i ./example_data/Biase_k3_FPKM_scRNA.csv --yaml ./config.yaml -o ./Biase_result -p 1000 --knn 5 --transpose
```

### Output files

1. `output/SecuerResult.txt` is the clustering result. 
2. `output/SecuerResult.h5ad` is the preprocessed data with the clustering result.

## Run Seucer-consensus (usage)

#### Essential parameters

To run Secuer-consensus with default parameters, you only need to specify:

- -i 

  two-dimensional data (observations by features) file for clustering. 

* --yaml 

  The parameters of data preprocessing. see [config.yaml](https://github.com/nanawei11/Secuer/blob/main/config.yaml) for more details.

#### options
You can also specify the following options:

- -p         

  The number of anchors, default by 1000.

- -o 

  Output file directory and file name, default by outputCon.

* -M 

  The times to run secuer.

* --knn 

  The number of k nearest neighbors anchors, default by 7.

- --transpose
  Require it if your data is a .csv, .txt or tsv file with genes by cells, default by False.

Example for run Secuer-consensus:
```sh
$ Secuer C -i ./example_data/Biase_k3_FPKM_scRNA.csv --yaml ./config.yaml -o ./Biase_conresult  -p 900 --knn 5 -M 7 --transpose
```

### Output files

1. `output/SecuerConsensusResult.txt` is the clustering result. 
2. `output/SecuerConsensusResult.h5ad` is the preprocessed data with the clustering result.

## Or run Secuer in Python

```python
import scanpy as sc
import secuer as sr
data = sc.read('example_data/Biase_k3_FPKM_scRNA.csv').T
# data preprocessing
sc.pp.filter_genes(data, min_counts=1)
sc.pp.filter_cells(data, min_counts=1)
sc.pp.normalize_total(data, target_sum=1e4)
sc.pp.log1p(data)
sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5)
data = data[:, data.var.highly_variable]
sc.pp.scale(data, max_value=10)
sc.tl.pca(data)

# run secuer
fea = data.obsm['X_pca']
res = sr.secuer(fea= fea,
                Knn=5,
                multiProcessState=True,
                num_multiProcesses=4)

# run secuer-consensus
resC = sr.secuerconsensus(run_secuer=True,
                          fea= fea,
                          Knn=5,
                          M=5,
                          multiProcessState=True,
                          num_multiProcesses=4)
```



## Citation

Wei N, Nie Y, Liu L, Zheng X, Wu H-J (2022) Secuer: Ultrafast, scalable and accurate clustering of single-cell RNA-seq data. PLOS Computational Biology 18(12): e1010753. https://doi.org/10.1371/journal.pcbi.1010753.

