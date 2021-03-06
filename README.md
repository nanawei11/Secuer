# Secuer: ultrafast, scalable and accurate clustering of single-cell RNA-seq data

`Secuer` is a superfast and scalable clustering algorithm for (ultra-)large scRNA-seq data analysis based on spectral clustering.  `Secuer-consensus` is a consensus clustering algorithm with Secuer as a subroutine. In addition, `Secuer` can also be applied to other large-scale omics data with two-dimenational (features * obsevation). For more details see [secuer](https://arxiv.org/abs/2205.12432v2).  

The workflow of Secuer:

<img src="https://github.com/nanawei11/Secuer/raw/main/Figures/Figure1.png" style="zoom: 33%;" />

## Installation

`Secuer` is available in [python](https://www.python.org). 

```python
pip install secuer
```

## Run Seucer (usage)

#### Essential parameters

To run Secuer with default parameters, you only need to specify:

- -i INPUTFILE 

  scRNA-seq data (cells by genes) file for clustering. 

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
$ Secuer S -i ./example_data/deng-rpkms.csv -o ./deng_result -p 1000 --knn 5 --transpose
```

## Run Seucer-consensus (usage)

#### Essential parameters

To run Secuer-consensus with default parameters, you only need to specify:

- -i INPUTFILE 

  two-dimensional data (observations by features) file for clustering. 

#### options
You can also specify the following options:

- -p         

  The number of anchors, default by 1000.

- -o 

  Output file directory and file name, default by outputCon.

* -M 

???	The times to run secuer.

* --knn 

???		The number of k nearest neighbors anchors, default by 7.

- --transpose
  Require it if your data is a .csv, .txt or tsv file with genes by cells, default by False.

Example for run Secuer-consensus:
```sh
$ Secuer C -i ./example_data/deng-rpkms.csv -o ./deng_conresult  -p 900 --knn 5 -M 7 --transpose
```



## Output files

1. `output/output.txt` is the clustering. 

## Citation
