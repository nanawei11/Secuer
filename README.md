# Secuer: ultrafast, scalable and accurate clustering of single-cell RNA-seq data

`Secuer` is a superfast and scalable clustering algorithm for (ultra-)large scRNA-seq data analysis based on spectral clustering.  `Secuer-consensus` is a consensus clustering algorithm with Secuer as a subroutine. In addition, `Secuer` can also be applied to other large-scale omics data with two-dimenational (observations by features). For more details see xxx. Secuer is available in Python. 

The workflow of Secuer:

<img src="https://github.com/nanawei11/Secuer/raw/main/Figures/Figure1.png" style="zoom: 33%;" />

## Installation
`Secuer` requires [python](https://www.python.org)  to run. 

```python
pip install Secuer
```

## Run Seucer (usage)

#### Essential paramters

To run Secuer with default parameters, you only need to give:

- -i INPUTFILE 

  two-dimensional data (observations by features) file for clustering. 
 
#### options
You can also specify the following options:

- -p         

  The number of anchors, default by 1000.

- -o 

  Output file directory and file name, default by output.

- --knn 

  The number of k nearest neighbors anchors, default by 7.

- --distance

  The metrics measuring the dissimmlarity between cells or anchors, default by euclidean.
  
- --transpose
  Require it if your data is a .csv, .txt or tsv file with features by observations.

- --eskMethod

Example for run Secuer with all default parameters:

```sh
$ Secuer S -i ${inputpath} 
```
Example for run Secuer with custom parameters:

## Run Seucer-consensus (usage)

#### Essential paramters

To run Secuer-consensus with default parameters, you only need to give:

- -i INPUTFILE 

  two-dimensional data (observations by features) file for clustering. 
 
#### options
You can also specify the following options:

- -p         

  The number of anchors, default by 1000.

- -o 

  Output file directory and file name, default by outputCon.

- --knn 

  The number of k nearest neighbors anchors, default by 7.
  
- --transpose
  Require it if your data is a .csv, .txt or tsv file with features by observations.

```sh
$ Secuer S -i ./example_data/deng-rpkms.csv -o ./result -p 900 --knn 5 --transpose
```

Example for run Secuer-consensus with all default parameters:
```sh
$ Secuer C -i ${inputpath}
```





## Output files

1. `output/output.txt` is the clustering reslut from seucer.
2. `output/umap.pdf`   
3. `output/tsne.pdf`  

## Citation

