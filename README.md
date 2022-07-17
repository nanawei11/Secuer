# Secuer: ultrafast, scalable and accurate clustering of single-cell RNA-seq data

`Secuer` is a superfast and scalable clustering algorithm for (ultra-)large scRNA-seq data analysis based on spectral clustering.  `Secuer-consensus` is a consensus clustering algorithm with Secuer as a subroutine. In addition, `Secuer` can also be applied to other large-scale omics data with two-dimenational (features * obsevation). For more details see xxx. Secuer is available in Python. 

The workflow of Secuer:

<img src="D:\My_data\Allproject\Secuer\USPEC\Secuer\Figures\Figure1.png" style="zoom: 33%;" />

## Installation
`Secuer` requires [python](https://www.python.org)  to run. 

```python
pip install Secuer
```

## Run Seucer (usage)

#### Essential paramters

To run Secuer with default parameters, you only need to give:

- -i INPUTFILE 

  two-dimensional data (features by observations) file for clustering

Example for run Secuer with all default parameters:

```sh
$ Secuer S -i ${inputpath} 
```


Example for run Secuer-consensus with all default parameters:
```sh
$ Secuer C -i ${inputpath}
```


#### options
You can also specify the following options:

- -p         

  The number of anchors, default by 1000.

* -o OUTFILE

  Output file directory and file name, default by output.

- --knn
  The number of k nearest neighbors anchors, default by 7.

- --distance

  The metrics measuring the dissimmlarity between cells or anchors, default by euclidean.

## Output files

1. `output/output.txt` is the clustering reslut:
     - The clustering label.
     

## Citation

