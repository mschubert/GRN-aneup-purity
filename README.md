# Gene networks in cancer are biased by aneuploidies and sample impurities

> Gene regulatory network inference is a standard technique for obtaining
> structured regulatory information from, for instance, gene expression
> measurements. Methods performing this task have been extensively evaluated on
> synthetic, and to a lesser extent real data sets. In contrast to these test
> evaluations, applications to gene expression data of human cancers are often
> limited by fewer samples and more potential regulatory links, and are biased
> by copy number aberrations as well as cell mixtures and sample impurities.
> Here, we take networks inferred from TCGA cohorts as an example to show that
> (1) transcription factor annotations are essential to obtain reliable
> networks, and (2) even for state of the art methods, we expect that between
> 20 and 80% of edges are caused by copy number changes and cell mixtures
> rather than transcription factor regulation.]

This is the analysis code for our paper at https://doi.org/10.1016/j.bbagrm.2019.194444

It was used to generate all figures in the article ([report](report) directory).

If you want to run this code, you will require:

* All R packages loaded in the scripts
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) for workflow execution
* [ebits](https://github.com/mschubert/ebits) and [data](https://github.com/mschubert/data)
  repositories with set up [module](https://github.com/klmr/modules) paths
