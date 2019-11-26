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
> rather than transcription factor regulation.

This is the analysis code for our paper at https://doi.org/10.1016/j.bbagrm.2019.194444

It was used to generate all figures in the article ([report](report) directory).
The main conclusions are outlined below.

## Findings

### DREAM evaluations differ from cancer data

![](https://drive.google.com/file/d/1gj95RzDUoGmC5DWTbK1WRCRz4nKCF1nr/preview)

Gene network inference methods have been extensively validated using consortium
benchmarks such as the [DREAM challenges](http://dreamchallenges.org/). However,
the cancer data these methods are often used on have many more genes and fewer
samples than those in the evaluation.

### Cancers are aneuploid and impure, but this has not been modeled

![](https://drive.google.com/file/d/1N-0uVIXWVL7tGuJqxiuI0lOgWsuQvfUw/preview)

Previous evaluations have also not considered specific biological properties of
cancer data: they show chromosome copy number changes that influences gene
expression in a coordinated fashion, and they are mixtures of different cell
types. These confounding factors are pervasive and not limited to individual
samples.

### Copy number alterations influence inferred gene networks

![](https://drive.google.com/file/d/1o1gaQiCkCInlJSGMOkW223ol1nijod2y/preview)

Regulatory links are often inferred for genes in the same currently amplified
region. They do, however, only comprise few genes so the effect on the total
network is small. Aneuploidies (whole chromosome copy number changes) introduce
a smaller fraction of false positive links, but this has a major influence on
the inferred network as a whole.

### Sample mixtures influence inferred gene networks

![](https://drive.google.com/file/d/1B3gHawzP7nXH95E90Es_LtatYQ1W2StD/preview)

Sample mixtures also show a strong effect on false positive regulatory links.
Looking at only cancer vs. non-cancer (purity of the sample), we find that
FP links are also often inferred between genes that correlate with sample
purity.

## Running the code

If you want to run this code, you will require:

* All R packages loaded in the scripts
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) for workflow execution
* [ebits](https://github.com/mschubert/ebits) and [data](https://github.com/mschubert/data)
  repositories with set up [module](https://github.com/klmr/modules) paths

Given that you have set up these packages and data, you can run the analyses:

```sh
cd set_enrichment
snakemake
```

And generate the figures:

```sh
cd report
make
```
