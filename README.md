[//]: # (TODO: Clean and update code snippets)
[//]: # (TODO: Update GetData docs)
[//]: # (TODO: Turn off syntax highlighting in code snippets (that is incorrect)

# MethylMix

[![Build Status](https://www.bioconductor.org/shields/build/release/bioc/MethylMix.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/MethylMix/)
[![Downloads](http://bioconductor.org/shields/downloads/MethylMix.svg)](http://bioconductor.org/packages/stats/bioc/MethylMix/)
[![License](https://img.shields.io/badge/license-GPL-yellow.svg)](https://opensource.org/licenses/GPL-2.0)
[![Version](https://img.shields.io/badge/version-2.10.0-lightgrey.svg)](https://www.bioconductor.org/packages/release/bioc/html/MethylMix.html)

An R package for identifying DNA methylation driven genes.

## Table of Contents

- [Getting Started](#getting-started)
- [Introduction](#introduction)
- [Data Access and Preprocessing](#data-access-and-preprocessing)
- [Data Input for MethylMix](#data-input-for-methylmix)
- [Running MethylMix](#running-methylmix)
- [Plotting Outputs](#plotting-outputs)
- [References](#references)
- [Author Information](#author-information)
- [Citation](#citation)
- [License](#license)


## Getting Started

__Installation__

To install the _MethylMix_ package, the easiest way is through Bioconductor:

```{r, eval=FALSE}
source("http://bioconductor.org/biocLite.R")
biocLite(MethylMix)
```

Alternatively, MethylMix can be installed by downloading the appropriate file for 
your platform from the [Bioconductor website](https://www.bioconductor.org/packages/release/bioc/html/MethylMix.html). 
For Windows, start R and select the `Packages` menu, then `Install package from local zip file`.
Find and highlight the location of the zip file and click on `open`. For Linux/Unix, use the usual command 
`R CMD INSTALL` or install from bioconductor.

__Loading__

To load the `MethylMix` package in your R session, type `library(MethylMix)`.

__Help__ 

Detailed information on `MethylMix` package functions can be obtained in the help files. For example, to view the help file 
for the function `MethylMix` in a R session, use `?MethylMix`.

## Introduction

DNA methylation is a mechanism whereby a methyl-group is added onto a CpG site. Methylation of these CpG sites is associated with gene silencing and is an important mechanism for normal tissue development and is often involved in diseases such as cancer. Recently, many high throughput data has been generated profiling CpG site methylation on a genome wide basis. This has created large amounts of data on DNA methylation for many diseases. Computational analysis of DNA methylation data is required to identify potentially aberrant DNA methylation compared to normal tissue. MethylMix ([Gevaert 2015](https://academic.oup.com/bioinformatics/article/31/11/1839/2364827); [Gevaert, Tibshirani & Plevritis 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0579-8)) was developed to tackle this question using a computational approach. MethylMix identifies differential and functional DNA methylation by using a beta mixture model to identify subpopulations of samples with different DNA methylation compared to normal tissue. Functional DNA methylation refers to a significant negative correlation based on matched gene expression data. MethylMix outputs hyper and hypomethylated genes, called MethylMix drivers, which can be used for downstream analysis. MethylMix was designed for cis-regulated promoter differential methylation and works best when specific CpG sites profiled are associated with a gene. For example, when using data from the 27k and 450k Illumina Infinium platforms. 

## Data Access and Preprocessing

The `MethylMix` package provides functions to access and preprocess data from [The Cancer Genome Atlas (TCGA)](https://cancergenome.nih.gov) portal. Given one of TCGA codified cancer types and a path to save downloaded files, the downloading and preprocessing of the data can be executed using the `GetData()` function. For example, for ovarian cancer, the data can be downloaded and preprocessed as follows:

```{r, eval=FALSE}
cancerSite <- "OV"
targetDirectory <- paste0(getwd(), "/")
GetData(cancerSite, targetDirectory)
```

All functions in the `MethylMix` package can be run in parallel, if the user provides a parallel setup, like the following:

```{r, eval=FALSE}
cancerSite <- "OV"
targetDirectory <- paste0(getwd(), "/")

library(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)
GetData(cancerSite, targetDirectory)
stopCluster(cl)
```

The `GetData()` function downloads DNA methylation data and gene expression data. The methylation data is provided using 27k or 450k Illumina platforms. If both 27k and 450k files are found, the data is carefully combined. For gene expression, either microarray (Agilent), RNA sequencing data or both are available. The `MethylMix` package downloads RNA sequencing data for all cancer sites except for ovarian and gliobastoma cancer sites (due to few RNA-seq samples available). For the preprocessing of the data, we perform missing value estimation and batch correction (using [Combat](https://bioconductor.org/packages/release/bioc/html/sva.html)). Finally, as in the TCGA case, when only probe level Illumina data is available, mapping probes to genes is recommended before building mixture models. This allows the algorithm to prioritize cis-regulated differential methylation by only focusing on differential methylation of CpG sites to their closest gene transcripts. Both the 27k and 450k Illumina Infinium platforms have database R packages that provide the necessary mapping information. We use the annotation to map probes to genes, before clustering the probes within each gene. This whole process can take a long time, depending of the size of the data.

It is also possible to perform each one of these tasks individually using other functions in the `MethylMix` package, as in the following example:


```{r, eval=FALSE, tidy=TRUE}
cancerSite <- "OV"
targetDirectory <- paste0(getwd(), "/")

cl <- makeCluster(5)
registerDoParallel(cl)

# Downloading methylation data
METdirectories <- Download_DNAmethylation(cancerSite, targetDirectory)
# Processing methylation data
METProcessedData <- Preprocess_DNAmethylation(cancerSite, METdirectories)
# Saving methylation processed data
saveRDS(METProcessedData, file = paste0(targetDirectory, "MET_", cancerSite, "_Processed.rds"))

# Downloading gene expression data
GEdirectories <- Download_GeneExpression(cancerSite, targetDirectory)
# Processing gene expression data
GEProcessedData <- Preprocess_GeneExpression(cancerSite, GEdirectories)
# Saving gene expression processed data
saveRDS(GEProcessedData, file = paste0(targetDirectory, "GE_", cancerSite, "_Processed.rds"))

# Clustering probes to genes methylation data
METProcessedData <- readRDS(paste0(targetDirectory, "MET_", cancerSite, "_Processed.rds"))
res <- ClusterProbes(METProcessedData[[1]], METProcessedData[[2]])

# Putting everything together in one file
toSave <- list(METcancer = res[[1]], METnormal = res[[2]], GEcancer = GEProcessedData[[1]], GEnormal = GEProcessedData[[2]], ProbeMapping = res$ProbeMapping)
saveRDS(toSave, file = paste0(targetDirectory, "data_", cancerSite, ".rds"))

stopCluster(cl)
```

For user supplied data not from TCGA, the user can provide DNA methylation beta-values (normal and cancer) and gene expression data in the form of data.matrix objects in R with the rows corresponding to the genes and the columns corresponding to the samples. Once the three matrices: METcancer, METnormal, and GEnormal are created, the last step is to cluster the methylation data:

```{r, eval=FALSE, tidy=TRUE}
METcancer = matrix(data = methylation_data, nrow = nb_of_genes, ncol = nb_of_samples)
METnormal = matrix(data = methylation_data, nrow = nb_of_genes, ncol = nb_of_samples)
GEcancer = matrix(data = expression_data, nrow = nb_of_genes, ncol = nb_of_samples)
ClusterProbes(MET_Cancer, MET_Normal, CorThreshold = 0.4)
```

If the data contains batches, the user must provide numeric batch data within the matrices. MethylMix can be applied on all Illumina arrays, including the Epic platform, and any array that outputs beta values. At the moment there are no restrictions for input of sequencing-based methylation data, as long as the data is formatted in the proper proportions. However, as mixture modeling is computationally expensive, MethylMix may require more time to finish in these cases.

## Data Input for MethylMix

To run MethylMix, three data sets of a particular disease are required. The first one is the methylation data for the disease samples, `METcancer`, which allows the identification of methylation states associated with a disease for each gene of interest. The second is a set of appropriate normal or baseline methylation data, `METnormal`, which is used to distinguish between hyper or increased methylation versus hypo or decreased methylation. Finally, the third data set is matched gene expression data for the disease samples, `GEcancer`, which is used to identify functional differential methylation by focusing only on differential methylation that has a significant inversely correlated effect on gene expression. 

Each of these three data sets are matrix objects with genes in the 
rows with unique rownames (e.g. gene symbols) and samples or patients in the 
columns with unique patient names. The `GetData` function described before
saves an R object which contains these matrices in the correct format.

For example, small data sets from TCGA of glioblastoma samples are avaliable in the `MethylMix` package:

```{r, tidy=TRUE}
library(MethylMix)
library(doParallel)
data(METcancer)
data(METnormal)
data(GEcancer)
head(METcancer[, 1:4])
head(METnormal)
head(GEcancer[, 1:4])
```

## Running MethylMix

Using the example glioblastoma data provided in the package, `MethylMix` can be used to identify hypomethylated and hypermethylated genes as follows:

```{r, tidy=TRUE, warning=F}
MethylMixResults <- MethylMix(METcancer, GEcancer, METnormal)
```

Or in parallel:

```{r, tidy=TRUE, eval=FALSE}
library(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)
MethylMixResults <- MethylMix(METcancer, GEcancer, METnormal)
stopCluster(cl)
```

The output from the `MethylMix` function is a list with the following elements:

* `MethylationDrivers`: Genes identified as transcriptionally predictive and differentially methylated by MethylMix.
* `NrComponents`: An integer of the number of methylation states found for each driver gene.
* `MixtureStates`: A list with the DM-values for each driver gene.
* `MethylationStates`: A matrix with DM-values for all driver genes (rows) and all samples (columns).
* `Classifications`: A matrix with integers indicating to which mixture component each cancer sample was assigned to, for each gene.
* `Models`: Beta mixture model parameters for each driver gene.

Differential Methylation values (DM-values) are defined as the difference between 
the methylation mean in one mixture component of cancer samples and the methylation mean in the normal samples, for a given gene. These outputs can be viewed in R as follows:

```{r, tidy=TRUE}
MethylMixResults$MethylationDrivers
MethylMixResults$NrComponents
MethylMixResults$MixtureStates
MethylMixResults$MethylationStates[, 1:5]
MethylMixResults$Classifications[, 1:5]
MethylMixResults$Models
```

## Plotting Outputs

The `MethylMix` package also provides functions to visually represent the findings:

```{r, tidy=TRUE, eval=F}
# Plot the most famous methylated gene for glioblastoma
plots <- MethylMix_PlotModel("MGMT", MethylMixResults, METcancer)
plots$MixtureModelPlot
```

<img src="vignettes/plot1.png" width="500">

```{r, tidy=TRUE, eval=F}
# Plot MGMT also with its normal methylation variation
plots <- MethylMix_PlotModel("MGMT", MethylMixResults, METcancer, METnormal = METnormal)
plots$MixtureModelPlot
```

<img src="vignettes/plot2.png" width="500">

```{r, tidy=TRUE, eval=F}
# Plot a MethylMix model for another gene
plots <- MethylMix_PlotModel("ZNF217", MethylMixResults, METcancer, METnormal = METnormal)
plots$MixtureModelPlot
```

<img src="vignettes/plot3.png" width="500">

```{r, tidy=TRUE, eval=F}
# Also plot the inverse correlation with gene expression (creates two separate plots)
plots <- MethylMix_PlotModel("MGMT", MethylMixResults, METcancer, GEcancer, METnormal)
plots$MixtureModelPlot
plots$CorrelationPlot
```

<img src="vignettes/plot2.png" width="300"> <img src="vignettes/plot4.png" width="300">

```{r, eval = FALSE, tidy=TRUE}
# Plot all functional and differential genes
for (gene in MethylMixResults$MethylationDrivers) {
    MethylMix_PlotModel(gene, MethylMixResults, METcancer, METnormal = METnormal)
}
```

## References

>Gevaert O. MethylMix: an R package for identifying DNA methylation-driven genes. Bioinformatics (Oxford, England). 2015;31(11):1839-41. doi:10.1093/bioinformatics/btv020 \n

>Gevaert O, Tibshirani R, Plevritis SK. Pancancer analysis of DNA methylation-driven genes using MethylMix. Genome Biology. 2015;16(1):17. doi:10.1186/s13059-014-0579-8 ",


## Author Information

| Olivier Gevaert 
|------------------
| Stanford Center for Biomedical Informatics 
| Department of Medicine 
| 1265 Welch Road 
| Stanford CA, 94305-5479 

## Citation

If you use MethylMix in your work, please cite:

> Gevaert, O. (2017). MethylMix: Identifying methylation driven cancer genes. R package version 2.10.0.

## License

MethylMix is licensed under the GPL-2 license. See the [LICENSE](https://github.com/gevaertlab/MethylMix/LICENSE) for more information.
