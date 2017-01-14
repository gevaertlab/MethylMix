#' BatchData data set
#'
#' Data set with batch number for TCGA samples.
#' @name BatchData
#' @docType data
#' @keywords datasets
#'
NULL

#' ProbeAnnotation data set
#'
#' Data set with annotation from Illumina methylatin arrays mapping CpG sites to genes.
#' @name ProbeAnnotation
#' @docType data
#' @keywords datasets
#'
NULL

#' SNPprobes data set
#'
#' Vector with probes with SNPs.
#' @name SNPprobes
#' @docType data
#' @keywords datasets
#'
NULL

#' @name GEcancer
#' @docType data
#' @title Cancer Gene expression data of glioblastoma patients from the TCGA project
#' @description
#' Cancer Gene expression data of glioblastoma patients from the TCGA project. 
#' A set of 14 genes that have been shown in the literature to be involved in 
#' differential methylation in glioblastoma were selected as an example to try out MethylMix.
#' @usage data(GEcancer)
#' @format A numeric matrix with 14 rows (genes) and 251 columns (samples).
#' @references 
#' Cancer Genome Atlas Research Network. Comprehensive genomic characterization
#' defines human glioblastoma genes and core pathways. Nature. 2008 Oct 23;
#' 455(7216):1061-8. doi: 10.1038/nature07385. Epub 2008 Sep 4. Erratum in: 
#' Nature. 2013 Feb 28;494(7438):506. PubMed PMID: 18772890; PubMed Central PMCID: PMC2671642.
#' @seealso TCGA: The Cancer Genome Atlas: \url{http://cancergenome.nih.gov/}
#' @keywords datasets
#'
NULL

#' @name METcancer
#' @docType data
#' @title DNA methylation data from cancer tissue from glioblastoma patients from the TCGA project
#' @description
#' Cancer Gene expression data of glioblastoma patients from the TCGA project.
#' A set of 14 genes that have been shown in the literature to be involved in 
#' differential methylation in glioblastoma were selected as an example to try 
#' out MethylMix.
#' @usage data(METcancer)
#' @format A numeric matrix with 14 rows (genes) and 251 columns (samples).
#' @references 
#' Cancer Genome Atlas Research Network. Comprehensive genomic characterization
#' defines human glioblastoma genes and core pathways. Nature. 2008 Oct 23;
#' 455(7216):1061-8. doi: 10.1038/nature07385. Epub 2008 Sep 4. Erratum in: 
#' Nature. 2013 Feb 28;494(7438):506. PubMed PMID: 18772890; PubMed Central PMCID: PMC2671642.
#' @seealso TCGA: The Cancer Genome Atlas: \url{http://cancergenome.nih.gov/}
#' @keywords datasets
#'
NULL

#' @name METnormal
#' @docType data
#' @title DNA methylation data from normal tissue from glioblastoma patients
#' @description
#' Normal tissue DNA methylation data of glioblastoma patients. These normal 
#' tissue samples were run on the same platform and are described in the      
#' publication referenced below.
#' @usage data(METnormal)
#' @format A numeric matrix with 14 rows (genes) and 4 columns (samples).
#' @references 
#' Noushmehr H, Weisenberger DJ, Diefes K, Phillips HS, Pujara K, Berman BP, 
#' Pan F, Pelloski CE, Sulman EP, Bhat KP, Verhaak RG, Hoadley KA, Hayes DN, 
#' Perou CM, Schmidt HK, Ding L, Wilson RK, Van Den Berg D, Shen H, 
#' Bengtsson H, Neuvial P, Cope LM, Buckley J, Herman JG, Baylin SB, Laird 
#' PW, Aldape K; Cancer Genome Atlas Research Network. Identification of a 
#' CpG island methylator phenotype that defines a distinct subgroup of 
#' glioma. Cancer Cell. 2010 May 18;17(5):510-22. doi: 
#' 10.1016/j.ccr.2010.03.017. Epub 2010 Apr 15. PubMed PMID: 20399149; PubMed
#' Central PMCID: PMC2872684
#' @keywords datasets
#'
NULL