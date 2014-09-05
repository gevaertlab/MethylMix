\name{MethylMix}
\alias{MethylMix}
\title{Mixture model for DNA methylation data in cancer.}
\description{
     MethylMix identifies DNA methylation driven genes by modeling DNA 
     methylation data in cancer vs. normal and looking for homogeneous 
     subpopulations. In addition matched gene expression data (e.g. from 
     microarray technology or RNA sequencing) is used to identify functional 
     DNA methylation events by requiring a negative correlation between 
     methylation and gene expression of a particular gene. 
}
\usage{
     MethylMix(METcancer,METnormal,MAcancer,OutputRoot,Parallel)
}
\arguments{
     \item{METcancer}{This a matrix with the methylation data of cancer tissue 
     with genes in rows and samples in columns}
     \item{METnormal}{This is a matrix with the normal methylation data of the 
     same genes as in METcancer. Again genes in rows and samples in columns. The
     samples do not have to match with the cancer data.}
     \item{MAcancer}{This is the matched gene expression data for the same 
     samples as in METcancer.}
     \item{OutputRoot}{Path to store the MethylMix results object.}
     \item{Parallel}{If true MethylMix will be run using parallel processing.}
}

\value{
     MethylMixResults is a list with the following components:
	\item{MethylationStates}{Matrix with for all genes the Methylation states 
     using DM-value (i.e. Differential methylation values) that are defined as 
     the methylation value with respect to the average normal methylation for a 
     gene.}
	\item{NrComponents}{The number of methylation states for each gene.}
	\item{Models}{Beta mixture model parameters for each gene.}
	\item{MethylationDrivers}{Genes identified as functional and differential 
     by MethylMix.}
	\item{MixtureStates}{A list with the DM-values for each gene that is 
     functional and differential.}
}

\examples{

# load the three data sets needed for MethylMix
data(METcancer)
data(METnormal)
data(MAcancer)

# run methylmix on a small set of example data
MethylMixResults=MethylMix(METcancer,METnormal,MAcancer)

# try the parallel toolbox to speed up MethylMix modeling
MethylMixResults=MethylMix(METcancer,METnormal,MAcancer,Parallel=TRUE)

}
