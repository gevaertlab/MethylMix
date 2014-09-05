\name{MethylMix_PlotModel}
\alias{MethylMix_PlotModel}
\title{Plotting a mixture model for a gene.}
\description{
     MethylMix identifies DNA methylation driven genes by modeling DNA 
     methylation data in cancer vs. normal and looking for homogeneous 
     subpopulations. In addition matched gene expression data (e.g. from 
     microarray technology or RNA sequencing) is used to identify functional DNA
     methylation events by requiring a negative correlation between methylation 
     and gene expression of a particular gene. 
}
\usage{
     MethylMix_PlotModel(GeneName,METdata,MixtureModelResults,MAdata=0,
     METnormal=0,FileName="")
}
\arguments{
     \item{GeneName}{Name of the gene for which to create a MethylMix plot.}
     \item{METdata}{This a matrix with the methylation data of cancer tissue 
     with genes in rows and samples in columns}
     \item{MixtureModelResults}{The results object from a MethylMix run.}
     \item{METnormal}{This is a matrix with the normal methylation data of the 
     same genes as in METcancer. Again genes in rows and samples in columns. The 
     samples do not have to match with the cancer data.}
     \item{MAdata}{This is the matched gene expression data for the same samples 
     as in METcancer.}
     \item{FileName}{Filename to export the figure. If empty figure is shown in 
     console.}
}

\examples{
     
	# load the three data sets needed for MethylMix
	data(METcancer)
	data(METnormal)
	data(MAcancer)

	# run methylmix on a small set of example data
	MethylMixResults=MethylMix(METcancer,METnormal,MAcancer)

	# Plot the most famous methylated gene for glioblastoma
	MethylMix_PlotModel('MGMT',METcancer,MethylMixResults)
	
	# plot MGMT also with its normal methylation variation
	MethylMix_PlotModel('MGMT',METcancer,MethylMixResults,MAdata=0,METnormal)

	# plot a MethylMix model for another gene
	MethylMix_PlotModel('ZNF217',METcancer,MethylMixResults,MAdata=0,METnormal)

	# also plot the inverse correlation with gene expression 
	# this creates two separate plots
	MethylMix_PlotModel('MGMT',METcancer,MethylMixResults,MAdata=MAcancer)
	
	# plot all functional and differential genes
	for (i in 1:length(MethylMixResults$MethylationDrivers)) {
     	MethylMix_PlotModel(MethylMixResults$MethylationDrivers[i],METcancer,
          MethylMixResults,MAdata=0,METnormal)
	}	
}