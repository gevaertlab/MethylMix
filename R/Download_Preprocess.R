#' The GetData function
#' 
#' This function wraps the functions for downloading and pre-processing DNA methylation and gene expression data,
#' as well as for clustering CpG probes.
#' @param cancerSite character of length 1 with TCGA cancer code.
#' @param targetDirectory character with directory where a folder for downloaded files will be created.
#' @export
#' @keywords download preprocess cluster
#' @details
#' Pre-process of DNA methylation data includes eliminating samples and genes with too many NAs, imputing NAs, and doing Batch correction.
#' If there is both 27k and 450k data, and both data sets have more than 50 samples, we combine the data sets, by reducing the 450k data to the probes present in the 27k data, and bath correction is performed again to the combined data set.
#' If there are samples with both 27k and 450k data, the 450k data is used and the 27k data is discarded, before the step mentioned above.
#' If the 27k or the 450k data does not have more than 50 samples, we use the one with the greatest number of samples, we do not combine the data sets.
#' 
#' For gene expression, this function downloads RNAseq data (file tag "mRNAseq_Preprocess.Level_3"), with the exception for OV and GBM, for which micro array data is
#' downloaded since there is not enough RNAseq data. Pre-process of gene expression data includes eliminating samples and genes with too many NAs, imputing NAs, and doing Batch correction.
#' 
#' For the clustering of the CpG probes, this function uses the annotation for Illumina methylation arrays to map each probe to a gene. Then, for each gene,
#' it clusters all its CpG sites using hierchical clustering and Pearson correlation as distance and complete linkage. 
#' If data for normal samples is provided, only overlapping probes between cancer and normal samples are used. 
#' Probes with SNPs are removed. 
#' 
#' This function is prepared to run in parallel if the user registers a parallel structure, otherwise it runs sequentially.
#' 
#' This function also cleans up the sample names, converting them to the 12 digit format.
#' @return The following files will be created in target directory:
#' \itemize{
#'  \item \code{gdac}: a folder with the raw data downloaded from TCGA.
#'  \item \code{MET_CancerSite_Processed.rds}: processed methylation data at the CpG sites level (not clustered).
#'  \item \code{GE_CancerSite_Processed.rds}: processed gene expression data.
#'  \item \code{data_CancerSite.rds}: list with both gene expression and methylation data. Methylation data is clustered and presented at the gene level. A matrix with the mapping from CpG sites to genes is included.
#' }
#' @examples 
#' \dontrun{
#' # Get data for ovarian cancer
#' cancerSite <- "OV"
#' targetDirectory <- paste0(getwd(), "/")
#' GetData(cancerSite, targetDirectory)
#' 
#' # Optional register cluster to run in parallel
#' library(doParallel)
#' cl <- makeCluster(5)
#' registerDoParallel(cl)
#' 
#' cancerSite <- "OV"
#' targetDirectory <- paste0(getwd(), "/")
#' GetData(cancerSite, targetDirectory)
#' 
#' stopCluster(cl)
#' }
#'
GetData <- function(cancerSite, targetDirectory) {
    
    # Methylation
    cat("Downloading methylation data for:", cancerSite, "\n")
    METdirectories <- Download_DNAmethylation(cancerSite, targetDirectory, TRUE)
    cat("Processing methylation data for:", cancerSite, "\n")
    METProcessedData <- Preprocess_DNAmethylation(cancerSite, METdirectories)
    cat("Saving methylation processed data for:", cancerSite, "\n")
    saveRDS(METProcessedData, file = paste0(targetDirectory, "MET_", cancerSite, "_Processed.rds"))
    
    # Gene expression
    cat("Downloading gene expression data for:", cancerSite, "\n")
    GEdirectories <- Download_GeneExpression(cancerSite, targetDirectory, TRUE)
    cat("Processing gene expression data for:", cancerSite, "\n")
    GEProcessedData <- Preprocess_GeneExpression(cancerSite, GEdirectories)
    cat("Saving gene expression processed data for:", cancerSite, "\n")
    saveRDS(GEProcessedData, file = paste0(targetDirectory, "GE_", cancerSite, "_Processed.rds"))
    
    # Clustering probes to genes
    cat("Clustering methylation data for:", cancerSite, "\n")
    res <- ClusterProbes(METProcessedData[[1]], METProcessedData[[2]])
    
    # Putting everything together
    cat("Saving data\n")
    toSave <- list(METcancer = res[[1]], METnormal = res[[2]], GEcancer = GEProcessedData[[1]], GEnormal = GEProcessedData[[2]], ProbeMapping = res$ProbeMapping)
    saveRDS(toSave, file = paste0(targetDirectory, "data_", cancerSite, ".rds"))
}

#' The Download_DNAmethylation function
#' 
#' Downloads DNA methylation data from TCGA.
#' @param CancerSite character of length 1 with TCGA cancer code.
#' @param TargetDirectory character with directory where a folder for downloaded files will be created.
#' @param downloadData logical indicating if data should be downloaded (default: TRUE). If false, the url of the desired data is returned.
#' @return list with paths to downloaded files for both 27k and 450k methylation data.
#' @export
#' @keywords download
#' @examples 
#' \dontrun{
#' 
#' # Optional register cluster to run in parallel
#' library(doParallel)
#' cl <- makeCluster(5)
#' registerDoParallel(cl)
#' 
#' # Methylation data for ovarian cancer
#' cancerSite <- "OV"
#' targetDirectory <- paste0(getwd(), "/")
#' 
#' # Downloading methylation data
#' METdirectories <- Download_DNAmethylation(cancerSite, targetDirectory, TRUE)
#' 
#' # Processing methylation data
#' METProcessedData <- Preprocess_DNAmethylation(cancerSite, METdirectories)
#' 
#' # Saving methylation processed data
#' saveRDS(METProcessedData, file = paste0(targetDirectory, "MET_", cancerSite, "_Processed.rds"))
#' 
#' # Clustering methylation data
#' res <- ClusterProbes(METProcessedData[[1]], METProcessedData[[2]])
#' 
#' # Saving methylation clustered data
#' toSave <- list(METcancer = res[[1]], METnormal = res[[2]], ProbeMapping = res$ProbeMapping)
#' saveRDS(toSave, file = paste0(targetDirectory, "MET_", cancerSite, "_Clustered.rds"))
#' 
#' stopCluster(cl)
#' }
#'  
Download_DNAmethylation <- function(CancerSite, 
                                    TargetDirectory, 
                                    downloadData = TRUE) {    
    
    dir.create(TargetDirectory,showWarnings=FALSE)
    
    # download the 27k data
    dataType='stddata'
    dataFileTag='Merge_methylation__humanmethylation27'
    cat('Searching 27k MET data for:',CancerSite,'\n')
    METdirectory27k=get_firehoseData(downloadData,TargetDirectory,CancerSite,dataType,dataFileTag)
    
    # download the 450k data
    dataFileTag='Merge_methylation__humanmethylation450'
    cat('Searching 450k MET data for:',CancerSite,'\n')
    METdirectory450k=get_firehoseData(downloadData,TargetDirectory,CancerSite,dataType,dataFileTag)
    
    return(METdirectories=list(METdirectory27k=METdirectory27k,METdirectory450k=METdirectory450k))
}

#' The get_firehoseData function
#' 
#' Gets data from TCGA's firehose.
#' @param downloadData logical indicating if data should be downloaded (default: TRUE). If false, the url of the desired data is returned.
#' @param saveDir path to directory to save downloaded files.
#' @param TCGA_acronym_uppercase TCGA's cancer site code.
#' @param dataType type of data in TCGA (default: "stddata").
#' @param dataFileTag name of the file to be downloaded (the default is to download RNAseq data, but this can be changed to download other data).
#' @param FFPE logical indicating if FFPE data should be downloaded (default: FALSE).
#' @param fileType type of downloaded file (default: "fileType", other type not admitted at the moment).
#' @param gdacURL gdac url.
#' @param untarUngzip logical indicating if the gzip file downloaded should be untarred (default: TRUE).
#' @param printDisease_abbr if TRUE data is not downloaded but all the possible cancer sites codes are shown (default: FALSE).
#' @return DownloadedFile path to directory with downloaded files.
#' @keywords internal
#' 
get_firehoseData <- function(downloadData=TRUE,
                             saveDir = "./",
                             TCGA_acronym_uppercase = "LUAD",
                             dataType="stddata",
                             dataFileTag = "mRNAseq_Preprocess.Level_3",
                             FFPE=FALSE,
                             fileType= "tar.gz",
                             gdacURL= "http://gdac.broadinstitute.org/runs/",
                             untarUngzip=TRUE,
                             printDisease_abbr=FALSE){  
    
    # Cases Shipped by BCR  # Cases with Data*  Date Last Updated (mm/dd/yy)
    cancers <- c("Acute Myeloid Leukemia [LAML] \n","Adrenocortical carcinoma [ACC]	\n",
                 "Bladder Urothelial Carcinoma [BLCA] \n",	"Brain Lower Grade Glioma [LGG] \n",
                 "Breast invasive carcinoma [BRCA] \n","Cervical squamous cell carcinoma and endocervical adenocarcinoma [CESC] \n",
                 "Cholangiocarcinoma [CHOL] \n",	"Colon adenocarcinoma [COAD] \n",	"Esophageal carcinoma [ESCA] \n",
                 "Glioblastoma multiforme [GBM] \n",	"Head and Neck squamous cell carcinoma [HNSC]	\n",
                 "Kidney Chromophobe [KICH]	\n","Kidney renal clear cell carcinoma [KIRC]	\n",
                 "Kidney renal papillary cell carcinoma [KIRP]	\n","Liver hepatocellular carcinoma [LIHC]	\n",
                 "Lung adenocarcinoma [LUAD]	\n", "Lung squamous cell carcinoma [LUSC] \n",
                 "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma [DLBC]	\n","Mesothelioma [MESO] \n",
                 "Ovarian serous cystadenocarcinoma [OV]	\n","Pancreatic adenocarcinoma [PAAD]	\n",
                 "Pheochromocytoma and Paraganglioma [PCPG] \n","Prostate adenocarcinoma [PRAD] \n",
                 "Rectum adenocarcinoma [READ]	\n","Sarcoma [SARC]	\n","Skin Cutaneous Melanoma [SKCM]	\n",
                 "Stomach adenocarcinoma [STAD] \n","Testicular Germ Cell Tumors [TGCT] \n","Thymoma [THYM] \n",
                 "Thyroid carcinoma [THCA]	\n","Uterine Carcinosarcoma [UCS]	 \n",
                 "Uterine Corpus Endometrial Carcinoma [UCEC]	\n","Uveal Melanoma [UVM] \n", "Colorectal Adenocarcinoma [COADREAD] \n")
    
    if(printDisease_abbr){      
        return(cat("here are the possible TCGA database disease acronyms. \nRe-run this function with printDisease_abbr=FALSE to then run an actual query.\n\n",cancers));      
    }
    gdacURL_orig <- gdacURL
    
    # New code to handle dates - Marcos
    gdacURLnew <- paste0(gdacURL_orig, dataType, "__latest/data/", TCGA_acronym_uppercase, "/")
    urlDataNew <- RCurl::getURL(gdacURLnew)
    urlDataNew <- limma::strsplit2(urlDataNew, "href=\\\"") #regular expressions: need \ to have R recognize any " or \ that's actually in our text
    getlatestdate <- urlDataNew[grep("^201[0-9][01][0-9][0123][0-9]", urlDataNew)]
    getlatestdate <- substring(getlatestdate, 1, 8)
    gdacURLnew <- paste0(gdacURLnew, getlatestdate, "/")
    urlData <- RCurl::getURL(gdacURLnew)
    urlData <- limma::strsplit2(urlData, "href=\\\"") #regular expressions: need \ to have R recognize any " or \ that's actually in our text
    lastDateCompress <- lastDate <- getlatestdate # for compatibility with the rest of old code
    gdacURL <- gdacURLnew # for compatibility with the rest of old code
    # end New code
    
    #remove any FFPE datasets, or only keep those depending on user inputs.
    if (FFPE) { 
        urlData <- urlData[grep("FFPE",urlData)]	
        if(length(urlData)==0){		
            stop("\nNo FFPE data found for this query. Try FFPE=FALSE.\n")		
        }	  
    } else {	
        #we DON'T want FFPE data.
        #but if no FFPE data to begin with: don't subset on this.
        if(length(grep("FFPE",urlData))>0){		
            urlData <- urlData[-grep("FFPE",urlData)]		
        }
        if(length(urlData)==0){		
            stop("\nNo non-FFPE data found for this query. Try FFPE=TRUE.\n")		
        }
    }
    #now get full dataset name.
    fileName <- urlData[grep(dataFileTag,urlData)]
    
    if(length(fileName)==0){	  
        #warnMessage <- paste0("\nNot returning any viable url data paths after searching by date for disease ",TCGA_acronym_uppercase," for data type ",dataFileTag ,".No data was downloaded.\n")
        #warning(warnMessage)
        cat("\tThere is no",dataFileTag,"data for",TCGA_acronym_uppercase,"\n")
        return(NA)	  
    }
    #some redundancy..but that' OK because we'll add back on the unique tar.gz file tag.
    #first file is one we want - not md5 file.
    fileName <- limma::strsplit2(fileName,"tar.gz")[1,1]
    fileName <- paste(fileName,fileType,sep="")
    
    #final download url
    gdacURL <- paste(gdacURL,fileName,sep="")
    # Directory for downloads
    saveDir <- paste(saveDir,"gdac_",lastDateCompress,'/',sep="")
    
    if (!grepl("Windows", Sys.info()['sysname'])) {
        # Not Windows
        tarfile=paste0(saveDir,fileName)
        finalDir <-  strsplit(tarfile, paste0(".", fileType))[[1]][1]
        if(downloadData){		
            cat("\tDownloading",dataFileTag,"data, version:",lastDate,"\n")				
            cat("\tThis may take 10-60 minutes depending on the size of the data set.\n")
            dir.create(saveDir,showWarnings=FALSE)
            # download file		
            setwd(saveDir)				
            download.file(gdacURL,fileName,quiet=FALSE,mode="wb")
            #this assumes a tar.gz file.
            if(fileType=="tar.gz" && untarUngzip) {					
                cat("\tUnpacking data.\n")
                tarfile=paste0(saveDir,fileName)
                untar(tarfile)
                #remove tarred file
                fileToRemove <- limma::strsplit2(gdacURL,"/")[ ,ncol(limma::strsplit2(gdacURL,"/"))]
                removed <- file.remove(paste0(saveDir,fileToRemove))
            } else if(untarUngzip) {		
                warning("File expansion/opening only built in for tar.gz files at the moment.\n")		
            }		
            cat("\tFinished downloading",dataFileTag,"data to",finalDir,"\n")
        } else {
            cat("\tdownload data url is :\n ",gdacURL,'\n')
        }
        DownloadedFile=paste0(finalDir,'/')
        return(DownloadedFile)
    } else {
        # new code to handle long names in windows - Marcos
        # WINDOWS: name of file can be too long (and it's repeated in the folder and the file) and windows
        # internally will use another representation for the name of the folder, so then when
        # we want to load the file it says it doesn't exist. So I'm changing the name of the folder
        # to prevent this. We can't change the name of the file as it's used in the Preprocess functions
        idx <- which(sapply(c("methylation27", "methylation450", "mRNAseq", "transcriptome"), grepl, dataFileTag))[1]
        newtag <- ifelse(idx == 1, "meth27", ifelse(idx == 2, "meth450", "geneexpr"))
        nameForFolder <- paste(TCGA_acronym_uppercase, dataType, newtag, sep = "_")
        nameForDownloadedFile <- paste0(nameForFolder, ".", fileType)
        nameForDownloadedFileFullPath <- paste0(saveDir, nameForDownloadedFile)
        finalDir <- paste0(saveDir, nameForFolder)
        if (downloadData) {
            cat("\tDownloading",dataFileTag,"data, version:",lastDate,"\n")				
            cat("\tThis may take 10-60 minutes depending on the size of the data set.\n")
            dir.create(saveDir, showWarnings = FALSE)
            
            # Create a virtual drive to overcome long names issue in Windows
            saveDir2 <- gsub("\\\\", "/", saveDir)
            saveDir2 <- substr(saveDir2, 1, nchar(saveDir2) - 1)
            system(paste("subst x:", saveDir2))
            
            download.file(gdacURL, destfile = paste0("x://", nameForDownloadedFile), quiet = FALSE, mode = "wb")
            #this assumes a tar.gz file.
            if(fileType == "tar.gz" && untarUngzip) {					
                cat("\tUnpacking data.\n")
                untar(nameForDownloadedFileFullPath, exdir = saveDir)
                # untar(paste0("x://", nameForDownloadedFile), exdir = "x://") # doesn't work because it calls system and system doesnt know about the virtual drive
                removed <- file.remove(paste0("x://", nameForDownloadedFile))
                # Anyway I change folder name to make it shorter
                changed <- file.rename(from = paste0("x://", gsub(".tar.gz", "", fileName)), to = paste0("x://", nameForFolder))
                system("subst x: /D") # stop the virtual drive
            } else if(untarUngzip) {		
                warning("File expansion/opening only built in for tar.gz files at the moment.\n")		
            }
            cat("\tFinished downloading", dataFileTag, "data to", finalDir,"\n")
        } else {
            cat("\tdownload data url is :\n ", gdacURL, '\n')
        }
        DownloadedFile = paste0(finalDir, '/')
        return(DownloadedFile)
        # end new code
    }
}


#' The Preprocess_DNAmethylation function
#' 
#' Pre-processes DNA methylation data from TCGA.
#' @param CancerSite character of length 1 with TCGA cancer code.
#' @param METdirectories character vector with directories with the downloaded data. It can be the object returned by the Download_DNAmethylation function.
#' @param MissingValueThreshold threshold for removing samples or genes with missing values.
#' @details 
#' Pre-process includes eliminating samples and genes with too many NAs, imputing NAs, and doing Batch correction.
#' If there is both 27k and 450k data, and both data sets have more than 50 samples, we combine the data sets, by reducing the 450k data to the probes present in the 27k data, and bath correction is performed again to the combined data set.
#' If there are samples with both 27k and 450k data, the 450k data is used and the 27k data is discarded, before the step mentioned above.
#' If the 27k or the 450k data does not have more than 50 samples, we use the one with the greatest number of samples, we do not combine the data sets.
#' @return List with the pre-processed data matrix for cancer and normal samples.
#' @export
#' @keywords preprocess
#' @examples 
#' \dontrun{
#' 
#' # Optional register cluster to run in parallel
#' library(doParallel)
#' cl <- makeCluster(5)
#' registerDoParallel(cl)
#' 
#' # Methylation data for ovarian cancer
#' cancerSite <- "OV"
#' targetDirectory <- paste0(getwd(), "/")
#' 
#' # Downloading methylation data
#' METdirectories <- Download_DNAmethylation(cancerSite, targetDirectory, TRUE)
#' 
#' # Processing methylation data
#' METProcessedData <- Preprocess_DNAmethylation(cancerSite, METdirectories)
#' 
#' # Saving methylation processed data
#' saveRDS(METProcessedData, file = paste0(targetDirectory, "MET_", cancerSite, "_Processed.rds"))
#' 
#' # Clustering methylation data
#' res <- ClusterProbes(METProcessedData[[1]], METProcessedData[[2]])
#' 
#' # Saving methylation clustered data
#' toSave <- list(METcancer = res[[1]], METnormal = res[[2]], ProbeMapping = res$ProbeMapping)
#' saveRDS(toSave, file = paste0(targetDirectory, "MET_", cancerSite, "_Clustered.rds"))
#' 
#' stopCluster(cl)
#' }
#' 
Preprocess_DNAmethylation <- function(CancerSite, METdirectories, MissingValueThreshold = 0.2) {    
    
    cat("\tProcessing data for",CancerSite,"\n")
    ProcessedData27k=c()
    ProcessedData450k=c()
    if (!is.na(METdirectories$METdirectory27k)) {
        cat('\tLoading data for 27k.\n')
        ProcessedData27k=Preprocess_CancerSite_Methylation27k(CancerSite,METdirectories$METdirectory27k, MissingValueThreshold)
    }
    
    if (!is.na(METdirectories$METdirectory450k)) {
        cat('\tLoading data for 450k.\n')
        ProcessedData450k=Preprocess_CancerSite_Methylation450k(CancerSite,METdirectories$METdirectory450k, MissingValueThreshold)
    }
    
    # check if we want to combine 27k and 450k
    if (length(ProcessedData27k)!=0 & length(ProcessedData450k)!=0 && CancerSite !="LAML") {
        # only do it when enough samples are in both
        if (ncol(ProcessedData27k$MET_Data_Cancer)>50 & ncol(ProcessedData450k$MET_Data_Cancer)>50) {               
            Mode='450kon27k'               
            cat("\tCombining 450k and 27k by mapping 450k probes to 27k array.\n")
            # Check if there are duplicate samples, remove the 27k ones. 
            OverlapSamplesCancer=intersect(colnames(ProcessedData27k$MET_Data_Cancer),colnames(ProcessedData450k$MET_Data_Cancer))
            if (length(OverlapSamplesCancer)>0) {
                cat("\tCancer sample overlap is not empty: ",length(OverlapSamplesCancer))
                ProcessedData27k$MET_Data_Cancer=ProcessedData27k$MET_Data_Cancer[,-OverlapSamplesCancer,drop=FALSE]
            }
            OverlapSamplesNormal=intersect(colnames(ProcessedData27k$MET_Data_Normal),colnames(ProcessedData450k$MET_Data_Normal))
            if (length(OverlapSamplesNormal)>0) {
                cat("\tNormal sample overlap is not empty: ",length(OverlapSamplesNormal))
                ProcessedData27k$MET_Data_Normal=ProcessedData27k$MET_Data_Normal[,-OverlapSamplesNormal,drop=FALSE]
            }
            
            # Overlap the probes
            ProcessedData=list(MET_Data_Cancer=c(),MET_Data_Normal=c())
            OverlapProbesCancer=intersect(rownames(ProcessedData27k$MET_Data_Cancer),rownames(ProcessedData450k$MET_Data_Cancer))               
            ProcessedData$MET_Data_Cancer=cbind(ProcessedData27k$MET_Data_Cancer[OverlapProbesCancer,,drop=FALSE],ProcessedData450k$MET_Data_Cancer[OverlapProbesCancer,,drop=FALSE])
            if ( length(colnames(ProcessedData27k$MET_Data_Normal))>0 & length(colnames(ProcessedData450k$MET_Data_Normal))>0 ) {
                OverlapProbesNormal=intersect(rownames(ProcessedData27k$MET_Data_Normal),rownames(ProcessedData450k$MET_Data_Normal))               
                ProcessedData$MET_Data_Normal=cbind(ProcessedData27k$MET_Data_Normal[OverlapProbesNormal,,drop=FALSE],ProcessedData450k$MET_Data_Normal[OverlapProbesNormal,,drop=FALSE])
            } else if ( length(colnames(ProcessedData27k$MET_Data_Normal))>0 ) {
                ProcessedData$MET_Data_Normal=ProcessedData27k$MET_Data_Normal                    
            } else if ( length(colnames(ProcessedData450k$MET_Data_Normal))>0 ) {
                ProcessedData$MET_Data_Normal=ProcessedData450k$MET_Data_Normal                    
            } 
            
            # Batch correction on combined Tumor data.
            Batch=matrix(1,length(colnames(ProcessedData$MET_Data_Cancer)),1)
            Batch[1:length(colnames(ProcessedData27k$MET_Data_Cancer)),1]=2
            BatchData=data.frame(ArrayName=colnames(ProcessedData$MET_Data_Cancer),SampleName=colnames(ProcessedData$MET_Data_Cancer),Batch=Batch)
            ProcessedData$MET_Data_Cancer=TCGA_BatchCorrection_MolecularData(ProcessedData$MET_Data_Cancer,BatchData,0)
            
            if (length(colnames(ProcessedData$MET_Data_Normal))>0) {
                # Batch correction on combined Normal data.
                Batch=matrix(1,length(colnames(ProcessedData$MET_Data_Normal)),1)
                Batch[1:length(colnames(ProcessedData27k$MET_Data_Normal)),1]=2
                BatchData=data.frame(ArrayName=colnames(ProcessedData$MET_Data_Normal),SampleName=colnames(ProcessedData$MET_Data_Normal),Batch=Batch)
                ProcessedData$MET_Data_Normal=TCGA_BatchCorrection_MolecularData(ProcessedData$MET_Data_Normal,BatchData,5)
            }               
            
        } else if (ncol(ProcessedData27k$MET_Data_Cancer)>ncol(ProcessedData450k$MET_Data_Cancer)) { 
            cat("\tNot enough 450k samples, only using the 27k (need min 50 samples).\n")
            Mode='27k'
            ProcessedData=ProcessedData27k          
        } else {
            cat("\tNot enough 27k samples, only using the 450k (need min 50 samples).\n")
            Mode='450k'
            ProcessedData=ProcessedData450k     
        }
    } else if (CancerSite == "LAML") {
        cat("\tLAML is a special case, only using 450k data.\n")
        OverlapSamplesCancer=intersect(colnames(ProcessedData27k$MET_Data_Cancer),colnames(ProcessedData450k$MET_Data_Cancer))
        cat("\tOverlap length is:",length(OverlapSamplesCancer),".\n")
        Mode='450k'          
        ProcessedData=ProcessedData450k          
    } else if (length(ProcessedData27k)!=0) {
        cat("\tOnly 27k samples.\n")
        Mode='27k'
        ProcessedData=ProcessedData27k          
    } else {       
        cat("\tOnly 450k samples.\n")
        Mode='450k'
        ProcessedData=ProcessedData450k          
    }     
    return(ProcessedData=ProcessedData)
}

#' The Preprocess_CancerSite_Methylation27k function
#' 
#' Internal. Pre-processes DNA methylation data from TCGA from Illymina 27k arrays.
#' @param CancerSite character of length 1 with TCGA cancer code.
#' @param METdirectory character with directory where a folder for downloaded files will be created. Can be the object returned by the Download_DNAmethylation function.
#' @param MissingValueThreshold threshold for removing samples or genes with missing values.
#' @return List with pre processed methylation data for cancer and normal samples.
#' @keywords internal
#'
Preprocess_CancerSite_Methylation27k <- function(CancerSite, METdirectory, MissingValueThreshold = 0.2) {
    
    # Settings
    get("BatchData")
    MinPerBatch=5
    
    if (grepl("Windows", Sys.info()['sysname'])) {
        # If Windows I'll create a virtual drive to handle the long file names issue
        # Create a virtual drive to overcome long names issue in Windows
        virtualDir <- METdirectory
        virtualDir <- gsub("\\\\", "/", virtualDir)
        virtualDir <- substr(virtualDir, 1, nchar(virtualDir) - 1)
        system(paste("subst x:", virtualDir))
        
        # Load data
        METfiles <- dir("x:")
        MatchedFilePosition <- grep('methylation__humanmethylation27', METfiles)             
        Filename <- paste0("x://", METfiles[MatchedFilePosition])          
        MET_Data <- TCGA_GENERIC_LoadIlluminaMethylationData(Filename)
        
        system("subst x: /D") #stop virtual drive
    } else {
        # Not windows
        # Load data
        METfiles=dir(METdirectory)
        MatchedFilePosition=grep('methylation__humanmethylation27',METfiles)             
        Filename=paste0(METdirectory,METfiles[MatchedFilePosition])          
        MET_Data=TCGA_GENERIC_LoadIlluminaMethylationData(Filename)
    }
    
    # Split up normal and cancer data
    Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(MET_Data))
    if (CancerSite=='LAML') {
        MET_Data_Cancer=MET_Data[,Samplegroups$PeripheralBloodCancer,drop=FALSE]
    } else {
        MET_Data_Cancer=MET_Data[,Samplegroups$Primary,drop=FALSE]
    }
    MET_Data_Cancer=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Cancer,15)
    if (CancerSite=='LAML') {
        MET_Data_Normal=MET_Data[,Samplegroups$BloodNormal,drop=FALSE]
    } else {
        MET_Data_Normal=MET_Data[,Samplegroups$SolidNormal,drop=FALSE]
    }
    if (length(MET_Data_Normal)>0) {
        MET_Data_Normal=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Normal,15)
    }
    cat("There are",length(colnames(MET_Data_Cancer)),"cancer samples and",length(colnames(MET_Data_Normal)),"normal samples in",CancerSite,"\n")
    
    # Missing value estimation
    cat("\tMissing value estimation for the cancer samples.\n")
    MET_Data_Cancer=TCGA_Process_EstimateMissingValues(MET_Data_Cancer,MissingValueThreshold)
    if (length(MET_Data_Normal)>0) {
        cat("\tMissing value estimation for the normal samples.\n")
        MET_Data_Normal=TCGA_Process_EstimateMissingValues(MET_Data_Normal,MissingValueThreshold)
    }
    
    # Batch correction for cancer and normal. 
    cat("\tBatch correction for the cancer samples.\n")
    BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Cancer,BatchData)
    MET_Data_Cancer=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer,BatchData,MinPerBatch)
    
    if (length(MET_Data_Normal)>0) {
        cat("\tBatch correction for the normal samples.\n")
        BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Normal,BatchData)
        MET_Data_Normal=TCGA_BatchCorrection_MolecularData(MET_Data_Normal,BatchData,MinPerBatch)
    } else {
        MET_Data_Normal=c()
    }
    
    # Reducing to 12 ids. 
    #MET_Data_Cancer=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Cancer,12)
    #MET_Data_Normal=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Normal,12)     
    
    MET_Data_Cancer[MET_Data_Cancer<0]=0
    MET_Data_Cancer[MET_Data_Cancer>1]=1
    if (length(MET_Data_Normal)>0) {
        MET_Data_Normal[MET_Data_Normal<0]=0
        MET_Data_Normal[MET_Data_Normal>1]=1
    }
    
    return(list(MET_Data_Cancer=MET_Data_Cancer,MET_Data_Normal=MET_Data_Normal))
}

#' The Preprocess_CancerSite_Methylation450k function
#' 
#' Internal. Pre-processes DNA methylation data from TCGA from Illymina 450k arrays.
#' @param CancerSite character of length 1 with TCGA cancer code.
#' @param METdirectory character with directory where a folder for downloaded files will be created. Can be the object returned by the Download_DNAmethylation function.
#' @param MissingValueThreshold threshold for removing samples or genes with missing values.
#' @return List with pre processed methylation data for cancer and normal samples.
#' @keywords internal
#'
Preprocess_CancerSite_Methylation450k <- function(CancerSite, METdirectory, MissingValueThreshold = 0.2) {
    
    get("BatchData")
    MinPerBatch=5
    
    if (grepl("Windows", Sys.info()['sysname'])) {
        # If Windows I'll create a virtual drive to handle the long file names issue
        # Create a virtual drive to overcome long names issue in Windows
        virtualDir <- METdirectory
        virtualDir <- gsub("\\\\", "/", virtualDir)
        virtualDir <- substr(virtualDir, 1, nchar(virtualDir) - 1)
        system(paste("subst x:", virtualDir))
        
        # Load data
        METfiles <- dir("x:")
        MatchedFilePosition <- grep('methylation__humanmethylation450', METfiles)             
        Filename <- paste0("x://", METfiles[MatchedFilePosition])          
        MET_Data <- TCGA_GENERIC_LoadIlluminaMethylationData(Filename)
        
        system("subst x: /D") #stop virtual drive
    } else {
        # Not windows
        # Load data
        METfiles=dir(METdirectory)
        MatchedFilePosition=grep('methylation__humanmethylation450',METfiles)             
        Filename=paste0(METdirectory,METfiles[MatchedFilePosition])          
        MET_Data=TCGA_GENERIC_LoadIlluminaMethylationData(Filename)
    }
    
    # Split up normal and cancer data
    Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(MET_Data))
    if (CancerSite=='LAML') {
        MET_Data_Cancer=MET_Data[,Samplegroups$PeripheralBloodCancer]
    } else {
        MET_Data_Cancer=MET_Data[,Samplegroups$Primary]
    }
    MET_Data_Cancer=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Cancer,15)
    if (CancerSite=='LAML') {
        MET_Data_Normal=MET_Data[,Samplegroups$BloodNormal]
    } else {
        MET_Data_Normal=MET_Data[,Samplegroups$SolidNormal]
    }     
    if (length(MET_Data_Normal)>0) {
        MET_Data_Normal=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Normal,15)
    }
    cat("There are",length(colnames(MET_Data_Cancer)),"cancer samples and",length(colnames(MET_Data_Normal)),"normal samples in",CancerSite,"\n")
    
    # Clear space
    rm(MET_Data); gc()
    
    # Missing value estimation
    cat("\tMissing value estimation for the cancer samples.\n")
    MET_Data_Cancer=TCGA_Process_EstimateMissingValues(MET_Data_Cancer,MissingValueThreshold)     
    if (length(MET_Data_Normal)>0) {
        cat("\tMissing value estimation for the normal samples.\n")
        MET_Data_Normal=TCGA_Process_EstimateMissingValues(MET_Data_Normal,MissingValueThreshold)
    }     
    
    # Split up cancer data set
    MET_Data_Cancer1=MET_Data_Cancer[1:50000,]
    MET_Data_Cancer2=MET_Data_Cancer[50001:100000,]
    MET_Data_Cancer3=MET_Data_Cancer[100001:150000,]
    MET_Data_Cancer4=MET_Data_Cancer[150001:200000,]
    MET_Data_Cancer5=MET_Data_Cancer[200001:250000,]
    MET_Data_Cancer6=MET_Data_Cancer[250001:300000,]
    MET_Data_Cancer7=MET_Data_Cancer[300001:350000,]
    MET_Data_Cancer8=MET_Data_Cancer[350001:nrow(MET_Data_Cancer),]
    
    # clearing some memory
    rm(MET_Data_Cancer); gc()
    
    # Split up normal data set
    if (length(MET_Data_Normal)>0) {
        MET_Data_Normal1=MET_Data_Normal[1:100000,,drop=FALSE]
        MET_Data_Normal2=MET_Data_Normal[100001:200000,,drop=FALSE]
        MET_Data_Normal3=MET_Data_Normal[200001:300000,,drop=FALSE]
        MET_Data_Normal4=MET_Data_Normal[300001:nrow(MET_Data_Normal),,drop=FALSE]
        # clearing some memory
        rm(MET_Data_Normal); gc()
    } else {
        MET_Data_Normal1=c()
    }
    
    cat("\tBatch correction for the cancer samples.\n")
    MET_Data_Cancer1=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer1,BatchData,MinPerBatch)
    MET_Data_Cancer2=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer2,BatchData,MinPerBatch)
    MET_Data_Cancer3=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer3,BatchData,MinPerBatch)
    MET_Data_Cancer4=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer4,BatchData,MinPerBatch)
    MET_Data_Cancer5=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer5,BatchData,MinPerBatch)
    MET_Data_Cancer6=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer6,BatchData,MinPerBatch)
    MET_Data_Cancer7=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer7,BatchData,MinPerBatch)
    MET_Data_Cancer8=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer8,BatchData,MinPerBatch)
    
    if (length(MET_Data_Normal1)>0) {
        cat("\tBatch correction for the normal samples.\n")
        MET_Data_Normal1=TCGA_BatchCorrection_MolecularData(MET_Data_Normal1,BatchData,2)
        MET_Data_Normal2=TCGA_BatchCorrection_MolecularData(MET_Data_Normal2,BatchData,2)
        MET_Data_Normal3=TCGA_BatchCorrection_MolecularData(MET_Data_Normal3,BatchData,2)
        MET_Data_Normal4=TCGA_BatchCorrection_MolecularData(MET_Data_Normal4,BatchData,2)
    }
    
    # Combine batch corrected data
    # It's possible that samples get deleted due to missings in one part and not another.
    OverlapSamples=Reduce(intersect,list(colnames(MET_Data_Cancer1),colnames(MET_Data_Cancer2),colnames(MET_Data_Cancer3),colnames(MET_Data_Cancer4)
                                         ,colnames(MET_Data_Cancer5),colnames(MET_Data_Cancer6),colnames(MET_Data_Cancer7),colnames(MET_Data_Cancer8)))
    MET_Data_Cancer=rbind(MET_Data_Cancer1[,OverlapSamples],MET_Data_Cancer2[,OverlapSamples],MET_Data_Cancer3[,OverlapSamples]
                          ,MET_Data_Cancer4[,OverlapSamples],MET_Data_Cancer5[,OverlapSamples],MET_Data_Cancer6[,OverlapSamples]
                          ,MET_Data_Cancer7[,OverlapSamples],MET_Data_Cancer8[,OverlapSamples])
    
    rm(MET_Data_Cancer1);rm(MET_Data_Cancer2);rm(MET_Data_Cancer3);rm(MET_Data_Cancer4);rm(MET_Data_Cancer5);rm(MET_Data_Cancer6);rm(MET_Data_Cancer7);rm(MET_Data_Cancer8);gc()
    
    if (length(MET_Data_Normal1)>0) {
        OverlapSamples=Reduce(intersect,list(colnames(MET_Data_Normal1),colnames(MET_Data_Normal2),colnames(MET_Data_Normal3),colnames(MET_Data_Normal4)))
        MET_Data_Normal=rbind(MET_Data_Normal1[,OverlapSamples,drop=FALSE],MET_Data_Normal2[,OverlapSamples,drop=FALSE],
                              MET_Data_Normal3[,OverlapSamples,drop=FALSE],MET_Data_Normal4[,OverlapSamples,drop=FALSE])
        rm(MET_Data_Normal1);rm(MET_Data_Normal2); rm(MET_Data_Normal3);rm(MET_Data_Normal4);gc()
    } else {
        MET_Data_Normal=c()
    }
    
    
    
    # Set values <0 to 0 and >1 to 1, because of batch correction
    MET_Data_Cancer[MET_Data_Cancer<0]=0
    MET_Data_Cancer[MET_Data_Cancer>1]=1
    if (length(MET_Data_Normal)>0) {
        MET_Data_Normal[MET_Data_Normal<0]=0
        MET_Data_Normal[MET_Data_Normal>1]=1
    }
    
    return(list(MET_Data_Cancer=MET_Data_Cancer,MET_Data_Normal=MET_Data_Normal))
}

#' The TCGA_GENERIC_LoadIlluminaMethylationData function
#' 
#' Internal. Read in an illumina methylation file with the following format: header row with sample labels, 
#' 2nd header row with 4 columns per sample: beta-value, geneSymbol, chromosome and GenomicCoordinate. 
#' The first column has the probe names. 
#' @param Filename name of the file with the data.
#' @return methylation data.
#' @keywords internal
#'
TCGA_GENERIC_LoadIlluminaMethylationData <- function(Filename) {
    
    # read in an illumina methylation file with the following format: 
    # header row with sample labels
    # 2nd header row with 4 columns per sample: beta-value, geneSymbol, chromosome and GenomicCoordinate
    # The first column has the probe names. 
    MET_Data<-data.table::fread(Filename)
    MET_Data=as.matrix(MET_Data)
    Probes=MET_Data[,1]
    rownames(MET_Data)=Probes
    MET_Data=MET_Data[,-1]
    MET_Data=MET_Data[-1,]
    MET_Data=MET_Data[,seq(1,ncol(MET_Data),4)]     
    class(MET_Data)='numeric'
    
    return(MET_Data)
}

#' The TCGA_GENERIC_GetSampleGroups function
#' 
#' Internal. Looks for the group of the samples (normal/cancer).
#' @param SampleNames vector with sample names.
#' @return a list.
#' @keywords internal
#'
TCGA_GENERIC_GetSampleGroups <-function(SampleNames) {
    
    # First replace any . with - so the sample groups are uniform. 
    SampleGroups=list()
    
    #1: Primary Tumor
    Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]01[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
    SampleGroups$Primary=SampleNames[Matches==1]
    
    #2: Recurrent tumor
    Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]02[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
    SampleGroups$Recurrent=SampleNames[Matches==1]
    
    #3: Primary blood derived cancer - peripheral blood
    Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]03[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
    SampleGroups$PeripheralBloodCancer=SampleNames[Matches==1]     
    
    #10: Blood derived normal
    Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]10[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
    SampleGroups$BloodNormal=SampleNames[Matches==1]     
    
    #11: Solid tissue derived normal
    Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]11[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
    SampleGroups$SolidNormal=SampleNames[Matches==1]     
    
    #20 Cellines
    Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]20[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
    SampleGroups$CellLines=SampleNames[Matches==1]    
    
    #06 Cellines
    Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]06[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
    SampleGroups$Metastatic=SampleNames[Matches==1]         
    
    return(SampleGroups)
}

#' The TCGA_GENERIC_CleanUpSampleNames function
#' 
#' Internal. Cleans the samples IDs into the 12 digit format and removes doubles.
#' @param GEN_Data data matrix.
#' @param IDlength length of samples ID.
#' @return data matrix with cleaned sample names.
#' @keywords internal
#'
TCGA_GENERIC_CleanUpSampleNames <-function(GEN_Data, IDlength = 12) {     
    SampleNames=colnames(GEN_Data)
    SampleNamesShort=as.character(apply(as.matrix(SampleNames),2,substr,1,IDlength))
    if (length(SampleNamesShort)!=length(unique(SampleNamesShort))) {
        # remove the doubles           
        Counts=table(SampleNamesShort)
        Doubles=rownames(Counts)[which(Counts>1)]
        
        cat("Removing doubles for",length(Doubles),"samples.\n")
        for(i in 1:length(Doubles)) {                         
            CurrentDouble=Doubles[i]          
            pos=grep(CurrentDouble,SampleNames)
            #GEN_Data[1:10,pos]
            #cor(GEN_Data[,pos])
            GEN_Data=GEN_Data[,-pos[2:length(pos)]]     
            SampleNames=colnames(GEN_Data) # need to update samplenames because pos is relative to this
        }
        SampleNames=colnames(GEN_Data)
        SampleNamesShort=as.character(apply(as.matrix(SampleNames),2,substr,1,IDlength))
        
        # now set the samplenames
        colnames(GEN_Data)=SampleNamesShort
    } else {
        colnames(GEN_Data)=SampleNamesShort     
    }     
    return(GEN_Data)
}

#' The TCGA_Process_EstimateMissingValues function
#' 
#' Internal. Removes patients and genes with more missing values than the MissingValueThreshold, and imputes remaining missing values using Tibshirani's KNN method.
#' @param MET_Data data matrix.
#' @param MissingValueThreshold threshold for removing samples and genes with too many missing values.
#' @return the data set with imputed values and possibly some genes or samples deleted.
#' @keywords internal
#'
TCGA_Process_EstimateMissingValues <- function(MET_Data, MissingValueThreshold = 0.2) {
    
    # FIRST REMOVING BAD PATIENTS
    # removing patients with too many missings values     
    NrMissingsPerSample=apply(MET_Data,2,function(x) sum(is.na(x)))/nrow(MET_Data)
    cat("Removing",sum(NrMissingsPerSample>MissingValueThreshold),"patients with more than",MissingValueThreshold*100,"% missing values.\n")
    if (sum(NrMissingsPerSample>MissingValueThreshold)>0) MET_Data=MET_Data[,NrMissingsPerSample<MissingValueThreshold,drop=FALSE]
    
    # removing clones with too many missing values
    NrMissingsPerGene=apply(MET_Data,1,function(x) sum(is.na(x)))/ncol(MET_Data)
    cat("Removing",sum(NrMissingsPerGene>MissingValueThreshold),"genes with more than",MissingValueThreshold*100,"% missing values.\n")
    if (sum(NrMissingsPerGene>MissingValueThreshold)>0) MET_Data=MET_Data[NrMissingsPerGene<MissingValueThreshold,,drop=FALSE]
    
    # knn impute using Tibshirani's method
    if (length(colnames(MET_Data))>1) {
        k=15
        KNNresults=impute::impute.knn(as.matrix(MET_Data),k)
        MET_Data_KNN=KNNresults$data
        
        # cleaning up sample names     
        return(MET_Data_KNN)
        
    } else {
        # when only 1 sample,need to make a matrix again
        #MET_Data=as.matrix(MET_Data)
        
        return(MET_Data)    
    }     
}

#' The TCGA_GENERIC_CheckBatchEffect function
#' 
#' Internal. Checks if batch correction is needed.
#' @param GEN_Data matrix with data to be corrected for batch effects.
#' @param BatchData Batch data.
#' @return list with results.
#' @keywords internal
#'
TCGA_GENERIC_CheckBatchEffect <-function(GEN_Data, BatchData) {
    
    # first match the samples to the batch
    Order=match(colnames(GEN_Data),BatchData[,1])
    BatchDataSelected=BatchData[Order,]
    BatchDataSelected$Batch <- factor(BatchDataSelected$Batch) 
    
    # PCA analysis
    # alternatively use fast.prcomp from package gmodels, but tests do not show this is faster
    PCAanalysis=prcomp(t(GEN_Data))
    PCdata=PCAanalysis$x
    #plot(PCdata[,1]~BatchDataSelected[,3])
    
    if (length(unique(BatchDataSelected$Batch[!is.na(BatchDataSelected$Batch)]))>1) {
        tmp=aov(PCdata[,1]~BatchDataSelected[,3])          
        return(list(Pvalues=summary(tmp),PCA=PCdata,BatchDataSelected=BatchDataSelected))
    } else {
        return(-1)
    }
}

#' The TCGA_BatchCorrection_MolecularData function
#' 
#' Internal. Wrapper to perform batch correction.
#' @param GEN_Data matrix with data to be corrected for batch effects.
#' @param BatchData Batch data.
#' @param MinInBatch minimum number of samples per batch.
#' @return The corrected data matrix.
#' @keywords internal
#'
TCGA_BatchCorrection_MolecularData <- function (GEN_Data,BatchData,MinInBatch) {
    
    # Remove samples with batch number 0
    if (length(-which(BatchData[,3]==0))>0) {
        BatchData=BatchData[-which(BatchData[,3]==0),]
    }
    
    # remove batches that are too small     
    #MinInBatch=5
    PresentSamples=is.element(BatchData[,1],colnames(GEN_Data))
    
    # changes this April 2014, such that the BatchDataSelected only deals with samples in the current GEN_Data
    BatchDataSelected=BatchData[PresentSamples,] 
    if (sum(PresentSamples) != length(colnames(GEN_Data))) BatchDataSelected=BatchData[-which(PresentSamples==FALSE),]
    BatchDataSelected$Batch <- factor(BatchDataSelected$Batch)
    
    NrPerBatch=table(BatchDataSelected$Batch)
    SmallBatches=NrPerBatch<MinInBatch
    BatchesToBeRemoved=names(SmallBatches)[which(SmallBatches==TRUE)]
    SamplesToBeRemoved=as.character(BatchDataSelected[which(BatchDataSelected$Batch %in% BatchesToBeRemoved),1])
    
    if (length(colnames(GEN_Data))-length(which(colnames(GEN_Data) %in% SamplesToBeRemoved))>5) { # just checking if we have enough samples after removing the too small batches
        if (length(which(colnames(GEN_Data) %in% SamplesToBeRemoved))>0) {
            cat("Removing",length(which(colnames(GEN_Data) %in% SamplesToBeRemoved)),"samples because their batches are too small.\n")
            GEN_Data=GEN_Data[,-which(colnames(GEN_Data) %in% SamplesToBeRemoved)]
        }          
        # batch correction with Combat, incorporate check for only 1 batch
        BatchCheck=TCGA_GENERIC_CheckBatchEffect(GEN_Data,BatchData)
        
        if (is.list(BatchCheck)) {
            GEN_Data_Corrected=TCGA_GENERIC_BatchCorrection(GEN_Data,BatchData)
            BatchCheck=TCGA_GENERIC_CheckBatchEffect(GEN_Data_Corrected,BatchData)
            return(GEN_Data_Corrected)
        } else {
            cat("Only one batch, no batch correction possible.\n")
            return(GEN_Data)
        }
    } else {
        cat("The nr of samples becomes to small, no batch correction possible.\n")
        return(GEN_Data)
    }
}

#' The TCGA_GENERIC_BatchCorrection function
#' 
#' Internal. Performs batch correction.
#' @param GEN_Data matrix with data to be corrected for batch effects.
#' @param BatchData Batch data.
#' @return The corrected data matrix.
#' @keywords internal
#'
TCGA_GENERIC_BatchCorrection <-function(GEN_Data,BatchData) {
    
    # select only samples with batch, others get deleted
    WithBatchSamples=is.element(colnames(GEN_Data),BatchData[,1])
    if (length(which(WithBatchSamples==FALSE))>0) GEN_Data=GEN_Data[,-which(WithBatchSamples==FALSE)]
    
    # select only the batch data that is present in the current data set, remove others (remember, the batch data is for all of TCGA)
    PresentSamples=is.element(BatchData[,1],colnames(GEN_Data))
    BatchDataSelected=BatchData
    if (sum(PresentSamples) != length(colnames(GEN_Data))) BatchDataSelected=BatchData[-which(PresentSamples==FALSE),]
    BatchDataSelected$Batch <- factor(BatchDataSelected$Batch)
    BatchDataSelected$ArrayName <- factor(BatchDataSelected$ArrayName)
    
    # reordening samples (not really necessary as Combat does this too)
    order <- match(colnames(GEN_Data),BatchDataSelected[,1])
    BatchDataSelected=BatchDataSelected[order,]
    BatchDataSelected$Batch <- factor(BatchDataSelected$Batch)
    
    # running combat
    CombatResults=ComBat_NoFiles(GEN_Data,BatchDataSelected)
    
    GEN_Data_Corrected=CombatResults[,-1]
    class(GEN_Data_Corrected) <- "numeric"
    return(GEN_Data_Corrected)
}

#' The Download_GeneExpression function
#' 
#' Downloads gene expression data from TCGA.
#' @param CancerSite character of length 1 with TCGA cancer code.
#' @param TargetDirectory character with directory where a folder for downloaded files will be created.
#' @param downloadData logical indicating if data should be downloaded (default: TRUE). If false, the url of the desired data is returned.
#' @return list with paths to downloaded files for both 27k and 450k methylation data.
#' @details This function downloads RNAseq data (file tag "mRNAseq_Preprocess.Level_3"), with the exception for OV and GBM, for which micro array data is
#' downloaded since there is not enough RNAseq data
#' @export
#' @keywords download
#' @examples 
#' \dontrun{
#' 
#' # Optional register cluster to run in parallel
#' library(doParallel)
#' cl <- makeCluster(5)
#' registerDoParallel(cl)
#' 
#' # Gene expression data for ovarian cancer
#' cancerSite <- "OV"
#' targetDirectory <- paste0(getwd(), "/")
#' 
#' # Downloading gene expression data
#' GEdirectories <- Download_GeneExpression(cancerSite, targetDirectory, TRUE)
#' 
#' # Processing gene expression data
#' GEProcessedData <- Preprocess_GeneExpression(cancerSite, GEdirectories)
#' 
#' # Saving gene expression processed data
#' saveRDS(GEProcessedData, file = paste0(targetDirectory, "GE_", cancerSite, "_Processed.rds"))
#' 
#' stopCluster(cl)
#' }
#' 
Download_GeneExpression <- function(CancerSite,TargetDirectory,downloadData=TRUE) {    
    
    dir.create(TargetDirectory,showWarnings=FALSE)
    
    # Settings
    TCGA_acronym_uppercase=toupper(CancerSite)
    
    # get RNA seq data (GBM does not have much RNAseq data.)
    dataType='stddata'	
    dataFileTag='mRNAseq_Preprocess.Level_3'	 
    
    #special case for GBM and OV, not enough RNAseq data, so using the microarray data instead
    if (CancerSite=="GBM") { 	             
        dataFileTag=c('Merge_transcriptome__agilentg4502a_07_1__unc_edu__Level_3__unc_lowess_normalization_gene_level__data','Merge_transcriptome__agilentg4502a_07_2__unc_edu__Level_3__unc_lowess_normalization_gene_level__data')        	         
    } else if(CancerSite=="OV") {	               
        dataFileTag='Merge_transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data'        
    }	
    cat('Searching MA data for:',CancerSite,"\n")
    if (length(dataFileTag)==1) {	  
        MAdirectories=get_firehoseData(downloadData,saveDir=TargetDirectory,TCGA_acronym_uppercase=TCGA_acronym_uppercase,dataFileTag=dataFileTag)    	
    } else {	    # a few data sets have multiple gene expression data sets. 
        MAdirectories=c()	  
        for (i in 1:length(dataFileTag)) {
            MAdirectories=c(MAdirectories,get_firehoseData(downloadData,saveDir=TargetDirectory,TCGA_acronym_uppercase=TCGA_acronym_uppercase,dataFileTag=dataFileTag[i]))	 
        }        
    }
    
    return(MAdirectories=MAdirectories)
}

#' The Preprocess_GeneExpression function
#' 
#' Pre-processes gene expression data from TCGA.
#' @param CancerSite character of length 1 with TCGA cancer code.
#' @param MAdirectories character vector with directories with the downloaded data. It can be the object returned by the Download_DNAmethylation function.
#' @param MissingValueThresholdGene threshold for missing values per gene. Genes with a percentage of NAs greater than this threshold are removed. Default is 0.3.
#' @param MissingValueThresholdSample threshold for missing values per sample. Samples with a percentage of NAs greater than this threshold are removed. Default is 0.1.
#' @details 
#' Pre-process includes eliminating samples and genes with too many NAs, imputing NAs, and doing Batch correction.
#' @return List with the pre-processed data matrix for cancer and normal samples.
#' @export
#' @keywords preprocess
#' @examples 
#' \dontrun{
#' 
#' # Optional register cluster to run in parallel
#' library(doParallel)
#' cl <- makeCluster(5)
#' registerDoParallel(cl)
#' 
#' # Gene expression data for ovarian cancer
#' cancerSite <- "OV"
#' targetDirectory <- paste0(getwd(), "/")
#' 
#' # Downloading gene expression data
#' GEdirectories <- Download_GeneExpression(cancerSite, targetDirectory, TRUE)
#' 
#' # Processing gene expression data
#' GEProcessedData <- Preprocess_GeneExpression(cancerSite, GEdirectories)
#' 
#' # Saving gene expression processed data
#' saveRDS(GEProcessedData, file = paste0(targetDirectory, "GE_", cancerSite, "_Processed.rds"))
#' 
#' stopCluster(cl)
#' }
#'
Preprocess_GeneExpression <- function(CancerSite,MAdirectories,MissingValueThresholdGene=0.3,MissingValueThresholdSample=0.1) {    
    
    get("BatchData")
    MinPerBatchCancer=5    
    MinPerBatchNormal=2    
    
    # Processing MA data, special case for OV and GBM where no RNA seq data is available
    if (CancerSite=="OV" || CancerSite=="GBM") { 
        MAstring='transcriptome__agilent'
    } else if (CancerSite=="STAD" || CancerSite=="ESCA") { # for these cancers RSEM data does not exist. 
        MAstring='mRNAseq_RPKM_log2.txt'
    } else {
        MAstring='mRNAseq_RSEM_normalized_log2.txt'
    }
    
    if (grepl("Windows", Sys.info()['sysname'])) {
        # If Windows I'll create a virtual drive to handle the long file names issue
        # Create a virtual drive to overcome long names issue in Windows
        MAdirectoriesOrig <- MAdirectories
        virtualDir <- MAdirectories
        virtualDir <- gsub("\\\\", "/", virtualDir)
        virtualDir <- substr(virtualDir, 1, nchar(virtualDir) - 1)
        system(paste("subst x:", virtualDir))
        MAdirectories <- "x://"
    }
    
    #cat("Loading mRNA data.\n")
    if (length(MAdirectories)>1) {      		
        cat("\tFound multiple MA data sets.\n")
        DataSetsCancer=list()
        GeneListsCancer=list()
        SampleListsCancer=list()                
        MetaBatchDataCancer=data.frame()
        
        DataSetsNormal=list()
        GeneListsNormal=list()
        SampleListsNormal=list()                
        MetaBatchDataNormal=data.frame()
        for (i in 1:length(MAdirectories)) {        
            cat("\tProcessing data set",i,"\n")
            MAfiles=dir(MAdirectories[i])
            MatchedFile=grep(MAstring,MAfiles)        
            if (length(MatchedFile)>0) {        
                # Getting the cancer data first
                DataSetsCancer[[i]]=Preprocess_MAdata_Cancer(CancerSite,MAdirectories[i],MAfiles[MatchedFile],MissingValueThresholdGene,MissingValueThresholdSample)
                GeneListsCancer[[i]]=rownames(DataSetsCancer[[i]])
                SampleListsCancer[[i]]=colnames(DataSetsCancer[[i]])                    
                currentBatchCancer=matrix(i,length(colnames(DataSetsCancer[[i]])),1) # growing a batch data object
                currentBatchDataCancer=data.frame(ArrayName=colnames(DataSetsCancer[[i]]),SampleName=colnames(DataSetsCancer[[i]]),Batch=currentBatchCancer)
                MetaBatchDataCancer=rbind(MetaBatchDataCancer,currentBatchDataCancer)
                
                # Getting the normal data as well.
                DataSetsNormal[[i]]=Preprocess_MAdata_Normal(CancerSite,MAdirectories[i],MAfiles[MatchedFile],MissingValueThresholdGene,MissingValueThresholdSample)
                GeneListsNormal[[i]]=rownames(DataSetsNormal[[i]])
                SampleListsNormal[[i]]=colnames(DataSetsNormal[[i]])                    
                currentBatchNormal=matrix(i,length(colnames(DataSetsNormal[[i]])),1) # growing a batch data object
                currentBatchDataNormal=data.frame(ArrayName=colnames(DataSetsNormal[[i]]),SampleName=colnames(DataSetsNormal[[i]]),Batch=currentBatchNormal)
                MetaBatchDataNormal=rbind(MetaBatchDataNormal,currentBatchDataNormal)
                
            } else {
                cat("MA file not found for this cancer.\n")
            }           
        }
        # combine data sets with Combat. 
        cat("Combining data sets.\n")
        OverlapProbesCancer=Reduce(intersect,GeneListsCancer)
        OverlapProbesNormal=Reduce(intersect,GeneListsNormal)
        OverlapSamplesCancer=Reduce(intersect,SampleListsCancer)    
        OverlapSamplesNormal=Reduce(intersect,SampleListsNormal)    
        if (length(OverlapSamplesCancer)>0 | length(OverlapSamplesNormal)>0) {
            cat('This should not happen. There is overlap between cancer or normal samples. No solution yet.\n')           
        }        
        
        for (i in 1:length(MAdirectories)) {
            DataSetsCancer[[i]]=DataSetsCancer[[i]][OverlapProbesCancer,]
            DataSetsNormal[[i]]=DataSetsNormal[[i]][OverlapProbesNormal,]
        }
        # combat on cancer data sets. 
        MA_TCGA_Cancer=Reduce(cbind,DataSetsCancer)
        MA_TCGA_Cancer=TCGA_BatchCorrection_MolecularData(MA_TCGA_Cancer,MetaBatchDataCancer,MinPerBatchCancer)    
        
        # combat on normal data sets. 
        MA_TCGA_Normal=Reduce(cbind,DataSetsNormal)
        MA_TCGA_Normal=TCGA_BatchCorrection_MolecularData(MA_TCGA_Normal,MetaBatchDataNormal,MinPerBatchNormal)    
        
    } else {
        
        MAfiles=dir(MAdirectories)
        MatchedFile=grep(MAstring,MAfiles)        
        if (length(MatchedFile)>0) {                  
            MA_TCGA_Cancer=Preprocess_MAdata_Cancer(CancerSite,MAdirectories,MAfiles[MatchedFile],MissingValueThresholdGene,MissingValueThresholdSample)
            MA_TCGA_Normal=Preprocess_MAdata_Normal(CancerSite,MAdirectories,MAfiles[MatchedFile],MissingValueThresholdGene,MissingValueThresholdSample)               
            cat("There are",length(colnames(MA_TCGA_Cancer)),"cancer samples and",length(colnames(MA_TCGA_Normal)),"normal samples in",CancerSite,"\n")
        } else {               
            stop("MA file not found for this cancer.\n")
        }           
    }
    MA_TCGA_Cancer = TCGA_GENERIC_CleanUpSampleNames(MA_TCGA_Cancer, 12)
    if (ncol(MA_TCGA_Normal) == 0) {
        MA_TCGA_Normal = NULL
    } else {
        MA_TCGA_Normal = TCGA_GENERIC_CleanUpSampleNames(MA_TCGA_Normal, 12)
    }
    
    if (grepl("Windows", Sys.info()['sysname'])) system("subst x: /D") #stop virtual drive
    
    return(list(GEcancer = MA_TCGA_Cancer, GEnormal = MA_TCGA_Normal))
}

#' The Preprocess_MAdata_Cancer function
#' 
#' Internal. Pre-process gene expression data for cancer samples.
#' @param CancerSite TCGA code for the cancer site.
#' @param Directory Directory.
#' @param File File.
#' @param MissingValueThresholdGene threshold for missing values per gene. Genes with a percentage of NAs greater than this threshold are removed. Default is 0.3.
#' @param MissingValueThresholdSample threshold for missing values per sample. Samples with a percentage of NAs greater than this threshold are removed. Default is 0.1.
#' @return The data matrix.
#' @keywords internal
#'
Preprocess_MAdata_Cancer <- function(CancerSite,Directory,File,MissingValueThresholdGene=0.3,MissingValueThresholdSample=0.1) {    
    
    get("BatchData")
    
    MinPerBatch=5   
    cat("Loading cancer mRNA data.\n")
    cat("\tMissing value estimation.\n")
    MA_TCGA=TCGA_Load_MolecularData(paste(Directory,File,sep=''), MissingValueThresholdGene, MissingValueThresholdSample)        
    Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(MA_TCGA))        
    if (CancerSite =='LAML') {
        MA_TCGA=MA_TCGA[,Samplegroups$PeripheralBloodCancer,drop=F]
    } else {
        MA_TCGA=MA_TCGA[,Samplegroups$Primary,drop=F]
    }          
    cat("\tBatch correction.\n")
    MA_TCGA=TCGA_BatchCorrection_MolecularData(MA_TCGA,BatchData,MinPerBatch)
    
    cat("\tProcessing gene ids and merging.\n")
    Genes=rownames(MA_TCGA)
    SplitGenes=limma::strsplit2(Genes,'\\|')
    rownames(MA_TCGA)=SplitGenes[,1]        
    MA_TCGA=MA_TCGA[!rownames(MA_TCGA) %in% '?',,drop=F]        
    MA_TCGA=TCGA_GENERIC_MergeData(unique(rownames(MA_TCGA)),MA_TCGA)  
    
    return(MA_TCGA=MA_TCGA)
}

#' The Preprocess_MAdata_Normal function
#' 
#' Internal. Pre-process gene expression data for normal samples.
#' @param CancerSite TCGA code for the cancer site.
#' @param Directory Directory.
#' @param File File.
#' @param MissingValueThresholdGene threshold for missing values per gene. Genes with a percentage of NAs greater than this threshold are removed. Default is 0.3.
#' @param MissingValueThresholdSample threshold for missing values per sample. Samples with a percentage of NAs greater than this threshold are removed. Default is 0.1.
#' @return The data matrix.
#' @keywords internal
#'
Preprocess_MAdata_Normal <- function(CancerSite,Directory,File,MissingValueThresholdGene=0.3,MissingValueThresholdSample=0.1) {    
    
    get("BatchData")
    
    MinPerBatch=2 # less samples in one batch when dealing with normals.
    
    cat("Loading normal mRNA data.\n")
    cat("\tMissing value estimation.\n")
    MA_TCGA=TCGA_Load_MolecularData(paste(Directory,File,sep=''), MissingValueThresholdGene, MissingValueThresholdSample)        
    Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(MA_TCGA))        
    if (CancerSite =='LAML') {
        MA_TCGA=MA_TCGA[,Samplegroups$BloodNormal,drop=F]
    } else {
        # MA_TCGA=MA_TCGA[,Samplegroups$SolidNormal]
        MA_TCGA=MA_TCGA[,Samplegroups$SolidNormal,drop=F]
    }          
    cat("\tBatch correction.\n")
    MA_TCGA=TCGA_BatchCorrection_MolecularData(MA_TCGA,BatchData,MinPerBatch)
    
    cat("\tProcessing gene ids and merging.\n")
    Genes=rownames(MA_TCGA)
    SplitGenes=limma::strsplit2(Genes,'\\|')
    rownames(MA_TCGA)=SplitGenes[,1]        
    # MA_TCGA=MA_TCGA[!rownames(MA_TCGA) %in% '?',]        
    MA_TCGA=MA_TCGA[!rownames(MA_TCGA) %in% '?', , drop=F]
    MA_TCGA=TCGA_GENERIC_MergeData(unique(rownames(MA_TCGA)),MA_TCGA)  
    
    return(MA_TCGA=MA_TCGA)
}

#' The TCGA_Load_MolecularData function
#' 
#' Internal. Reads in gene expressiondata. Deletes samples and genes with more NAs than the respective thresholds. Imputes other NAs values.
#' @param Filename name of the file with the data.
#' @param MissingValueThresholdGene threshold for missing values per gene. Genes with a percentage of NAs greater than this threshold are removed. Default is 0.3.
#' @param MissingValueThresholdSample threshold for missing values per sample. Samples with a percentage of NAs greater than this threshold are removed. Default is 0.1.
#' @return gene expression data.
#' @keywords internal
#'
TCGA_Load_MolecularData <- function(Filename, MissingValueThresholdGene = 0.3, MissingValueThresholdSample = 0.1) {     
    
    # loading the data in R
    Matching=regexpr('.mat$',Filename,perl=TRUE)
    if (Matching[[1]]>0) {     
        # MET_Data=TCGA_GENERIC_ReadDataMatrixMatFile(Filename)
        DataList=R.matlab::readMat(Filename)
        MATdata=as.matrix(DataList$RawData)
        rownames(MATdata)=DataList[[2]]
        colnames(MATdata)=DataList[[3]]
        
    } else {
        MET_Data=read.csv(Filename,sep="\t",row.names=1,header=TRUE,na.strings=c("NA","null"))
    }
    if (rownames(MET_Data)[1]=='Composite Element REF') {
        cat("Removing first row with text stuff.\n")
        MET_Data=MET_Data[-1,]         
        Genes=rownames(MET_Data)
        MET_Data=apply(MET_Data,2,as.numeric)
        rownames(MET_Data)=Genes
    }
    
    SampleNames=colnames(MET_Data)
    SampleNames=gsub('\\.','-',SampleNames)
    colnames(MET_Data)=SampleNames
    
    # removing clones with too many missing values
    NrMissingsPerGene=apply(MET_Data,1,function(x) sum(is.na(x))) /ncol(MET_Data)
    cat("Removing",sum(NrMissingsPerGene>MissingValueThresholdGene),"genes with more than",MissingValueThresholdGene*100,"% missing values.\n")
    if (sum(NrMissingsPerGene>MissingValueThresholdGene)>0) MET_Data=MET_Data[NrMissingsPerGene<MissingValueThresholdGene,]
    
    # removing patients with too many missings values     
    NrMissingsPerSample=apply(MET_Data,2,function(x) sum(is.na(x))) /nrow(MET_Data)
    cat("Removing",sum(NrMissingsPerSample>MissingValueThresholdSample),"patients with more than",MissingValueThresholdSample*100,"% missing values.\n")
    if (sum(NrMissingsPerSample>MissingValueThresholdSample)>0) MET_Data=MET_Data[,NrMissingsPerSample<MissingValueThresholdSample]
    
    # knn impute using Tibshirani's method     
    if (length(colnames(MET_Data))>1) {
        k=15
        KNNresults=impute::impute.knn(as.matrix(MET_Data),k)
        MET_Data_KNN=KNNresults$data
        
        # cleaning up sample names
        MET_Data_KNN_Clean=TCGA_GENERIC_CleanUpSampleNames(MET_Data_KNN,15)
        return(MET_Data_KNN_Clean)
        
    } else {
        # when only 1 sample,need to make a matrix again
        #MET_Data=as.matrix(MET_Data)
        MET_Data_Clean=TCGA_GENERIC_CleanUpSampleNames(MET_Data,15)
        return(MET_Data_Clean)    
    }
    
}

#' The TCGA_GENERIC_MergeData function
#' 
#' Internal.
#' @param NewIDListUnique unique rownames of data.
#' @param DataMatrix data matrix.
#' @return data matrix.
#' @keywords internal
#'
TCGA_GENERIC_MergeData <-function(NewIDListUnique, DataMatrix) {
    
    NrUniqueGenes=length(NewIDListUnique)
    MergedData=matrix(0,NrUniqueGenes,length(colnames(DataMatrix)))
    for (i in 1:NrUniqueGenes) {
        currentID=NewIDListUnique[i]          
        # tmpData=DataMatrix[which(rownames(DataMatrix) %in% currentID),]
        tmpData=DataMatrix[which(rownames(DataMatrix) %in% currentID), , drop=F]
        if (length(rownames(tmpData)) >1) {
            MergedData[i,]=colMeans(tmpData)
        } else {
            MergedData[i,]=tmpData
        }
    }
    rownames(MergedData)=NewIDListUnique
    colnames(MergedData)=colnames(DataMatrix)
    
    return(MergedData)
}

# #' The TCGA_GENERIC_ReadDataMatrixMatFile function
# #' 
# #' Internal. Reads a MAT file structure from a connection or a file. 
# #' @param Filename name of the file with the data.
# #' @return matrix with the data.
# #' @keywords internal
# #'
# TCGA_GENERIC_ReadDataMatrixMatFile <- function(Filename) {
# 
#     DataList=R.matlab::readMat(Filename)
#     
#     MATdata=as.matrix(DataList$RawData)
#     rownames(MATdata)=DataList[[2]]
#     colnames(MATdata)=DataList[[3]]
#     
#     return(MATdata)
# }

#' The ClusterProbes function
#' 
#' This function uses the annotation for Illumina methylation arrays to map each probe to a gene. Then, for each gene,
#' it clusters all its CpG sites using hierchical clustering and Pearson correlation as distance and complete linkage. 
#' If data for normal samples is provided, only overlapping probes between cancer and normal samples are used. 
#' Probes with SNPs are removed. 
#' This function is prepared to run in parallel if the user registers a parallel structure, otherwise it runs sequentially.
#' This function also cleans up the sample names, converting them to the 12 digit format.
#' @param MET_Cancer data matrix for cancer samples.
#' @param MET_Normal data matrix for normal samples.
#' @param CorThreshold correlation threshold for cutting the clusters.
#' @return List with the clustered data sets and the mapping between probes and genes.
#' @export
#' @keywords cluter_probes
#' @importFrom foreach %dopar%
ClusterProbes <- function(MET_Cancer, MET_Normal, CorThreshold = 0.4) {
    
    # Top level function that implements an equivalent cluster algorithm but using hierarchical clustering with complete linkage. 
    
    # overlapping cancer & normal probes
    if (!is.null(MET_Normal)) {
        OverlapProbes=intersect(rownames(MET_Cancer),rownames(MET_Normal))
        MET_Cancer=MET_Cancer[OverlapProbes,]
        MET_Normal=MET_Normal[OverlapProbes,]
    }
    
    #Get probe information
    get("ProbeAnnotation")
    
    # remove probes with SNPs
    get("SNPprobes")
    
    GoodProbes=setdiff(rownames(MET_Cancer),SNPprobes)
    NrProbesToRemove=length(rownames(MET_Cancer))-length(GoodProbes)
    cat("Removing",NrProbesToRemove,"probes with SNPs.\n")
    MET_Cancer=MET_Cancer[GoodProbes,]
    if (!is.null(MET_Normal)) MET_Normal=MET_Normal[GoodProbes,]
    
    ###### only iterating over genes that have probes present
    # Getting the positions relative to probe annotation of the probes present in this data set. 
    PresentProbes=match(ProbeAnnotation[,1],rownames(MET_Cancer))
    UniqueGenes=sort(unique(ProbeAnnotation[!is.na(PresentProbes),2]))
    UniqueGenes=UniqueGenes[which(UniqueGenes != "")] # there is one empty one
    
    # create large matrix, delete zeros in the end.      
    MET_Cancer_C=matrix(0,length(rownames(MET_Cancer)),length(colnames(MET_Cancer)))
    colnames(MET_Cancer_C)=colnames(MET_Cancer)
    if (!is.null(MET_Normal)) {
        MET_Normal_C=matrix(0,length(rownames(MET_Cancer)),length(colnames(MET_Normal)))
        colnames(MET_Normal_C)=colnames(MET_Normal)
    }
    ProbeMapping=matrix(0,length(rownames(MET_Cancer)),2)
    ProbeCounter=1
    METmatrixCounter=1
    ClusteredRownames=c()
    
    # alternative to do it in parallel
    cat("Clustering",length(rownames(MET_Cancer)),"probes in CpG site clusters.\n")
    # cluster needs to know the function because it is defined only here, so it does not know other functions.
    tmpClusterResults = foreach::foreach(i=1:length(UniqueGenes), .export='TCGA_GENERIC_MET_ClusterProbes_Helper_ClusterGenes_with_hclust') %dopar% {
        TCGA_GENERIC_MET_ClusterProbes_Helper_ClusterGenes_with_hclust(UniqueGenes[i],ProbeAnnotation,MET_Cancer,MET_Normal,CorThreshold)
    }
    
    for ( i in 1:length(UniqueGenes) ) {            
        MET_Cancer_C[METmatrixCounter:(METmatrixCounter-1+nrow(tmpClusterResults[[i]][[1]])),]=tmpClusterResults[[i]][[1]]               
        if (!is.null(MET_Normal)) MET_Normal_C[METmatrixCounter:(METmatrixCounter-1+nrow(tmpClusterResults[[i]][[2]])),]=tmpClusterResults[[i]][[2]]
        METmatrixCounter=METmatrixCounter+nrow(tmpClusterResults[[i]][[1]])
        ClusteredRownames=c(ClusteredRownames,rownames(tmpClusterResults[[i]][[1]]))
        ProbeMapping[ProbeCounter:(ProbeCounter-1+nrow(tmpClusterResults[[i]][[3]])),]=tmpClusterResults[[i]][[3]]
        ProbeCounter=ProbeCounter+nrow(tmpClusterResults[[i]][[3]])
    }
    # remove excessively large matrix
    MET_Cancer_C=MET_Cancer_C[1:length(ClusteredRownames),]
    rownames(MET_Cancer_C)=ClusteredRownames
    MET_Cancer_C = TCGA_GENERIC_CleanUpSampleNames(MET_Cancer_C, 12)
    if (!is.null(MET_Normal)) {
        MET_Normal_C=MET_Normal_C[1:length(ClusteredRownames),]
        rownames(MET_Normal_C)=ClusteredRownames
        MET_Normal_C = TCGA_GENERIC_CleanUpSampleNames(MET_Normal_C, 12)
    } else {
        MET_Normal_C = NULL
    }
    
    cat("\nFound",length(rownames(MET_Cancer_C)),"CpG site clusters.\n")
    return(list(MET_Cancer_Clustered=MET_Cancer_C,MET_Normal_Clustered=MET_Normal_C,ProbeMapping=ProbeMapping))
}

#' The TCGA_GENERIC_MET_ClusterProbes_Helper_ClusterGenes_with_hclust function
#' 
#' Internal. Cluster probes into genes.
#' @param Gene gene.
#' @param ProbeAnnotation data set matching probes to genes.
#' @param MET_Cancer data matrix for cancer samples.
#' @param MET_Normal data matrix for normal samples.
#' @param CorThreshold correlation threshold for cutting the clusters.
#' @return List with the clustered data sets and the mapping between probes and genes.
#' @keywords internal
#'
TCGA_GENERIC_MET_ClusterProbes_Helper_ClusterGenes_with_hclust <- function(Gene,ProbeAnnotation,MET_Cancer,MET_Normal=NULL,CorThreshold=0.4) {
    
    # first lookup the probes matching a single gene. DO NOT USE grep, it does not do exact matching, but looks for the pattern anywhere !!!
    Probes=ProbeAnnotation[which(ProbeAnnotation[,2] == Gene),1]
    Probes=Probes[which(Probes %in% rownames(MET_Cancer))]
    
    METcancer_Clustered=matrix(0,0,length(colnames(MET_Cancer)))
    if (!is.null(MET_Normal)) METnormal_Clustered=matrix(0,0,length(colnames(MET_Normal)))
    Clusternames=c()
    InverseCorrelationThreshold=1-CorThreshold
    GeneClustersForProbeMapping=array(dim=length(Probes))
    if (length(Probes)>1) {               
        ProbeCorrelation=cor(t(MET_Cancer[Probes,]),method='pearson')
        ClusterResults=hclust(as.dist(1-ProbeCorrelation), method = "complete", members = NULL)
        #plot(ClusterResults)
        Clusters=cutree(ClusterResults,h=InverseCorrelationThreshold)
        for (i in 1:length(unique(Clusters)) ) {
            tmpGeneProbes=Probes[Clusters==i]
            if ( length(tmpGeneProbes)>1 ) {
                
                tmpAveragedProfile=colMeans(MET_Cancer[tmpGeneProbes,])
                METcancer_Clustered=rbind(METcancer_Clustered,tmpAveragedProfile)          
                # Same for normal
                if (!is.null(MET_Normal)) tmpAveragedProfile=colMeans(MET_Normal[tmpGeneProbes,])
                if (!is.null(MET_Normal)) METnormal_Clustered=rbind(METnormal_Clustered,tmpAveragedProfile)
            } else {          
                METcancer_Clustered=rbind(METcancer_Clustered,MET_Cancer[tmpGeneProbes,])
                if (!is.null(MET_Normal)) METnormal_Clustered=rbind(METnormal_Clustered,MET_Normal[tmpGeneProbes,])
            }
            Clusternames=c(Clusternames,paste(Gene,'---Cluster',i,sep=""))               
            pos=which(Probes %in% tmpGeneProbes)
            GeneClustersForProbeMapping[pos]=paste(Gene,"---Cluster",i,sep="")
        }
        rownames(METcancer_Clustered)=Clusternames
        if (!is.null(MET_Normal)) rownames(METnormal_Clustered)=Clusternames         
        
    } else {
        METcancer_Clustered=MET_Cancer[Probes,,drop=FALSE]
        if (!is.null(MET_Normal)) METnormal_Clustered=MET_Normal[Probes,,drop=FALSE]
        Clusternames=Gene
        GeneClustersForProbeMapping=Gene
        rownames(METcancer_Clustered)=Clusternames
        if (!is.null(MET_Normal)) rownames(METnormal_Clustered)=Clusternames
    }
    
    ProbeMapping=t(rbind(Probes,GeneClustersForProbeMapping))
    
    if (is.null(MET_Normal)) METnormal_Clustered = NULL
    return(list(METcancer_Clustered,METnormal_Clustered,ProbeMapping))
}

