########################################################################
# Modifications done on MethylMix code
########################################################################
#
# Added arguments in main function:
#   - lisfOfGenes: to specify for which genes run MethylMix. Can be useful when we want to look at only one gene or when we are doing analysis like the pancancer.
#   - filter: logical, default TRUE, perform the simple linear regression to pick negative correlated genes. Useful when doing analysis like the pancancer
# In function MethylMix_ModelGeneExpression: 
#   - Removed argument "Method" because wasn't being used.
#   - Changed the way CovariateData was being used so that now it is considered a factor for the model
#   - Now it's parallel
# In function MethylMix_ModelSingleGene:
#   - No longer use chi-square value to compare BIC values.
#   - Added arguments which we could now play with, but they are not available for the user, they are fixed in the main function in their traditional values.
#     Those arguments are: maxComp = 3, minSamplesPerGroup = -1,PvalueThreshold = 0.01, MeanDifferenceTreshold = 0.10
#   - minSamplesPerGroup: minimum number of samples in each component, to prevent creating mixture states with only 1 or 2 dots.
#     Default is 0 (same as in first package), but could be used to require that all the 
#     components have at least a certain number of samples. If you provide the value -1, we require at least 5% of the samples in a new mixture state to accept it.
# Parallelization:
#   - Removed all the building of parallel structures from the code as it is suggested (should be done by the user according to what's best for them)
#   - Unified the code in such a way that will run under parallelization or not, and will detect it automatically
#   - Made MethylMix_ModelGeneExpression run in parallel as well.
# Number of components to be evaluated:
#   - Can be changed with the maxComp argument internally, it's fixed for the user (3 as always).
#   - Changed the code so it's run in a for loop which doesn't depend on the number of components as it did before
# Documentation:
#   - Wrote roxygen2 documentation for each function to it's ready to use for building the package.
# In MethylMix_RemoveFlipOver:
#   - Before only took care of the case with 2 genes, now I added code to take care of a common pattern when there are 3 components.
# Plots:
#   - Moved to ggplot2 plots.
# Getting data from firehose:
#   - Dates handling modified
#   - Managed issue with long file names in Windows
# Predictions:
#   - Added function to predict mixture component for new data

# Fix for NOTEs from R CMD check (no visible binding for global variable)
if(getRversion() >= "2.15.1")  {
    utils::globalVariables(c('SNPprobes', 'ProbeAnnotation', 'BatchData'))
}

#' @importFrom grDevices dev.off pdf
NULL

#' @importFrom graphics lines par plot title
NULL

#' @importFrom stats anova aov as.dist cor cutree dbeta density dnorm hclust lm optim prcomp qqline qqnorm qqplot quantile rgamma t.test var wilcox.test
NULL

#' @importFrom utils download.file read.csv tail untar write.table
NULL

#' MethylMix: Mixture model for DNA methylation data in cancer.
#' 
#' MethylMix identifies DNA methylation driven genes by modeling DNA methylation data in cancer vs. normal and looking 
#' for homogeneous subpopulations. In addition matched gene expression data (e.g. from microarray technology or RNA sequencing) 
#' is used to identify functional DNA methylation events by requiring a negative correlation between methylation 
#' and gene expression of a particular gene. See references below.
#' @param METcancer	Matrix with the methylation data of cancer tissue with genes in rows and samples in columns.
#' @param GEcancer Gene expression data for cancer tissue with genes in rows and samples in columns.
#' @param METnormal	Matrix with the normal methylation data of the same genes as in METcancer. Again genes in rows and samples in columns. The samples do not have to match with the cancer data. If this argument is NULL, MethylMix will run without comparing to normal samples.
#' @param listOfGenes Vector with genes names to be evaluated, names must coincide with the names of the rows of METcancer.
#' @param filter Logical indicating if the linear regression to select genes with significative linear negative relation between methylation and gene expression should be performed (default: TRUE).
#' @param NoNormalMode Logical indicating if the methylation states found in the cancer samples should be compared to the normal samples (default: FALSE).
#' @param OutputRoot Path to store the MethylMix results object.
#' @return MethylMixResults is a list with the following components:
#' \item{MethylationDrivers}{Genes identified as transcriptionally predictive and differentially methylated by MethylMix.}
#' \item{NrComponents}{The number of methylation states found for each driver gene.}
#' \item{MixtureStates}{A list with the DM-values for each driver gene. 
#' Differential Methylation values (DM-values) are defined as the difference between 
#' the methylation mean in one mixture component of cancer samples and the methylation mean
#' in the normal samples, for a given gene.}
#' \item{MethylationStates}{Matrix with DM-values for all driver genes (rows) and all samples (columns).}
#' \item{Classifications}{Matrix with integers indicating to which mixture component each cancer sample was assigned to, for each gene.}
#' \item{Models}{Beta mixture model parameters for each driver gene.}
#' @export
#' @references 
#' Gevaert 0. \href{https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btv020}{MethylMix: an R package for identifying DNA methylation-driven genes}. Bioinformatics (Oxford, England). 2015;31(11):1839-41. doi:10.1093/bioinformatics/btv020.
#' 
#' Gevaert O, Tibshirani R, Plevritis SK. \href{http://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0579-8}{Pancancer analysis of DNA methylation-driven genes using MethylMix}. Genome Biology. 2015;16(1):17. doi:10.1186/s13059-014-0579-8.
#' @examples
#' # load the three data sets needed for MethylMix
#' data(METcancer)
#' data(METnormal)
#' data(GEcancer)
#' 
#' # run MethylMix on a small set of example data
#' MethylMixResults <- MethylMix(METcancer, GEcancer, METnormal)
#' 
#' \dontrun{
#' # run in parallel
#' library(doParallel)
#' cl <- makeCluster(5)
#' registerDoParallel(cl)
#' MethylMixResults <- MethylMix(METcancer, GEcancer, METnormal)
#' stopCluster(cl)
#' }
#' 
MethylMix <- function(METcancer, 
                      GEcancer,
                      METnormal = NULL,
                      listOfGenes = NULL, 
                      filter = TRUE,
                      NoNormalMode = FALSE, 
                      OutputRoot='') { 
    
    if (missing(METcancer)) stop("Need to provide METcancer matrix")
    if (missing(GEcancer)) stop("Need to provide GEcancer matrix")
    stopifnot(
        class(METcancer) == "matrix",
        class(GEcancer) == "matrix",
        class(METnormal) %in% c("matrix", "NULL"),
        class(listOfGenes) %in% c("character", "NULL"),
        class(filter) == "logical",
        class(NoNormalMode) == "logical",
        class(OutputRoot) == "character"
    )
    
    # Keep only the genes provided by user
    if (!is.null(listOfGenes)) {
        listOfGenes <- intersect(listOfGenes, rownames(METcancer))
        METcancer <- METcancer[listOfGenes, , drop = F]
    }
    
    ### Step 1: modeling the gene expression using cancer methylation data (beta values scale)
    if (filter) {
        FunctionalGenes <- MethylMix_ModelGeneExpression(METcancer, GEcancer, CovariateData = NULL)
    } else {
        # No regression filter, but keep genes present in both METcancer and GEcancer data, as the ones returned by the regression filter
        METcancer.split.names  <- sapply(strsplit(rownames(METcancer),  '---'), function(x) x[1])
        genes.to.keep.MET <- METcancer.split.names %in% rownames(GEcancer)
        FunctionalGenes <- rownames(METcancer)[genes.to.keep.MET]
    }
    
    ### Step 2: modeling the methylation data as a mixture of beta distributions 
    if (length(FunctionalGenes) > 0) {
        MixtureModelResults <- MethylMix_MixtureModel(METcancer, METnormal, FunctionalGenes, NoNormalMode)
    } else {
        cat("No transcriptionally predictive genes found.\n")
    }
    
    ### Step 3: write output to file
    if (OutputRoot != "") {
        saveRDS(MixtureModelResults, file = paste0(OutputRoot, "MethylMix_Results.rds"))
    }

    return(MixtureModelResults)
}

#' The MethylMix_ModelGeneExpression function
#' 
#' Model gene expression as a function of gene expression with a simple linear regression model. 
#' Genes with a significant negative linear association between DNA methylation and gene expression are returned.
#' @param METcancer matrix with methylation data for cancer samples (genes in rows, samples in columns).
#' @param GEcancer matrix with gene expression data for cancer samples (genes in rows, samples in columns).
#' @param CovariateData vector (numeric or character) indicating a covariate to be included in the model to adjust for it. Not used in an standard run of MethylMix.
#' It can be used if samples can from different tissue type, for example.
#' @return vector with the names of the genes for which there is a significant linear and negative association between methylation and gene expression.
#' @importFrom foreach %dopar%
#' @export
#' @examples 
#' # load data sets
#' data(METcancer)
#' data(GEcancer)
#' 
#' # model gene expression
#' MethylMixResults <- MethylMix_ModelGeneExpression(METcancer, GEcancer)
#' 
MethylMix_ModelGeneExpression <- function(METcancer, GEcancer, CovariateData = NULL) {      
    
    # overlapping samples     
    OverlapSamples = intersect(colnames(METcancer), colnames(GEcancer))
    cat("Found", length(OverlapSamples), "samples with both methylation and expression data.\n")
    GEcancer = GEcancer[, OverlapSamples, drop = FALSE]
    METcancer = METcancer[, OverlapSamples, drop = FALSE]
    if (!is.null(CovariateData)) CovariateData = CovariateData[OverlapSamples, , drop = FALSE]
    
    Rsquares = matrix(0, nrow = length(rownames(METcancer)), ncol = 1)    
    Genes = rownames(METcancer)  
    PvalueThreshold = 0.001
    RsquareThreshold = 0.1
    
    cat("Correlating methylation data with gene expression...\n")
    i <- NULL # to avoid "no visible binding for global variable" in R CMD check
    Rsquares = foreach::foreach(i = 1 : length(rownames(METcancer)), .combine = 'c') %dopar% {
        Rsq = 0
        tmpGene = unlist(strsplit(Genes[i], '---'))[1]
        pos = which(rownames(GEcancer) == tmpGene)
        if (length(pos) > 0) {
            if (!is.null(CovariateData)) {
                res = lm(GEcancer[pos, ] ~ METcancer[Genes[i], ] + factor(CovariateData))
                res.summary = summary(res)
                an = anova(res)
                if (res$coefficients[2] < 0 & res.summary$coefficients[2, 4] < PvalueThreshold ) { 
                    # Note: to also request  methylation effect to be bigger than tissue add above: & an$`Pr(>F)`[1] < an$`Pr(>F)`[2], but tissue effect is always bigger so now I remove it
                    Rsq = res.summary$r.squared
                }
            } else {           
                res = lm(GEcancer[pos, ] ~ METcancer[Genes[i], ])
                res.summary = summary(res)
                if (res$coefficients[2] < 0 & res.summary$coefficients[2, 4] < PvalueThreshold) {
                    Rsq= res.summary$r.squared 
                }
            }
        }
        Rsq
    }
    Rsquares = matrix(Rsquares, ncol = 1)
    
    # Rsquare threshold
    FunctionalGenes = Genes[Rsquares > RsquareThreshold]
    cat("\nFound", length(FunctionalGenes), "transcriptionally predictive genes.\n")
    return(FunctionalGenes)
}

#' The MethylMix_MixtureModel function
#' 
#' Internal. Prepares all the structures to store the results and calls in a foreach loop a function that fits the mixture model in each gene.
#' @param METcancer matrix with methylation data for cancer samples (genes in rows, samples in columns).
#' @param METnormal matrix with methylation data for normal samples (genes in rows, samples in columns). If NULL no comparison to normal samples will be done.
#' @param FunctionalGenes vector with genes names to be considered for the mixture models.
#' @param NoNormalMode logical, if TRUE no comparison to normal samples is performed. Defaults to FALSE.
#' @return MethylationStates matrix of DM values, with driver genes in the rows and samples in the columns.
#' @return NrComponents matrix with the number of components identified for each driver gene.
#' @return Models list with the mixture model fitted for each driver gene.
#' @return MethylationDrivers character vector with the genes found by MethylMix as differentially methylated and transcriptionally predictive (driver genes).
#' @return MixtureStates a list with a matrix for each driver gene containing the DM values.
#' @return Classifications a vector indicating to which component each sample was assigned.
#' @importFrom foreach %dopar%
#' @keywords internal
#' 
MethylMix_MixtureModel <- function(METcancer, METnormal = NULL, FunctionalGenes, NoNormalMode = FALSE) {
    
    # overlap of samples
    if (!is.null(METnormal)) {
        GeneOverlap = intersect(intersect(rownames(METcancer), rownames(METnormal)), FunctionalGenes)
        METcancer = METcancer[GeneOverlap, ,drop=FALSE]
        METnormal = METnormal[GeneOverlap, ,drop=FALSE]
    } else {
        GeneOverlap = intersect(rownames(METcancer), FunctionalGenes)
        METcancer = METcancer[GeneOverlap, ,drop=FALSE]
    }

    cat("\nStarting Beta mixture modeling.\n")
    
    MethylationStates = matrix(0, length(rownames(METcancer)), length(colnames(METcancer)))    
    rownames(MethylationStates) = rownames(METcancer)
    colnames(MethylationStates) = colnames(METcancer)  

    cat("Running Beta mixture model on",length(rownames(METcancer)),"genes and on",length(colnames(METcancer)),"samples.\n")
    
    i <- NULL # to avoid "no visible binding for global variable" in R CMD check
    if (!is.null(METnormal)) {
        res <- foreach::foreach(i = seq(rownames(METcancer)), .combine = "combineForEachOutput", .export = c("MethylMix_ModelSingleGene", "MethylMix_RemoveFlipOver", "blc_2", "betaEst_2")) %dopar% {
            MethylMix_ModelSingleGene(rownames(METcancer)[i], METcancer[i,], METnormal[i,], NoNormalMode)
        }
    } else {
        res <- foreach::foreach(i = seq(rownames(METcancer)), .combine = "combineForEachOutput", .export = c("MethylMix_ModelSingleGene", "MethylMix_RemoveFlipOver", "blc_2", "betaEst_2")) %dopar% {
            MethylMix_ModelSingleGene(rownames(METcancer)[i], METcancer[i,], NULL, TRUE)
        }
    }
    
    if (is.null(dim(res$Classifications))) res$Classifications <- matrix(res$Classifications, nrow = 1)
    if (is.null(dim(res$MethylationStates))) res$Classifications <- matrix(res$Classifications, nrow = 1)
    rownames(res$MethylationStates) <- rownames(res$Classifications) <- rownames(METcancer)
    colnames(res$MethylationStates) <- colnames(res$Classifications) <- colnames(METcancer)

    # Removing the genes without any differential methylation. 
    if (!NoNormalMode) {  
        # If NoNormalMode == T, no comparison to normal is made, and we don't remove genes with methylation states equal to 0
        # (for example, in pancancer analysis, running MethylMix only with normal samples, we use NoNormalMode = TRUE, and genes with only one state equal to 0 are kept.)
        NonZeroPositions = rowSums(res$MethylationStates) != 0
        res$NrComponents = res$NrComponents[NonZeroPositions, drop=FALSE]
        res$MixtureStates = res$MixtureStates[NonZeroPositions, drop=FALSE]
        res$Models = res$Models[NonZeroPositions, drop=FALSE]
        res$MethylationStates = res$MethylationStates[NonZeroPositions, , drop=FALSE]
        res$Classifications = res$Classifications[NonZeroPositions, , drop=FALSE]
    }
    
    # Adding names and removing things that we don't want to return
    res$MethylationDrivers = rownames(res$MethylationStates)
    res$FlipOverStates <- res$FunctionalGenes <- NULL
    names(res$NrComponents) = rownames(res$MethylationStates)
    names(res$MixtureStates) = rownames(res$MethylationStates)
    names(res$Models) = rownames(res$MethylationStates)
    
    return(res)
}

#' The MethylMix_ModelSingleGene function
#' 
#' Internal. For a given gene, this function fits the mixture model, selects the number of components and defines the respective methylation states. 
#' @param GeneName character string with the name of the gene to model
#' @param METdataVector vector with methylation data for cancer samples.
#' @param METdataNormalVector vector with methylation data for normal samples. It can be NULL and then no normal mode will be used.
#' @param NoNormalMode logical, if TRUE no comparison to normal samples is performed. Defaults to FALSE.
#' @param maxComp maximum number of mixture components admitted in the model (3 by default).
#' @param minSamplesPerGroup minimum number of samples required to belong to a new mixture component in order to accept it. Defaul is 1 (not used). If -1, each component has to have at least 5\% of all cancer samples.
#' @param PvalueThreshold threshold to consider results significant.
#' @param MeanDifferenceTreshold threshold in beta value scale from which two methylation means are considered different.
#' @return NrComponents number of components identified.
#' @return Models an object with the parameters of the model fitted.
#' @return MethylationStates vector with DM values for each sample.
#' @return MixtureStates vector with DMvalues for each component.
#' @return Classifications a vector indicating to which component each sample was assigned.
#' @return FlipOverState FlipOverState
#' @details maxComp, PvalueThreshold, METDiffThreshold, minSamplesPerGroup are arguments for this function but are fixed in their default values for the user
#' because they are not available in the main MethylMix function, to keep it simple. It would be easy to make them available to the user if we want to.
#' @keywords internal
#' 
MethylMix_ModelSingleGene <-function(GeneName, METdataVector, METdataNormalVector = NULL, NoNormalMode = FALSE, maxComp = 3,
                                     PvalueThreshold = 0.01, MeanDifferenceTreshold = 0.10, minSamplesPerGroup = 1) {
    
    # 1 component model
    w0.m = matrix(1, length(METdataVector), 1)
    mods = vector("list", maxComp + 1)
    mods[[1]] = blc_2(matrix(METdataVector, ncol = 1), w = w0.m, maxiter = 100, tol = 1e-06, verbose = FALSE)
    bic = numeric(maxComp + 1)
    bic[1] = -2 * mods[[1]]$llike + 2 * log(length(METdataVector))  
    
    # 2 to maxComp components model
    for (comp in 2:maxComp) {
        
        # Divide initial groups using quantiles
        prob = seq(1/comp, (comp-1)/comp, 1/comp)
        tmpQuantiles = quantile(METdataVector, prob)
        w0.m = matrix(0, nrow = length(METdataVector), ncol = comp)
        tmpQuantiles = c(tmpQuantiles, Inf)
        w0.m[METdataVector < tmpQuantiles[1], 1] <- 1
        for (i in 2:comp) {
            w0.m[METdataVector >= tmpQuantiles[i - 1] & METdataVector < tmpQuantiles[i], i] <- 1
        }
        
        # Fit beta mixture model
        mods[[comp]] <- blc_2(matrix(METdataVector, ncol=1), w = w0.m, maxiter = 100, tol = 1e-06, verbose = FALSE)
        if (sum(is.na(mods[[comp]]$mu)) > 0) mods[[comp]]$llike=0
        df = comp * 3 - 1
        bic[[comp]] = -2 * mods[[comp]]$llike + df * log(length(METdataVector))
        
        # See differences between model's means to compare them to threshold
        model.means = sort(mods[[comp]]$mu)
        different.means = ifelse(all(abs(diff(model.means)) > MeanDifferenceTreshold), T, F)
        
        # Check if smallest group has at least minSamplesPerGroup observations:
        if (minSamplesPerGroup < 0) {
            minSamplesPerGroup = max(5, 0.05 * length(METdataVector))
        }
        minOK = ifelse(min(table(apply(mods[[comp]]$w, 1, which.max))) >= minSamplesPerGroup, TRUE, FALSE)
        
        # We try adding another component if the following 2 conditions are satisfied:
        #   A: Adding one component reduces BIC
        #   B: All absolute differences between methylation means in model with one extra component are above the MeanDifferenceThreshold
        # But I also check C = smallest group has at least minSamplesPerGroup observations:
        # So, continue with another component if A & B & C, else not continue, which is the same as saying not continue if !A OR !B OR !C  
        # If the model was improved try one more, else break
        if (bic[[comp]] >= bic[[comp - 1]] || !different.means || !minOK) {
            NrComponents = comp - 1
            break
        } else {
            # It improved, try one more (unless comp already is maxComp, in that case we end here)
            NrComponents = comp
        }
        
    }
    
    Model = mods[[NrComponents]]
    MethylationState = matrix(0, 1, length(METdataVector))
    FlipOverState = 0
    res = list(p.value = 1)
    MixtureStates = matrix(0, NrComponents, 1)
    classification = apply(Model$w, 1, which.max)
    
    if (NrComponents == 1) {
        if (!is.null(METdataNormalVector)) {
            res = wilcox.test(METdataVector, METdataNormalVector)
            Difference = mean(METdataVector) - mean(METdataNormalVector)
        } else {
            res = list(p.value = 1)
            Difference = mean(METdataVector)
        }
        if ((res$p.value < PvalueThreshold & abs(Difference) > MeanDifferenceTreshold) | NoNormalMode) {
            MethylationState[1, ] = Difference
            MixtureStates[1, 1] = Difference
        }
        cat(c(GeneName,": 1 component is best.\n"))
    } else {
        for (comp in 1:NrComponents) {
            METdataVector_comp = METdataVector[classification == comp]
            if (!is.null(METdataNormalVector)) {
                if (length(METdataVector_comp) > 0) res = wilcox.test(METdataVector_comp, METdataNormalVector) else res$p.value = 1
                Difference = mean(METdataVector_comp) - mean(METdataNormalVector)
            } else {
                res = list(p.value = 1)
                Difference = mean(METdataVector_comp)
            }
            if ((res$p.value < PvalueThreshold & abs(Difference) > MeanDifferenceTreshold) | NoNormalMode) {
                MethylationState[1, classification == comp] = Difference            
                MixtureStates[comp, 1] = Difference
            }  
        }
        # Flipover correction (in first MethylMix package was done only for 2 mixture states)
        if (NrComponents == 2 || NrComponents == 3) {
            OrigOrder = order(METdataVector)
            FlipOverResults = MethylMix_RemoveFlipOver(OrigOrder, MethylationState, classification, METdataVector, NrComponents)
            MethylationState = FlipOverResults$MethylationState
            classification = FlipOverResults$classification
            FlipOverState = FlipOverResults$LearnedState         
        }
        cat(c(GeneName,": ", NrComponents, " components are best.\n"))
    }                 
    return(list(MethylationStates = MethylationState,
                NrComponents = NrComponents,
                Models = list(Model),
                MixtureStates = list(MixtureStates),
                Classifications = classification,
                FlipOverStates = FlipOverState))
}

#' The combineForEachOutput function
#' 
#' Internal. Function to combine results from the foreach loop.
#' @param out1 result from one foreach loop.
#' @param out2 result from another foreach loop.
#' @return List with the combined results.
#' @keywords internal
#'
combineForEachOutput <- function(out1, out2) {
    NrComponents <- c(out1$NrComponents, out2$NrComponents)
    MethylationStates <- rbind(out1$MethylationStates, out2$MethylationStates)
    MixtureStates <- append(out1$MixtureStates, out2$MixtureStates)
    Models <- append(out1$Models, out2$Models)
    FlipOverStates <- c(out1$FlipOverStates, out2$FlipOverStates)
    Classifications <- rbind(out1$Classifications, out2$Classifications)
    return(list(NrComponents = NrComponents,
                MethylationStates = MethylationStates,
                MixtureStates = MixtureStates,
                Models = Models,
                FlipOverStates = FlipOverStates,
                Classifications = Classifications))
}

#' The MethylMix_RemoveFlipOver function
#' 
#' Internal. The estimated densities for each beta component can overlap, generating samples that look like being separated from their group.
#' This function re classifies such samples.
#' @param OrigOrder order of sorted values in the methylation vector.
#' @param MethylationState methylation states for this gene.
#' @param classification vector with integers indicating to wich component each sample was classified into.
#' @param METdataVector vector with methylation values from the cancer samples.
#' @param NrComponents number of components in this gene.
#' @param UseTrainedFlipOver .
#' @param FlipOverState .
#' @return Corrected vectors with methylation states and classification.
#' @keywords internal
#' 
MethylMix_RemoveFlipOver <- function(OrigOrder, MethylationState, classification, METdataVector, NrComponents, UseTrainedFlipOver = FALSE, FlipOverState = 0) {
    
    Differences = diff(MethylationState[1, OrigOrder]) 
    DifferencesSel = Differences[Differences != 0]  
    LearnedState = 0
    
    # If the model has 2 components, we perform exactly what was done originally in the first MethylMix package.
    # The new addition is to hadle the flipover when the model has 3 components.
    if (NrComponents == 2) {
        # This handles the escenario where one of the components has observations at both sides of the other
        # If (length(DifferencesSel) < 2, there's no mixing
        if (length(DifferencesSel) == 2) {
            if (DifferencesSel[1] * -1 == DifferencesSel[2]) {
                
                posDiff1 = which(Differences == DifferencesSel[1])
                stateSize1 = posDiff1
                posDiff2 = which(Differences == DifferencesSel[2])
                stateSize2 = length(Differences) - posDiff2
                
                if (UseTrainedFlipOver == TRUE) {
                    if (FlipOverState == 1) {
                        # I changed OrigOrder[posDiff2 + 1:length(Differences)] into OrigOrder[(posDiff2 + 1) : (length(Differences) + 1)] so that index is not out of range and NAs are not introduced (it didn't affect anything anyway)
                        MethylationState[1, OrigOrder[(posDiff2 + 1) : (length(Differences) + 1)]] = MethylationState[1, OrigOrder[posDiff2]]
                        classification[OrigOrder[(posDiff2 + 1) : (length(Differences) + 1)]] = classification[OrigOrder[posDiff2]]
                    } else if (FlipOverState == 2) {
                        MethylationState[1, OrigOrder[1:posDiff1]] = MethylationState[1, OrigOrder[posDiff1 + 1]]
                        classification[OrigOrder[1:posDiff1]] = classification[OrigOrder[posDiff1 + 1]]
                    }                    
                } else {
                    if (stateSize2 > stateSize1) {                    
                        MethylationState[1, OrigOrder[1:posDiff1]] = MethylationState[1, OrigOrder[posDiff1 + 1]]
                        classification[OrigOrder[1:posDiff1]] = classification[OrigOrder[posDiff1 + 1]]
                        LearnedState = 2
                    } else if (stateSize1 > stateSize2) {                    
                        MethylationState[1, OrigOrder[(posDiff2 + 1) : (length(Differences) + 1)]] = MethylationState[1, OrigOrder[posDiff2]]
                        classification[OrigOrder[(posDiff2 + 1) : (length(Differences) + 1)]] = classification[OrigOrder[posDiff2]]
                        LearnedState = 1
                    }  
                }
            }
        }
    } 
    
    # For models with 3 components, we have this new implementation that handles escenarios where only one of the
    # components is divided and separated between the other 2.
    if (NrComponents == 3) {
        # If length(DifferencesSel) == 2, there is no overlapping and mixing
        if (length(DifferencesSel) > 2) {
            
            # seq.of.states: Shows the order in which the different states are being mixing
            # posDiff: for each sequence of the same state, which is the position of the last element
            seq.of.states = numeric(0)
            posDiff = numeric(0)
            seq.of.class = numeric(0)
            for (i in 1:length(DifferencesSel)) {
                posDiff[i] = which(Differences == DifferencesSel[i])
                seq.of.states[i] = MethylationState[1, OrigOrder][posDiff[i]]
                seq.of.class[i] = classification[OrigOrder][posDiff[i]]
            }
            # Catch the last group:
            seq.of.states = c(seq.of.states, tail(MethylationState[1, OrigOrder], 1))
            seq.of.class= c(seq.of.class, tail(classification[OrigOrder], 1))
            posDiff = c(0, posDiff, length(MethylationState))
            
            # Size of each group:
            size = diff(posDiff)
            
            # A table of seq.of.states will show which state is separated
            tab = table(seq.of.states)
            
            # I will deal only with the case where only one state is separated
            if (sum(tab > 1) == 1) {
                
                # Mean of each subgroup
                means = tapply(METdataVector[OrigOrder], rep(1:length(size), size), mean)
                
                # Identify which is the state that is divided into subgroups
                separated.state = round(as.numeric(names(tab)[tab > 1]), 4)
                
                # Subgroups that correspond to this separated state
                subgr = which(round(seq.of.states, 4) == separated.state)
                
                # Subgroups well defined and separated
                subgr.ok = which(round(seq.of.states, 4) != separated.state)
                
                # In the separated state, the largest subgroup will remain, and the others
                # will be allocated to one of the ok subgroups, to the one of closest mean
                subgr.remains = subgr[which.max(size[subgr])]
                subgr.allocate = subgr[!subgr %in% subgr.remains]
                for (gr in subgr.allocate) {
                    allocate.in.group = subgr.ok[which.min(abs(means[gr] - means[subgr.ok]))]
                    pos.to.change = (posDiff[gr] + 1) : posDiff[gr + 1]
                    MethylationState[1, OrigOrder[pos.to.change]] = seq.of.states[allocate.in.group]
                    classification[OrigOrder[pos.to.change]] = seq.of.class[allocate.in.group]
                }
                LearnedState = 3 # I chose a 3 just to distiguish from the 1 or 2 that the original flipover function returns
            }
        }
    }
    return(list(MethylationState = MethylationState, classification = classification, LearnedState = LearnedState))
}

#' The blc_2 function
#' 
#' Internal. Fits a beta mixture model for any number of classes. Adapted from RPMM package.
#' @param Y Data matrix (n x j) on which to perform clustering.
#' @param w Initial weight matrix (n x k) representing classification.
#' @param maxiter Maximum number of EM iterations.
#' @param tol Convergence tolerance.
#' @param weights Case weights.
#' @param verbose Verbose output.
#' @return A list of parameters representing mixture model fit, including posterior weights and log-likelihood.
#' @keywords internal
#' 
blc_2 <- function (Y, w, maxiter = 25, tol = 1e-06, weights = NULL, verbose = TRUE) {
    Ymn <- min(Y[Y > 0], na.rm = TRUE)
    Ymx <- max(Y[Y < 1], na.rm = TRUE)
    Y <- pmax(Y, Ymn/2)
    Y <- pmin(Y, 1 - (1 - Ymx)/2) 
    Yobs = !is.na(Y)
    J <- dim(Y)[2] 
    K <- dim(w)[2]
    n <- dim(w)[1]
    if (n != dim(Y)[1]) 
        stop("Dimensions of w and Y do not agree")
    if (is.null(weights)) 
        weights <- rep(1, n)
    mu <- a <- b <- matrix(Inf, K, J)
    crit <- Inf
    for (i in 1:maxiter) {
        warn0 <- options()$warn
        options(warn = -1)
        eta <- apply(weights * w, 2, sum)/sum(weights)
        mu0 <- mu
        for (k in 1:K) {
            for (j in 1:J) {  
                ab <- betaEst_2(Y[, j], w[, k], weights)
                a[k, j] <- ab[1]
                b[k, j] <- ab[2]
                mu[k, j] <- ab[1]/sum(ab)
            }
        }     
        ww <- array(0, dim = c(n, J, K))
        for (k in 1:K) {
            for (j in 1:J) {
                ww[Yobs[, j], j, k] <- dbeta(Y[Yobs[, j], j], 
                                             a[k, j], b[k, j], log = TRUE)
            }
        }          
        options(warn = warn0)
        w <- apply(ww, c(1, 3), sum, na.rm = TRUE)
        wmax <- apply(w, 1, max)
        for (k in 1:K) w[, k] <- w[, k] - wmax
        w <- t(eta * t(exp(w)))
        like <- apply(w, 1, sum)
        w <- (1/like) * w
        llike <- weights * (log(like) + wmax)
        crit <- max(abs(mu - mu0))
        if (verbose) 
            print(crit)          
        if (is.na(crit) || crit < tol) 
            break
    }
    return(list(a = a, b = b, eta = eta, mu = mu, w = w, llike = sum(llike)))
}

#' The betaEst_2 function
#' 
#' Internal. Estimates a beta distribution via Maximum Likelihood. Adapted from RPMM package.
#' @param Y data vector.
#' @param w posterior weights.
#' @param weights Case weights.
#' @return (a,b) parameters.
#' @keywords internal
#' 
betaEst_2 <-function (Y, w, weights) {
    y=Y
    yobs = !is.na(y)
    if (sum(yobs) <= 1) 
        return(c(1, 1))
    y = y[yobs]
    w = w[yobs]
    weights = weights[yobs]
    N <- sum(weights * w)
    p <- sum(weights * w * y)/N
    v <- sum(weights * w * y * y)/N - p * p
    logab <- log(c(p, 1 - p)) + log(pmax(1e-06, p * (1 - p)/v - 
                                             1))
    if (sum(yobs) == 2) 
        return(exp(logab))
    opt <- try(optim(logab, RPMM::betaObjf, ydata = y, wdata = w, weights = weights, 
                     method = "BFGS"), silent = TRUE)
    if (inherits(opt, "try-error")) 
        return(c(1, 1)) 
    exp(opt$par) # if using optimx, exp(as.numeric(opt$par))
}

#' The MethylMix_PlotModel function.
#' 
#' Produces plots to represent MethylMix's output.
#' @param GeneName Name of the gene for which to create a MethylMix plot.
#' @param MixtureModelResults List returned by MethylMix function.
#' @param METcancer	Matrix with the methylation data of cancer tissue with genes in rows and samples in columns.
#' @param GEcancer Gene expression data for cancer tissue with genes in rows and samples in columns (optional).
#' @param METnormal	Matrix with the normal methylation data of the same genes as in METcancer (optional). Again genes in rows and samples in columns.
#' @return a list with MethylMix plots, a histogram of the methylation data (MixtureModelPlot) and a scatterplot between DNA methylation and gene expression
#' (CorrelationPlot, is NULL if gene expression data is not provided). Both plots show the different mixture components identified.
#' @export
#' @examples
#' # Load the three data sets needed for MethylMix
#' data(METcancer)
#' data(METnormal)
#' data(GEcancer)
#' 
#' # Run methylmix on a small set of example data
#' MethylMixResults <- MethylMix(METcancer, GEcancer, METnormal)
#' 
#' # Plot the most famous methylated gene for glioblastoma
#' MethylMix_PlotModel("MGMT", MethylMixResults, METcancer)
#' 
#' # Plot MGMT also with its normal methylation variation
#' MethylMix_PlotModel("MGMT", MethylMixResults, METcancer, METnormal = METnormal)
#' 
#' # Plot a MethylMix model for another gene
#' MethylMix_PlotModel("ZNF217", MethylMixResults, METcancer, METnormal = METnormal)
#' 
#' # Also plot the inverse correlation with gene expression (creates two separate plots)
#' MethylMix_PlotModel("MGMT", MethylMixResults, METcancer, GEcancer, METnormal)
#' 
#' # Plot all functional and differential genes
#' for (gene in MethylMixResults$MethylationDrivers) {
#'      MethylMix_PlotModel(gene, MethylMixResults, METcancer, METnormal = METnormal)
#' }
#' 
MethylMix_PlotModel <- function(GeneName, MixtureModelResults, METcancer, GEcancer = NULL, METnormal = NULL) {
    
    met <- x <- dens <- comp <- ge <- group <- ..density.. <- NULL # to avoid "no visible binding for global variable" in R CMD check
    
    GeneNameMA <- unlist(strsplit(GeneName, '---'))[1]
    
    if (GeneName %in% MixtureModelResults$MethylationDrivers && GeneName %in% rownames(METcancer)) {
        
        # Prepare data
        OverlapSamples <- colnames(METcancer)
        if (!is.null(GEcancer)) OverlapSamples <- intersect(colnames(METcancer), colnames(GEcancer))
        data <- data.frame(met = METcancer[GeneName, OverlapSamples, drop = F][1,],
                           group = factor(MixtureModelResults$Classification[GeneName, OverlapSamples], levels = 1:length(unique(MixtureModelResults$Classification[GeneName, OverlapSamples]))))
        if (!is.null(GEcancer)) data$ge <- GEcancer[GeneNameMA, OverlapSamples]
        data$met[data$met < 0] <- 0
        data$met[data$met > 1] <- 1
        cols <- RColorBrewer::brewer.pal(8, "Dark2")[1:length(unique(data$group))]
        names(cols) <- 1:length(unique(data$group))
        comps <- as.numeric(levels(data$group))

        # Densities    
        betaDensity <- numeric(0)
        # etas <- numeric(0)
        for (i in comps) {
            a <- MixtureModelResults$Models[[GeneName]]$a[i]
            b <- MixtureModelResults$Models[[GeneName]]$b[i]
            eta <- MixtureModelResults$Models[[GeneName]]$eta[i]
            betaDensity <- c(betaDensity, dbeta(seq(0,1,0.001), shape1 = a, shape2 = b) * eta)
            # etas <- c(etas, MixtureModelResults$Models[[GeneName]]$eta[i])
        }
        # maxDens <- max(betaDensity)
        # Factor <- maxHist / maxDens
        betaDensity <- data.frame(x = rep(seq(0,1,0.001), length(levels(data$group))),
                                  comp = rep(as.numeric(levels(data$group)), each = 1001),
                                  dens = betaDensity)
                
        # Histogram
        g2 <- 
            ggplot2::ggplot(data, ggplot2::aes(x = met)) + 
            ggplot2::geom_histogram(ggplot2::aes(y = ..density..), breaks = seq(0, 1, by = .02), col = "black", fill = "lightgrey") + 
            ggplot2::scale_y_continuous(name = "Density") +
            ggplot2::scale_x_continuous(name = "DNA Methylation", limits = 0:1) +
            ggplot2::geom_line(data = betaDensity, ggplot2::aes(x = x, y = dens, colour = factor(comp)), size = 1.5) +
            ggplot2::scale_color_manual(values = cols, name="Mixture\ncomponent") +
            ggplot2::theme(axis.text = ggplot2::element_text(size=16),
                           axis.title = ggplot2::element_text(size=16),
                           legend.title = ggplot2::element_text(size=14),
                           legend.text = ggplot2::element_text(size=14),
                           plot.title = ggplot2::element_text(size=18, face="bold")) +
            ggplot2::ggtitle(paste("Mixture model of", GeneName))
        
        # Add line with CI for normal samples
        if (!is.null(METnormal)) {
            g2build <- ggplot2::ggplot_build(g2)
            # maxY <- max(max(g2build$data[[1]]$y), max(betaDensity$dens)) 
            yheight <- g2build$layout$panel_ranges[[1]]$y.range[2] * 0.95
            metnormal <- METnormal[GeneName, , drop=FALSE]
            if (length(metnormal) > 1) {
                tmpTtest1 <- t.test(metnormal)
                g2 <- g2 + ggplot2::geom_segment(ggplot2::aes(x = tmpTtest1$conf.int[1], y = yheight, xend = tmpTtest1$conf.int[2], yend = yheight), color = "black", size = 1.5)
            }
        }
        
        # Scatterplot
        g <- NULL
        if (!is.null(GEcancer)) {
            g <- ggplot2::ggplot(data, ggplot2::aes(x = met, y = ge)) + 
                    ggplot2::geom_point(shape = 19, alpha = 0.4, size = 4, ggplot2::aes(color = group) ) +
                    ggplot2::scale_color_manual(values = cols, name="Mixture\ncomponent") + #, labels=c("Normal", 1:length(unique(data$group)))) +
                    ggplot2::scale_y_continuous(name = "Gene expression") +
                    ggplot2::scale_x_continuous(name = "DNA Methylation", limits = 0:1) + 
                    ggplot2::geom_smooth(method=lm, se=FALSE) +
                    ggplot2::theme(axis.text = ggplot2::element_text(size=16),
                                   axis.title = ggplot2::element_text(size=16),
                                   legend.title = ggplot2::element_text(size=14),
                                   legend.text = ggplot2::element_text(size=14),
                                   plot.title = ggplot2::element_text(size=18, face="bold")) +
                    ggplot2::ggtitle(paste("Correlation for", GeneName))
        }
        return(list(MixtureModelPlot = g2, CorrelationPlot = g))
    } else {
        cat(paste(GeneName, "not found.\n"))
    }
}

#' The MethylMix_Predict function
#' 
#' Given a new data set with methylation data, this function predicts the mixture
#' component for each new sample and driver gene. Predictions are based on posterior 
#' probabilities calculated with MethylMix'x fitted mixture model.
#' 
#' @param newBetaValuesMatrix Matrix with new observations for prediction, 
#' genes/cpg sites in rows, samples in columns. Although this new matrix can have
#' a different number of genes/cpg sites than the one provided as METcancer when
#' running MethylMix, naming of genes/cpg sites should be the same.
#' @param MethylMixResult Output object from MethylMix
#' 
#' @return A matrix with predictions (indices of mixture component), driver genes in rows, new samples in columns
#' @export
#' @examples
#' # load the three data sets needed for MethylMix
#' data(METcancer)
#' data(METnormal)
#' data(GEcancer)
#' 
#' # run MethylMix on a small set of example data
#' MethylMixResults <- MethylMix(METcancer, GEcancer, METnormal)
#' # toy example new data, of same dimension of original METcancer data
#' newMETData <- matrix(runif(length(METcancer)), nrow = nrow(METcancer))
#' rownames(newMETData) <- rownames(METcancer)
#' colnames(newMETData) <- paste0("sample", 1:ncol(METcancer))
#' predictions <- MethylMix_Predict(newMETData, MethylMixResults)
#' 
MethylMix_Predict <- function(newBetaValuesMatrix, MethylMixResult) {
    
    # From new data keep only driver genes
    drivers <- intersect(MethylMixResult$MethylationDrivers, rownames(newBetaValuesMatrix))
    
    # For each driver gene, predict class through all samples
    predictions <- matrix(NA, nrow = length(drivers), ncol = ncol(newBetaValuesMatrix), dimnames = list(drivers, colnames(newBetaValuesMatrix)))
    for (driver in drivers) {
        if (MethylMixResult$NrComponents[[driver]] > 1) {
            predictions[driver, ] <- predictOneGene(newBetaValuesMatrix[driver, ], MethylMixResult$Models[[driver]])
        } else {
            predictions[driver, ] <- 1  # driver gene has only one component, no prediction to make
        }
        
    }
    return(predictions)
}

#' The predictOneGene function
#' 
#' Auxiliar function. Given a new vector of beta values, this function calculates a matrix with posterior prob of belonging 
#' to each mixture commponent (columns) for each new beta value (rows), 
#' and return the number of the mixture component with highest posterior probabilit
#' 
#' @param newVector vector with new beta values
#' @param mixtureModel beta mixture model object for the gene being evaluated.
#' 
#' @return A matrix with predictions (indices of mixture component), driver genes in rows, new samples in columns
#' 
predictOneGene <- function(newVector, mixtureModel) {
    densities <- sapply(newVector, function(newObs) {
        mapply(dbeta, newObs, mixtureModel$a, mixtureModel$b) # density value in each mixture
    })
    numerator <- apply(densities, 2, `*`, mixtureModel$eta)
    denominator <- colSums(numerator)
    posteriorProb <- apply(numerator, 1, function(x) x / denominator)
    # If densities are all 0 (can happen for x near 0 or 1), denominator is 0, assign mixing prop as posterior prob
    idx <- which(denominator == 0)
    if (length(idx) > 0) posteriorProb[idx, ] <- matrix(mixtureModel$eta, nrow = length(idx), ncol = length(mixtureModel$eta), byrow = T)
    predictedMixtureComponent <- apply(posteriorProb, 1, which.max)
    return(predictedMixtureComponent)
}
