### MethylMix.R
### MethylMix identifies differential and functional methylation drivers. 
### Author: Olivier Gevaert
### Date: May 2, 2014

MethylMix <- function(METcancer,METnormal,MAcancer,OutputRoot='',Parallel=FALSE) { 
     
     # Step 1: modeling the gene expression using cancer methylation data
     FunctionalGenes=MethylMix_ModelGeneExpression(METcancer,MAcancer,Method="Regression",NULL)
     
     # Step 2: modeling the methylation data as a mixture of beta's.     
     if (length(FunctionalGenes)>0) {
          MixtureModelResults=MethylMix_MixtureModel(METcancer,METnormal,FunctionalGenes,Parallel,1,length(FunctionalGenes))
          # writing to file
          if (length(OutputRoot)>1) {
               MethylMix_WriteToFile(OutputRoot,MixtureModelResults)
          }
          return(MixtureModelResults)
     } else {
          cat("No functional genes found.\n")
     }
}

MethylMix_ModelGeneExpression <- function(METcancer,MAcancer,Method=c("Regression","Pearson","Spearman"),CovariateData=NULL) {      
     
     # overlapping samples     
     OverlapSamples=intersect(colnames(METcancer),colnames(MAcancer))
     cat("Found",length(OverlapSamples),"samples with both methylation and expression data.\n")
     #cat(paste("Found",length(OverlapSamples),"samples with both methylation and expression data."),sep=' ')
     MAcancer=MAcancer[,OverlapSamples,drop=FALSE]
     METcancer=METcancer[,OverlapSamples,drop=FALSE]
     if (!is.null(CovariateData)) CovariateData=as.matrix(CovariateData[OverlapSamples,])
     
     Rsquares=matrix(0,nrow=length(rownames(METcancer)),ncol=1)    
     Genes=rownames(METcancer)  
     PvalueThreshold=0.001  
     RsquareThreshold=0.1
     
     pb=txtProgressBar(1,length(rownames(METcancer)),style=3,width=100)
     cat("Correlating methylation data with gene expression.\n")
     for(i in 1:length(rownames(METcancer))) {
          #cat(paste(tmpGene,"\r"))
          setTxtProgressBar(pb,i)
          tmpGene=Genes[i]
          tmpGene=unlist(strsplit(tmpGene,'---'))
          tmpGene=tmpGene[1]  
          pos=which(rownames(MAcancer)==tmpGene)
          
          if (length(pos)>0) {
               if (!is.null(CovariateData)) {
                    res=lm(MAcancer[pos,]~METcancer[Genes[i],] + CovariateData)
               } else {           
                    res=lm(MAcancer[pos,]~METcancer[Genes[i],])
               }
               res.summary=summary(res)
               if (!is.null(CovariateData)) {
                    if (res$coefficients[2]<0 & res.summary$coefficients[2,4]<PvalueThreshold & abs(res.summary$coefficients[2,1])>abs(res.summary$coefficients[3,1]) ) { # methylation effect bigger than tissue
                         
                         Rsquares[i]=res.summary$r.squared
                    }
               } else if (res$coefficients[2]<0 & res.summary$coefficients[2,4]<PvalueThreshold) {
                    Rsquares[i]=res.summary$r.squared    
               }
          }
     }
     
     # Rsquare threshold
     FunctionalGenes=Genes[Rsquares>RsquareThreshold]
     cat("\nFound",length(FunctionalGenes),"functional genes.\n")
     return(FunctionalGenes)
}


MethylMix_MixtureModel <- function(METcancer,METnormal,FunctionalGenes,Parallel=FALSE,Start,Stop,NoNormalMode=FALSE) {
     
     # overlap of samples
     GeneOverlap=intersect(rownames(METcancer),rownames(METnormal))
     METcancer=METcancer[GeneOverlap,,drop=FALSE]
     METnormal=METnormal[GeneOverlap,,drop=FALSE]
     
     if (Stop>length(FunctionalGenes)) Stop=length(FunctionalGenes)
     cat("\nStarting Beta mixture modeling.\n")
     METcancer=METcancer[which(rownames(METcancer) %in% FunctionalGenes),,drop=FALSE]  
     METnormal=METnormal[which(rownames(METnormal) %in% FunctionalGenes),,drop=FALSE]
     
     MethylationStates=matrix(0,length(rownames(METcancer)),length(colnames(METcancer)))    
     rownames(MethylationStates)=rownames(METcancer)
     colnames(MethylationStates)=colnames(METcancer)
     
     AllNrComponents=matrix(1,length(rownames(METcancer)),1)
     AllFlipOverStates=matrix(1,length(rownames(METcancer)),1)
     AllMixtureStates=array(list(),length(rownames(METcancer)))
     Models=array(list(),length(rownames(METcancer)))
     GeneNames=rownames(MethylationStates)
     
     cat("Running Beta mixture model on",length(rownames(METcancer)),"functional genes and on",length(colnames(METcancer)),"samples.\n")
     #cat("Running genes from position",Start,"to",Stop,".\n")
     if (Parallel & .Platform$OS.type!="windows") {
#           tmp=tryCatch(
#                {
                    # registering an appropriate nr of cores
                    NrCores=detectCores()
                    NrCores=NrCores-2
                    if (NrCores<2) NrCores=1
                    registerDoParallel(cores=NrCores) # this will register nr of cores/threads
                    
                    # alternative to do it in parallel
                    tmpResults=foreach::foreach(i=Start:Stop) %dopar% {      #length(rownames(METcancer))               
                         MixtureModelResults_SingleGene=MethylMix_ModelSingleGene(GeneNames[i],METcancer[i,],METnormal[i,],NoNormalMode)
                    }
                    for (i in 1:(Stop-Start+1)) {   #length(rownames(METcancer))
                         MethylationStates[Start+i-1,]=tmpResults[[i]]$MethylationState
                         Models[[Start+i-1]]=tmpResults[[i]]$Model
                         AllNrComponents[Start+i-1,1]=tmpResults[[i]]$NrComponents
                         AllFlipOverStates[Start+i-1,1]=tmpResults[[i]]$FlipOverState
                         AllMixtureStates[[Start+i-1]]=tmpResults[[i]]$MixtureStates
#                     }  
#                }, finally = {
#                     if (.Platform$OS.type=="windows") {
#                          stopImplicitCluster()
#                     }
                }
#           )
     } else {
          #pb=txtProgressBar(1,length(rownames(METcancer)),style=3,width=100)
          for(i in Start:Stop) {      #length(rownames(METcancer))               
               #setTxtProgressBar(pb,i)
               MixtureModelResults_SingleGene=MethylMix_ModelSingleGene(GeneNames[i],METcancer[i,],METnormal[i,],NoNormalMode)
               MethylationStates[i,]=MixtureModelResults_SingleGene$MethylationState
               Models[[i]]=MixtureModelResults_SingleGene$Model
               AllNrComponents[i,1]=MixtureModelResults_SingleGene$NrComponents
               AllFlipOverStates[i,1]=MixtureModelResults_SingleGene$FlipOverState
               AllMixtureStates[[i]]=MixtureModelResults_SingleGene$MixtureStates
          }   
     }   
     
     # removing the genes without any differential methylation. 
     NonZeroPositions=rowSums(MethylationStates)!=0
     AllNrComponents=as.matrix(AllNrComponents[NonZeroPositions,])
     AllFlipOverStates=as.matrix(AllFlipOverStates[NonZeroPositions,])
     FunctionalGenes=FunctionalGenes[NonZeroPositions]
     rownames(AllNrComponents)=FunctionalGenes
     
     AllMixtureStates=AllMixtureStates[NonZeroPositions]
     Models=Models[NonZeroPositions]
     MethylationStates=MethylationStates[NonZeroPositions,,drop=FALSE]
     MethylationDrivers=rownames(MethylationStates)
     
     return(list(MethylationStates=MethylationStates,NrComponents=AllNrComponents,Models=Models,MethylationDrivers=MethylationDrivers,MixtureStates=AllMixtureStates))
}


MethylMix_ModelSingleGene <-function(GeneName,METdataVector,METdataNormalVector,NoNormalMode=FALSE) {
     
     OrigOrder=order(METdataVector)
     ChiSquareThreshold=qchisq(0.95,df=3)
     PvalueThreshold=0.01
     MeanDifferenceTreshold=0.10
     
     w0.m=matrix(1,length(METdataVector),1)
     Model1=blc_2(matrix(METdataVector,ncol=1),w=w0.m,maxiter = 100, tol = 1e-06,verbose=FALSE)
     BIC1=-2*Model1$llike+2*log(length(METdataVector))
     
     # modeling with mixture of two beta's
     #cat("Building mixture model with two components.")
     NrClusters=2
     w0.m=matrix(0,nrow=length(METdataVector),ncol=NrClusters)
     w0.m[which(METdataVector < median(METdataVector)),1] <- 1
     w0.m[which(METdataVector >= median(METdataVector)),2] <- 1        
     Model2 <- blc_2(matrix(METdataVector,ncol=1),w=w0.m,maxiter = 100, tol = 1e-06,verbose=FALSE)
     if (sum(is.na(Model2$mu))>0) Model2$llike=0
     BIC2=-2*Model2$llike+5*log(length(METdataVector))
     
     MethylationState=matrix(0,1,length(METdataVector))
     FlipOverState=0
     
     Model=list()
     NrComponents=1                   
     
     if ((BIC1-BIC2)>ChiSquareThreshold & abs(Model2$mu[1]-Model2$mu[2])>MeanDifferenceTreshold) {
          #cat("Building mixture model with three components.\n")
          NrClusters=3        
          tmpQuantiles=quantile(METdataVector,c(0.33,0.66))
          w0.m=matrix(0,nrow=length(METdataVector),ncol=NrClusters)
          w0.m[which(METdataVector < tmpQuantiles[1]),1] <- 1
          w0.m[which(METdataVector >= tmpQuantiles[1] & METdataVector<tmpQuantiles[2]),2] <- 1
          w0.m[which(METdataVector >= tmpQuantiles[2]),3] <- 1   
          
          Model3 <- blc_2(matrix(METdataVector,ncol=1),w=w0.m,maxiter = 100, tol = 1e-06,verbose=FALSE) 
          if (sum(is.na(Model3$mu))>0) Model3$llike=0
          BIC3=-2*Model3$llike+8*log(length(METdataVector))
          
          
          SortedMeans=sort(Model3$mu)
          if ( (BIC2-BIC3)>ChiSquareThreshold & abs(SortedMeans[2]-SortedMeans[1])>MeanDifferenceTreshold & abs(SortedMeans[3]-SortedMeans[2])>MeanDifferenceTreshold ) { # three components are the best
               res=list(p.value=1)
               MixtureStates=matrix(0,3,1)
               METdataVector_comp1=METdataVector[Model3$w[,1]>(1/3)]          
               if (length(METdataVector_comp1)>0) res=wilcox.test(METdataVector_comp1,METdataNormalVector) else res$p.value=1 
               Difference=mean(METdataVector_comp1)-mean(METdataNormalVector)
               if ((res$p.value<PvalueThreshold & abs(Difference)>MeanDifferenceTreshold) | NoNormalMode) {
                    MethylationState[1,Model3$w[,1]>(1/3)]=Difference            
                    MixtureStates[1,1]=Difference
               }
               
               METdataVector_comp2=METdataVector[Model3$w[,2]>(1/3)]
               if (length(METdataVector_comp2)>0) res=wilcox.test(METdataVector_comp2,METdataNormalVector) else res$p.value=1
               Difference=mean(METdataVector_comp2)-mean(METdataNormalVector)
               if ((res$p.value<PvalueThreshold & abs(Difference)>MeanDifferenceTreshold) | NoNormalMode) {
                    MethylationState[1,Model3$w[,2]>(1/3)]=Difference            
                    MixtureStates[2,1]=Difference
               }
               
               METdataVector_comp3=METdataVector[Model3$w[,3]>(1/3)]
               if (length(METdataVector_comp3)>0) res=wilcox.test(METdataVector_comp3,METdataNormalVector) else res$p.value=1
               Difference=mean(METdataVector_comp3)-mean(METdataNormalVector)
               if ((res$p.value<PvalueThreshold & abs(Difference)>MeanDifferenceTreshold) | NoNormalMode) {
                    MethylationState[1,Model3$w[,3]>(1/3)]=Difference            
                    MixtureStates[3,1]=Difference
               }
               cat(c(GeneName,": Three components are best.\n"))
               NrComponents=3
               Model=Model3
               
          } else { # two components are the best
               # need to compare both components to normal
               # check first state
               MixtureStates=matrix(0,2,1)
               METdataVector_comp1=METdataVector[Model2$w[,1]>0.5]
               res=wilcox.test(METdataVector_comp1,METdataNormalVector)
               Difference=mean(METdataVector_comp1)-mean(METdataNormalVector)
               if ((res$p.value<PvalueThreshold & abs(Difference)>MeanDifferenceTreshold) | NoNormalMode) {
                    MethylationState[1,Model2$w[,1]>0.5]=Difference            
                    MixtureStates[1,1]=Difference
               }
               
               METdataVector_comp2=METdataVector[Model2$w[,2]>0.5]
               res=wilcox.test(METdataVector_comp2,METdataNormalVector)
               Difference=mean(METdataVector_comp2)-mean(METdataNormalVector)
               if ((res$p.value<PvalueThreshold & abs(Difference)>MeanDifferenceTreshold) | NoNormalMode) {
                    MethylationState[1,Model2$w[,2]>0.5]=Difference            
                    MixtureStates[2,1]=Difference
               }
               cat(c(GeneName,": Two components are best.\n"))
               NrComponents=2
               Model=Model2
               
               # correcting the flipover effect
               FlipOverResults=MethylMix_RemoveFlipOver(OrigOrder,MethylationState)
               MethylationState=FlipOverResults$MethylationState
               FlipOverState=FlipOverResults$LearnedState     
               
          }
     } else { # 1 component is best
          # compare with normal
          MixtureStates=matrix(0,1,1)
          res=wilcox.test(METdataVector,METdataNormalVector)
          Difference=mean(METdataVector)-mean(METdataNormalVector)
          if ((res$p.value<PvalueThreshold & abs(Difference)>MeanDifferenceTreshold) | NoNormalMode) {
               MethylationState[1,]=Difference # was mean(METdataVector), but should be Difference????
               MixtureStates[1,1]=Difference
          }
          cat(c(GeneName,": One component is best.\n"))
          Model=Model1
     }                 
     return(list(MethylationState=MethylationState,NrComponents=NrComponents,Model=Model,MixtureStates=MixtureStates,FlipOverState=FlipOverState))
}

MethylMix_RemoveFlipOver <-function(OrigOrder,MethylationState,UseTrainedFlipOver=FALSE,FlipOverState=0) {
     
     Differences=diff(MethylationState[1,OrigOrder])
     DifferencesSel=Differences[Differences!=0]
     LearnedState=0
     if (length(DifferencesSel)==2) {
          if (DifferencesSel[1]*-1==DifferencesSel[2]) {     
               
               posDiff1=which(Differences==DifferencesSel[1])
               stateSize1=posDiff1
               posDiff2=which(Differences==DifferencesSel[2])
               stateSize2=length(Differences)-posDiff2
               
               if (UseTrainedFlipOver==TRUE) {
                    if (FlipOverState==1) {
                         MethylationState[1,OrigOrder[posDiff2+1:length(Differences)]]=MethylationState[1,OrigOrder[posDiff2]]
                    } else if (FlipOverState==2) {
                         MethylationState[1,OrigOrder[1:posDiff1]]=MethylationState[1,OrigOrder[posDiff1+1]]
                    }                    
               } else {
                    if (stateSize2>stateSize1) {                    
                         MethylationState[1,OrigOrder[1:posDiff1]]=MethylationState[1,OrigOrder[posDiff1+1]]
                         LearnedState=2
                    } else if (stateSize1>stateSize2) {                    
                         MethylationState[1,OrigOrder[posDiff2+1:length(Differences)]]=MethylationState[1,OrigOrder[posDiff2]]
                         LearnedState=1
                    }  
               }
          }
     }
     return(list(MethylationState=MethylationState,LearnedState=LearnedState))
}  

MethylMix_PlotModel <-function(GeneName,METdata,MixtureModelResults,MAdata=0,METnormal=0,FileName="") {
     
     GeneNameMA=GeneName
     GeneNameMA=unlist(strsplit(GeneNameMA,'---'))
     GeneNameMA=GeneNameMA[1]
     
     
     Pos=which(rownames(MixtureModelResults$MethylationStates) %in% GeneName)
     SecondPos=which(rownames(METdata) %in% GeneName)
     if (length(Pos)>0 & length(SecondPos)>0) {
          
          BetaModel=MixtureModelResults$Models[[Pos]]
          NrComponents=MixtureModelResults$NrComponents[Pos]
          
          Pos2=which(rownames(METdata) %in% GeneName)
          METdataVector=METdata[Pos2,]
          
          # start plotting
          #par(mfrow=c(2,1))          
          
          if (FileName!="") {
               File=paste(FileName,"MixtureModel.tif",sep='')
               tiff(filename = File, compression="none",pointsize=20,bg="white",width=600,height=600)                     
          }
          
          histogram=hist(METdataVector,50,plot=FALSE)
          HistHeight=max(histogram$counts)*1.2
          histogram=hist(METdataVector,50,xlim=c(0,1),ylim=c(0,1.4*HistHeight),main=paste("Mixture model of",GeneName),xlab="DNA methylation",ylab="Frequency")
          
          # plotting the normal 95% confidence interval, if the normal data is present
          if (length(METnormal)>1) {
               tmpTtest=t.test(METnormal[GeneName,])            
               NormalMean=mean(METnormal[GeneName,])
               lines(c(tmpTtest$conf.int[1],tmpTtest$conf.int[2]),c(1.2*HistHeight,1.2*HistHeight),ylim=c(0,1.2*HistHeight),lwd=5,type='l',col=1)
               points(NormalMean,1.2*HistHeight,ylim=c(0,1.2*HistHeight),lwd=5,type='l',col=1,pch=8)
          }   
          if(NrComponents==1) { 
               MaxHeight=max(dbeta(seq(0,1,0.001),BetaModel$a[1],BetaModel$b[1]))
               Factor=HistHeight/MaxHeight
               lines(seq(0,1,0.001),Factor*dbeta(seq(0,1,0.001),BetaModel$a[1],BetaModel$b[1]),ylim=c(0,1.2*HistHeight),lwd=5,type='l',col=2,xlab="DNA methylation",ylab="Frequency")
               
               if (length(MAdata)>1) {
                    OverlapSamples=intersect(colnames(METdata),colnames(MAdata))
                    if (FileName!="") {
                         dev.off()
                         File=paste(FileName,"Expression.tif",sep='')
                         tiff(filename = File, compression="none",pointsize=20,bg="white",width=600,height=600)                     
                    }
                    StateVector=MixtureModelResults$MethylationStates[GeneName,OverlapSamples]
                    States=MixtureModelResults$MixtureStates[[Pos]]
                    StateVector[]=States[1] 
                    plot(METdata[GeneName,StateVector==States[1]],MAdata[GeneNameMA,StateVector==States[1]],xlim=c(0,1),col=2,lty=1,xlab="DNA methylation",ylab="Gene expression",main=paste("Expression correlation for",GeneName))
               }
          } else if (NrComponents==2) {               
               MaxHeight=max(dbeta(seq(0,1,0.001),BetaModel$a[1],BetaModel$b[1]))
               MaxHeight=max(MaxHeight,max(dbeta(seq(0,1,0.001),BetaModel$a[2],BetaModel$b[2])))
               Factor=max(HistHeight/MaxHeight)
               lines(seq(0,1,0.001),Factor*BetaModel$eta[1]*dbeta(seq(0,1,0.001),BetaModel$a[1],BetaModel$b[1]),ylim=c(0,1.2*HistHeight),lwd=5,type='l',col=2,main="BetaMixtureModel",xlab="DNA methylation",ylab="Frequency")
               lines(seq(0,1,0.001),Factor*BetaModel$eta[2]*dbeta(seq(0,1,0.001),BetaModel$a[2],BetaModel$b[2]),ylim=c(0,1.2*HistHeight),lwd=5,type='l',col=3,main="BetaMixtureModel",xlab="DNA methylation",ylab="Frequency")
               
               #plotting the Gene expression values
               if (length(MAdata)>1) {
                    OverlapSamples=intersect(colnames(METdata),colnames(MAdata))
                    METdata=METdata[,OverlapSamples,drop=FALSE]
                    MAdata=MAdata[,OverlapSamples,drop=FALSE]
                    if (FileName!="") {
                         dev.off()
                         File=paste(FileName,"Expression.tif",sep='')
                         tiff(filename = File, compression="none",pointsize=20,bg="white",width=600,height=600,)                     
                    }
                    StateVector=MixtureModelResults$MethylationStates[GeneName,OverlapSamples]
                    States=MixtureModelResults$MixtureStates[[Pos]]
                    plot(METdata[GeneName,StateVector==States[1]],MAdata[GeneNameMA,StateVector==States[1]],xlim=c(0,1),col=2,pch=1,bg=2,xlab="DNA methylation",ylab="Gene expression",main=paste("Expression correlation for",GeneName))
                    points(METdata[GeneName,StateVector==States[2]],MAdata[GeneNameMA,StateVector==States[2]],xlim=c(0,1),col=3,pch=1,bg=3)            
               }
               
          } else if (NrComponents==3) {
               
               
               MaxHeight=max(dbeta(seq(0,1,0.001),BetaModel$a[1],BetaModel$b[1]))
               MaxHeight=max(MaxHeight,max(dbeta(seq(0,1,0.001),BetaModel$a[2],BetaModel$b[2])))
               MaxHeight=max(MaxHeight,max(dbeta(seq(0,1,0.001),BetaModel$a[3],BetaModel$b[3])))            
               if (is.na(MaxHeight) | is.infinite(MaxHeight)) MaxHeight=30
               Factor=HistHeight/MaxHeight
               lines(seq(0,1,0.001),Factor*BetaModel$eta[1]*dbeta(seq(0,1,0.001),BetaModel$a[1],BetaModel$b[1]),ylim=c(0,1.2*HistHeight),lwd=5,type='l',col=2,bg=2,main="BetaMixtureModel",xlab="DNA methylation",ylab="Frequency")
               lines(seq(0,1,0.001),Factor*BetaModel$eta[2]*dbeta(seq(0,1,0.001),BetaModel$a[2],BetaModel$b[2]),ylim=c(0,1.2*HistHeight),lwd=5,type='l',col=3,bg=3,main="BetaMixtureModel",xlab="DNA methylation",ylab="Frequency")
               lines(seq(0,1,0.001),Factor*BetaModel$eta[3]*dbeta(seq(0,1,0.001),BetaModel$a[3],BetaModel$b[3]),ylim=c(0,1.2*HistHeight),lwd=5,type='l',col=4,bg=4,main="BetaMixtureModel",xlab="DNA methylation",ylab="Frequency")
               
               #plotting the Gene expression values
               if (length(MAdata)>1) {
                    OverlapSamples=intersect(colnames(METdata),colnames(MAdata))
                    METdata=METdata[,OverlapSamples,drop=FALSE]
                    MAdata=MAdata[,OverlapSamples,drop=FALSE]
                    if (FileName!="") {
                         dev.off()
                         File=paste(FileName,"Expression.tif",sep='')
                         tiff(filename = File, compression="none",pointsize=20,bg="white",width=600,height=600)                     
                    }
                    
                    StateVector=MixtureModelResults$MethylationStates[GeneName,OverlapSamples]
                    States=MixtureModelResults$MixtureStates[[Pos]]
                    plot(METdata[GeneName,StateVector==States[1]],MAdata[GeneNameMA,StateVector==States[1]],xlim=c(0,1),col=2,lty=1,xlab="DNA methylation",ylab="Gene expression",main=paste("Expression correlation for",GeneName))
                    points(METdata[GeneName,StateVector==States[2]],MAdata[GeneNameMA,StateVector==States[2]],xlim=c(0,1),col=3,lty=1)
                    points(METdata[GeneName,StateVector==States[3]],MAdata[GeneNameMA,StateVector==States[3]],xlim=c(0,1),col=4,lty=1)
               }
               
          }
     } else {
          cat("This gene does not exist.\n") 
     }
     if (FileName!="") {
          dev.off()
     }
}


MethylMix_WriteToFile <-function(OutputRoot,MethylMixResults) {
     OutputFile=paste(OutputRoot,'_MethylationStates.txt',sep="")
     write.table(MethylMixResults$MethylationStates,OutputFile,sep="\t")
     
     OutputFile=paste(OutputRoot,"_NrComponents.txt",sep="")
     write.table(MethylMixResults$NrComponents,OutputFile,sep="\t")
     
     # also saving an Robject
     OutputFile=paste(OutputRoot,"_MethylMixResults.RData",sep="")
     save(MethylMixResults,file=OutputFile)
}



blc_2 <- function (Y, w, maxiter = 25, tol = 1e-06, weights = NULL, verbose = TRUE) 
{
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


betaEst_2 <-function (Y, w, weights) 
{
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
     opt <- try(optim(logab, betaObjf, ydata = y, wdata = w, weights = weights, 
                      method = "BFGS"), silent = TRUE) 
     if (inherits(opt, "try-error")) 
          return(c(1, 1))
     exp(opt$par) # if using optimx, exp(as.numeric(opt$par))
}
