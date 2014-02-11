#!/usr/bin/env Rscript
# Modified from Nov14.R by Jeff Phillips, 1-20-2014.
# Includes a different mask threshold for sub003 than other subjects,
# to avoid SVD from failing to converge.
########################################################
# this assumes that we know the "block" stimuli   ######
# which allows us to majority vote over a block   ######
# to identify the final predictor for that time   ######
########################################################
library(caret)
library(randomForest)
library(ANTsR)
library(pheatmap)
library(e1071)
########################################################
# Global variables:
studydir<-"/Users/jeff/krns/haxby/openfmri"
########################################################
# Declare function fmriPredictorMatrix.
# Input: BOLD 4D dataset, mask file, motion parameters, 
fmriPredictorMatrix <- function( fmri, mask , motionin , selector = NA , ncompcor= 3 , nsvdcomp = 4 )
{
  mat<-timeseries2matrix( fmri, mask )
  if ( sum(is.na( selector))>0 ) selector<-rep(TRUE, 1:nrow( motionin ))
  mat<-mat[ selector, ]
  msvd <- svd(t(motionin[selector,3:ncol(motionin) ]))
  motion <-as.matrix(msvd$v[, 1:nsvdcomp])
  cmpc<-as.matrix( compcor( mat, ncompcor = ncompcor ))
  globsig<-apply( mat , FUN=mean , MARGIN=1 )
  if ( length(c(cmpc)) == nrow(mat)) mat<-residuals( lm( mat ~ as.matrix(motion)   )) else mat<-residuals( lm( mat ~ as.matrix(motion)  ))
#  if ( length(c(cmpc)) == nrow(mat)) mat<-residuals( lm( mat ~ as.matrix(motionin[selector,])   )) else mat<-residuals( lm( mat ~ as.matrix(motionin[selector,])  ))
  mydot<-F
  if ( mydot ) mat<-temporalwhiten( mat )
  return( list( mat = mat, cmpc=cmpc, globsig=globsig )  )
}
########################################################
majoritylabel <- function( groundtruth, myprediction )
  {
  y  <- myprediction
  yt <- groundtruth
  ct<-1
  mylabs<-unique(  yt  )
  myyruns<-rep(NA,length(yt))
  myyruns[1]<-ct
  for ( i in 1:(length(y)-1))
    {
    if ( yt[i] == yt[i+1] ) { myyruns[i+1]<-ct }  else { ct<-ct+1  ;   myyruns[i+1] <- ct }
    }
  myvotedlabels<-data.frame( groundtruth=rep(NA,ct) , voted=rep(NA,ct))
  for ( i in unique( myyruns ))
    {
    mylabs<-y[ myyruns == i ] 
    ff<-as.data.frame(table(mylabs))
    mylab<-ff$mylabs[ which.max( ff$Freq ) ]
    print( paste( i, mylab ,yt[myyruns == i][1]))
    myvotedlabels$groundtruth[i]<-as.character(yt[myyruns == i][1])
    myvotedlabels$voted[i]<-as.character(mylab)
    }
  return(myvotedlabels)
  }

# Create a single design file called "labels.txt" for each subject.
# Note that this is shifted forward in time relative to the PyMVPA version of this study's data.
#condkey<-read.table(paste(studydir,"/models/model001/condition_key.txt",sep=""),header=F,colClasses="character")
#names(condkey)<-c("task","cond","label")
#design<-data.frame(labels=rep("rest",1452),chunks=rep(0:11,each=121),stringsAsFactors=FALSE)
#for (cond in 1:8) {
#	for (run in 1:12) {
#		tmp<-read.table(paste("model/model001/onsets/task001_run",sprintf("%03d",run),"/cond",sprintf("%03d",cond),".txt",sep=""),header=F)
#		vols<-121*(run-1)+unique(floor(tmp$V1/2.5))
#		vols<-121*(run-1)+which(approx(x=seq(0,300,by=2),y=(seq(0,300,by=2) %in% tmp$V1),xout=2.5*(0:120))$y>0.25)[1:10]
#		h<-hist(tmp$V1/2.5,breaks=0:121,plot=FALSE)
#		vols<-121*(run-1)+which(h$counts>0)[1:9]
#		design$labels[vols]<-condkey$label[condkey$cond==paste("cond",sprintf("%03d",cond),sep="")]
#	}
#}
#design$labels<-as.factor(design$labels)
#write.table(file="labels.txt",design,row.names=F,sep="\t",quote=F)

if ( ! exists("myrates")) myrates<-rep(NA,12)
design<-read.table('labels.txt',header=T)
#unique(design$chunks)
runstotest<-unique(design$chunks)
runstotest<-runstotest[ runstotest < 12 ] 
if ( ! file.exists("AFFINE.nii.gz"))
  {
  print("FAILURE --- you need to be within a subject's directory")
  q()
  }
fmri<-antsImageRead("AFFINE.nii.gz",4)
fmriavg<-antsImageRead("AFFINE_avg.nii.gz",3)
motionin<-read.csv('AFFINEMOCOparams.csv')
maskFull<-getMask(fmriavg,125,1.e9,TRUE)
for ( wrun in runstotest )
{
#  print(paste("Processing run ", wrun, "...",sep=""))
  mask<-antsImageRead('mask.nii.gz',3)
  selector<-( as.numeric( design$chunks ) != wrun   )
  selector2<-( as.numeric( design$chunks ) == wrun  )
  subdesign<-subset( design, selector )
  ncc <- 4
  print(paste("milepost 1",wrun,sep=", "))
  fmripreds<-fmriPredictorMatrix( fmri, mask, motionin, selector, ncompcor = ncc )
  print(paste("milepost 2",wrun,sep=", "))  
  fmripredsFull<-fmriPredictorMatrix( fmri, maskFull, motionin, selector, ncompcor = ncc )
  print(paste("milepost 2.5",wrun,sep=", "))  
  mat<-residuals( lm( fmripreds$mat ~ fmripredsFull$cmpc ))
  myclasses <- levels( subdesign$labels )
  nclasses<-length(myclasses )
  myblocks<-matrix( rep(0,(nrow(mat))*nclasses ), nrow=( nrow(mat)  ))
  for ( i in 1:nclasses ) myblocks[,i]<-as.numeric(  subdesign$labels == myclasses[i] )
  mysblocks<-myblocks
  for ( i in 1:ncol(mysblocks)) mysblocks[,i]<-predict(smooth.spline(mysblocks[,i],df=100))$y
  mydesign<-cbind( myblocks, mysblocks )
  nv<-85
  ff<-svd( mat )
  mysccanimages<-t(ff$v[,1:nv])
  mysccanpreds <- ( mat  ) %*% t( mysccanimages )
  mydf         <- data.frame( factpreds = as.factor((subdesign$labels))  , imgs = mysccanpreds )
  my.rf        <- svm( factpreds ~ . , data=mydf, probability = TRUE  )
#######################
##### test phase ######
#######################

  print(paste("milepost 3",wrun,sep=", "))
  fmriTest     <- fmriPredictorMatrix( fmri, mask,     motionin, selector2, ncompcor = ncc )
  print(paste("milepost 4",wrun,sep=", "))
  fmriTestFull <- fmriPredictorMatrix( fmri, maskFull, motionin, selector2, ncompcor = ncc )
  mat2<-residuals( lm( fmriTest$mat ~ fmriTestFull$cmpc ))
  mysccanpreds2 <- ( mat2  ) %*% t( mysccanimages )
  mydf2<-data.frame(  imgs = mysccanpreds2 ) # , corrs = cor( t( mysccanpreds2 )) )
  mypred2<-predict( my.rf , newdata = mydf2 )
  subdesign2<-subset( design, selector2 )
  sublabels<-as.factor((subdesign2$labels))
  zz<-majoritylabel( sublabels , mypred2 )
  myrate<-100*(sum(zz$groundtruth==zz$voted)/length(zz$groundtruth))
  print(paste("CorrectClassify:",myrate,"%",getwd(),wrun+1))
  myrates[ wrun+1 ]<-myrate
} # wrun loop
############################################################################################
ratedf<-data.frame( RunNumber=c(0:11), CrossValidatedPredictionForRun=myrates )
write.csv(ratedf,'mypredictionresults.csv',row.names=F)
