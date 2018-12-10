data4TyersMake <- function(est_noctrls.annot, cond, fit.main, contrast){
  #est_noctrls.annot: Matriu est_nocontrols amb anotació de gens annotatetableC
  #cond: vector de condicions de cada una de les mostres utilitzada en el contrast
  #fit.main: model limma utilitzat al DE analysis
  #contrast: vector de vectors length dos amb els contrastos que es fan
  
  #Obtenim els valors dels contrastos amb limma
  ConList <- vector("list", length(contrast)) 
  for (i in 1:length(contrast)) {
    ConList[[i]] <- topTable(fit.main,n=Inf,coef=i, adjust="fdr")[,c("logFC","P.Value","adj.P.Val")]
    ConList[[i]] <- ConList[[i]][order(rownames(ConList[[i]])),]
  }
  
  #Escalem els valors de est_noctrls per al heatmap
  ia <- which(colnames(est_noctrls.annot) == "Chrom")
  est_noctrls <- est_noctrls.annot[, c(1:(ia-1))]
  est_noctrls_s<-est_noctrls[order(rownames(est_noctrls)),] 
  est_centered<-est_noctrls_s-apply(est_noctrls_s,1,mean)
  est_scaled<-est_centered/apply(est_noctrls_s,1,sd)
  #est_scaled.o <- est_scaled[, order(cond)]
  est_scaled.o <- est_scaled
  colnames(est_scaled.o) <- paste(colnames(est_scaled.o),"scaled",sep=".")
  est_scaled.o$AffyID <- rownames(est_noctrls_s)
  
  #Construim la matriu amb els mean per condici?
  u.cond <- unique(cond)
  mean.matrix <- matrix(data= NA, nrow=nrow(est_noctrls), ncol=length(u.cond))
  for (ic in 1:length(u.cond)) {
    uc <- u.cond[ic]
    mean.matrix[,ic] <- apply(est_noctrls_s[,cond==uc],1,mean)
  }
  colnames(mean.matrix) <- paste("mean",u.cond, sep=".")
  
  #Calculem el FC utilitzant els mean calculats anteriorment
  FC.matrix <- matrix(data= NA, nrow=nrow(est_noctrls), ncol=length(contrast))
  col.FC.Names <- vector()
  for (ic in 1:length(contrast)) {
    FC.matrix[,ic] <- 2^abs(ConList[[ic]]$logFC) * sign(ConList[[ic]]$logFC)
    col.FC.Names <- c(col.FC.Names, paste("FC",contrast[[ic]][1], "vs", contrast[[ic]][2], sep="."))
  }
  colnames(FC.matrix) <- col.FC.Names
  
  #Construim la matriu dels topDiff amb els c("FC", "logFC","P.Value","adj.P.Val")
  topDiff.mat <- matrix(data= NA, nrow=nrow(est_noctrls), ncol=length(contrast)*4)
  col.topDiff.Names <- vector()
  for (ic in 1:length(contrast)) {
    icc <- 1+(4*(ic-1))
    cont.name <- paste(contrast[[ic]][1], "vs", contrast[[ic]][2], sep=".")
    topDiff.mat[,icc] <- FC.matrix[,ic]
    topDiff.mat[,icc+1] <- ConList[[ic]][,1]
    topDiff.mat[,icc+2] <- ConList[[ic]][,2]
    topDiff.mat[,icc+3] <- ConList[[ic]][,3]
    col.topDiff.Names <- c(col.topDiff.Names, colnames(FC.matrix)[ic],
                           paste(colnames(ConList[[ic]]), cont.name, sep="."))
  }
  colnames(topDiff.mat) <- col.topDiff.Names
  
  #Ordenem i construim la matriu amb les anotacions
  topDiff.annot<-est_noctrls.annot[order(rownames(est_noctrls.annot)), 
                                   c("Symbol","mrna","UCSC_symbols","GO_biological_process",
                                     "GO_cellular_component", "GO_molecular_function", 
                                     "pathway", "Description","Chrom","Strand","Start","Stop")]
  colnames(topDiff.annot)[c(4:6)] <- c("GOBP", "GOCC", "GOMF")
  colnames(est_noctrls_s) <- paste(colnames(est_noctrls_s),"RMA",sep=".")
  
  #Construim la matriu definitiva
  data4Tyers<-data.frame(est_scaled.o, topDiff.annot,
                         topDiff.mat, mean.matrix,
                         est_noctrls_s)
  return(data4Tyers)
  #write.csv2(data4Tyers,file=paste(resultsDir,"Data4Tyers.csv",sep="/"),row.names=FALSE)
}
