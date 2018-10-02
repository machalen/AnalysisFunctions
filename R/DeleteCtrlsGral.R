DeleteCtrlsGral <-
function (estimates, AnnotfilePath){
  #estimates: expression matrix
  #AnnotfilePath: path to the affy annotation file in .csv format
  require(data.table)
  annotFile <- fread(AnnotfilePath,skip="transcript_cluster_id",header=TRUE,stringsAsFactors=FALSE)
  cjt_annot_sense_chr<-annotFile[annotFile$seqname=="---" | annotFile$seqname=="",]
  return(estimates[!(rownames(estimates) %in% cjt_annot_sense_chr$transcript_cluster_id),])
}
