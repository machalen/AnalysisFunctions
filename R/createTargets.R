createTargets <-
function(names,fileName="Targets.txt", DestDir=analysisDir){
  #180915 afegeixo el nom del fitxer com a parametre d'entrada
  df<-data.frame("FileName"=names)
  write.table(df, file=file.path(DestDir,fileName),dec=",",sep="\t",row.names=FALSE, quote=FALSE)
}
