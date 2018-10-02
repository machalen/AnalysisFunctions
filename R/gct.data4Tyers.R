gct.data4Tyers <-
function(data4Tyers, Symb.col="AffyID", gct){
  #data4Tyers: data4Tyers matrix used to create the gct file
  #Symb.col: Column with the used gene.id
  #gct: file path name where the results are stored + name of the file with .gct ending
  matriu <- data4Tyers[,grep(".RMA$",colnames(data4Tyers))]
  colnames(matriu) <- gsub(".RMA", "", colnames(matriu))
  rows <- nrow(matriu)
  cols <- ncol(matriu)
  #Generate the headers
  cat("#1.2", '\n', sep="", file = gct)
  cat(rows, '\t', sep="", file = gct, append=TRUE)
  cat(cols, '\n', sep="", file = gct, append=TRUE)
  
  #Generate the columns of data
  colnm <- c("NAME", "Description", colnames(matriu))
  col1 <- data4Tyers[,Symb.col]
  col2 <- rep("NA", rows)
  gctfile <- rbind(colnm, data.frame(col1, col2, matriu, stringsAsFactors=FALSE))
  #Append the data to the file
  write.table(gctfile, file=gct, append = T, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}
