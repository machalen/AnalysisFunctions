generate_gct <-
function(matriu, gct){
  #Matriu: est_noctrls object
  #gct: file path name where the results are stored + name of the file with .gct ending
  rows <- nrow(matriu)
  cols <- ncol(matriu)
  #Generate the headers
  cat("#1.2", '\n', sep="", file = gct)
  cat(rows, '\t', sep="", file = gct, append=TRUE)
  cat(cols, '\n', sep="", file = gct, append=TRUE)
  
  #Generate the columns of data
  colnm <- c("NAME", "Description", colnames(matriu))
  col1 <- rownames(matriu)
  col2 <- rep("NA", rows)
  gctfile <- rbind(colnm, data.frame(col1, col2, matriu, stringsAsFactors=FALSE))
  #Append the data to the file
  write.table(gctfile, file=gct, append = T, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}
