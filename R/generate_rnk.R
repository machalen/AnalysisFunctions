generate_rnk <-
function(RNKfile,col.names ,rnk) {
  #RNKfile: data frame with columns: Symbol, P.Value and logFC
  #col.names: Vector with the name of the columns: "Symbol", "P.Value" and "logFC"
  #rnk: Path complet and file name with .rnk termination
  distrRank <- (-log10(RNKfile[,col.names[2]]))*(RNKfile[,col.names[3]]/abs(RNKfile[,col.names[3]]))
  RNKfiler <- data.frame(Symbol=RNKfile[,col.names[1]],distrRank )
  dfbase2 <- aggregate(. ~ Symbol, data = RNKfiler, mean)  #Unificar els genesymbols duplicats i fer la mitjana
  #range(dfbase2[,2]) #-3.517880  3.617405
  write.table(dfbase2,file=rnk,quote = FALSE, sep="\t",
              row.names=FALSE, col.names = FALSE)
}
