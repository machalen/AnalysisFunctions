Unique.Sym.Estimates <-
function(Annot.Table, operator="Mean", colName="Symbol"){
  #Annot.Table: Numeric matrix with a unique column of characters corresponding to the gene symbols
  #Operator: Mean o Median to calculate for the intensities in the case of repeated genes
  #colName: Name of the column containing characters (gene symbols)
  cn <- which(colnames(Annot.Table) == colName)
  UniqueSym <- unique(Annot.Table[,cn])
  UniqueSym <- UniqueSym[!is.na(UniqueSym)]
  Annot.Table <- as.data.frame(Annot.Table, stringsAsFactors=F)
  res.Table <- data.frame()
  
  if (operator == "Mean"){
    for(i in 1:length(UniqueSym)){
      gene <- UniqueSym[i]
      redTab <- Annot.Table[(Annot.Table[,cn] == gene) & !is.na(Annot.Table[,cn]),]
      if(nrow(redTab) > 1){
        res.Table <- rbind(res.Table,
                           colMeans(redTab[,which(colnames(Annot.Table) != colName)]))
        colnames(res.Table) <- colnames(Annot.Table)[which(colnames(Annot.Table) != colName)]
      }else{
        res.Table <- rbind(res.Table, redTab[,which(colnames(Annot.Table) != colName)])
        #colnames(res.Table) <- colnames(Annot.Table)[which(colnames(Annot.Table) != colName)]    
      }
    }
  } 
  if (operator == "Median") {
    for(i in 1:length(UniqueSym)){
      gene <- UniqueSym[i]
      redTab <- Annot.Table[(Annot.Table[,cn] == gene) & !is.na(Annot.Table[,cn]),]
      if(nrow(redTab) > 1){
        res.Table <- rbind(res.Table,
                           apply(redTab[,which(colnames(Annot.Table) != colName)],
                                 2, FUN = median))
        colnames(res.Table) <- colnames(Annot.Table)[which(colnames(Annot.Table) != colName)]
      }else{
        res.Table <- rbind(res.Table,redTab[,which(colnames(Annot.Table) != colName)])
        #colnames(res.Table) <- colnames(Annot.Table)[which(colnames(Annot.Table) != colName)]
      }
    }
  }
  
  rownames(res.Table) <- UniqueSym
  return(res.Table)
}
