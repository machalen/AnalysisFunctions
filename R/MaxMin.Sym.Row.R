MaxMin.Sym.Row <-
function(Annot.Table, ColField, operator="Max", colName="Symbol"){
  #Annot.Table: matrix to filter with at least one numeric column. 
  #ColField: Numeric column to calculate max and min per gene symbol
  #Operator: Which operator is used to calculate, Max or Min
  #colName= Column with geneIds
  
  cn <- which(colnames(Annot.Table) == colName)
  UniqueSym <- unique(Annot.Table[,cn])
  Annot.Table <- as.matrix(Annot.Table)
  res.Table <- matrix(NA, nrow=length(UniqueSym), ncol=ncol(Annot.Table))
  colnames(res.Table) <- colnames(Annot.Table)
  if(operator == "Max"){
    for (i in 1:length(UniqueSym)) {
      Sym <- UniqueSym[i]
      TabMaxLog <- Annot.Table[Annot.Table[,cn] == Sym,]
      if (nrow(TabMaxLog) > 1){
        res.Table[i,] <- TabMaxLog[max(TabMaxLog[,ColField]) == TabMaxLog[,ColField], ]
      }else{
        res.Table[i,] <- TabMaxLog
      }
    }
  }
  
  if(operator == "Min"){
    for (i in 1:length(UniqueSym)) {
      Sym <- UniqueSym[i]
      TabMaxLog <- Annot.Table[Annot.Table[,cn] == Sym,]
      if (nrow(TabMaxLog) > 1){
        res.Table[i,] <- TabMaxLog[min(TabMaxLog[,ColField]) == TabMaxLog[,ColField], ]
      }else{
        res.Table[i,] <- TabMaxLog
      }
    }
  }
  
  return(res.Table)
}
