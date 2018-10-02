generate_cls <-
function(class, cls) {
  #class: pData(my.targets)$(columna amb els grups que es volen comparar)
  #cls: cls file path where the results are stored + name of the file with .cls ending
  #S'han d'eliminar els espais finals de cada row potser
  #Definir les rows
  class <- as.character(class)
  levels <- unique(class)
  firstr <- c(length(class), length(levels),1)
  secndr <- c("#",levels)
  #Les condicions es poden posar literalment o amb numeros, s'ha canviat 01/03/2016
  #numbers <- c(0:(length(levels)-1))
  #map = setNames(numbers, levels)
  #thirdr <- map[class]
  thirdr <- class
  #Crear el fixer
  cat(firstr, '\n', file = cls)
  cat(secndr, '\n', file = cls, append=TRUE)
  cat(thirdr, '\n', file = cls, append=TRUE)
  
}
