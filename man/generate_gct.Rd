\name{generate_gct}
\alias{generate_gct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate gct format files
}
\description{
This function generates .gct format files to use in the GSEA analysis.
}
\usage{
generate_gct(matriu, gct)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matriu}{
Numeric expression matrix with rows genes and columns samples
}
  \item{gct}{
Path and name of the output file with .gct termination
}
}
\details{
The .gct file contains information about the expression values of genes to be considered in the GSEA analysis. The order of samples in the columns of the matrix has to be coincident with the name of the conditions in the cls file. 
}
\value{
.gct format file
}
\references{
https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
}
\author{
Magdalena Arnal Segura <marnal@imim.es>
}
