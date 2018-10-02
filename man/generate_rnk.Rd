\name{generate_rnk}
\alias{generate_rnk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate rnk format file
}
\description{
This function generates .rnk format files to use in the GSEA analysis
}
\usage{
generate_rnk(RNKfile, col.names, rnk)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RNKfile}{
Data frame with values of Symbol, P.Value and logFC
}
  \item{col.names}{
Vector with the name of the columns corresponding to Symbol, P.Value and logFC
}
  \item{rnk}{
Path and name of the output file with .rnk termination
}
}
\details{
.rnk files have a column with geneID and a column with a ranking score calculated using the formula -log10(p.Val)*sign(logFC)
}
\value{
.rnk format file
}
\references{
https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
}
\author{
Magdalena Arnal Segura <marnal@imim.es>
}
