\name{gct.data4Tyers}
\alias{gct.data4Tyers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gct file for GSEA analysis using data4Tyers
}
\description{
Function to generate .gct files from the data4Tyers object
}
\usage{
gct.data4Tyers(data4Tyers, Symb.col = "AffyID", gct)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data4Tyers}{
data4Tyers generated during the microarray analysis pipeline. It can be recovered using read.csv2
}
  \item{Symb.col}{
Column with gene IDs used in the GSEA analysis
}
  \item{gct}{
file path and file name of the output gct.
}
}
\details{
This function uses the template data4Tyers in order to make easy the generation of the gct files for the GSEA analysis.
}
\value{
gct format file
}

\references{
https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
}

\author{
Magdalena Arnal Segura <marnal@imim.es>
}
