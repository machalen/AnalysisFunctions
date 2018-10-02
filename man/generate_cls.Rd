\name{generate_cls}
\alias{generate_cls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate cls format files
}
\description{
This function generates .cls format files to use in the GSEA analysis.
}
\usage{
generate_cls(class, cls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{class}{
Sample conditions to be analyzed. Same order as samples in in the gct file
}
  \item{cls}{
Path and name of the output file with .cls termination
}
}
\details{
.cls files contain the information of the sample annotation to be used in the GSEA analysis contrasts.
}
\value{
.cls format file
}
\references{
https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
}
\author{
Magdalena Arnal Segura <marnal@imim.es>
}
