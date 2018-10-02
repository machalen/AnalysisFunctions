\name{createTargets}
\alias{createTargets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Target file with sample anotation
}
\description{
This function creates the target file in tabulated (.txt) format with sample annotation to use in microarray and NGS analysis. The name of the samples in ordered as the one in the R objects to avoid possible mistakes.
}
\usage{
createTargets(names, fileName = "Targets.txt", DestDir = analysisDir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{names}{
Vector with sample names as they appear in the R matrix to analyze. Same name and same order.
}
  \item{fileName}{
Name of the output file, default is Targets.txt
}
  \item{DestDir}{
Name of the output directory
}
}
\details{
The output file has to be manually filled in a posteriori with the clinical data or sample information of interest.
}
\value{
Tabulated separated (.txt) file 
}

\author{
Magdalena Arnal Segura <marnal@imim.es>
}

