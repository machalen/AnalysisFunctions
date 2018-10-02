\name{DeleteCtrlsGral}
\alias{DeleteCtrlsGral}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Delete probeID that correspond to controls
}
\description{
This function filters out rows in a microarray expression matrix that correspond to control ProbeIDs. 
}
\usage{
DeleteCtrlsGral(estimates, AnnotfilePath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estimates}{
Numeric matrix with microarray estimates. Rownames are the transcript cluster IDs.
}
  \item{AnnotfilePath}{
File path were annotations are stored. Files have to be in a .csv format provided by Affymetrix.
}
}
\details{
Control probe IDs do not contain information about biological samples, therefore can be deleted from de expression matrix. To do that, the annotation files provided by Affymetrix in .csv format are used.
}
\value{
Expression matrix without control ProbeIDs.
}

\author{
Magdalena Arnal Segura <marnal@imim.es>
}
