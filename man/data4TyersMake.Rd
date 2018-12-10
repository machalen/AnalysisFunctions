\name{data4TyersMake}
\alias{data4TyersMake}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a data frame with the data4Tyers format from array data.
}
\description{
Generate a data frame with the data4Tyers format to run the macro in excel.
}
\usage{
data4TyersMake(est_noctrls.annot, cond, fit.main, contrast)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{est_noctrls.annot}{
Matrix with the intensity values of the microarray together with the annotation data made with the function annotatetableC.
}
  \item{cond}{
Vector with sample conditions used in the contrasts, sample order with the same order as the expression matrix.
}
  \item{fit.main}{
Model to use in limma analysis.
}
  \item{contrast}{
List of vectors with each contrast to use.
}

}
\details{
Build a data frame with the columns required to run the macro data4Tyers. This data frame contains information regarding gene anotation, statistics of the analysis and expression values of microarray data normalized with  Robust Multi-array Average (RMA).
}
\value{
data.frame object
}
\author{
Magdalena Arnal Segura <marnal@imim.es>
}
