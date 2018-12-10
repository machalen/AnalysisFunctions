\name{RNAseq.data4Tyers}
\alias{RNAseq.data4Tyers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a data frame with the data4Tyers format from RNAseq data.
}
\description{
Generate a data frame with the data4Tyers format to run the macro in excel.
}
\usage{
RNAseq.data4Tyers(expr.mat, annot.mat, cond, fit.main, contrast, specie="human", GOndKEGG=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr.mat}{
Matrix with expression values in CPM (E object provided by voom), TMM or counts with rownames equal to the "Geneid"" column in annot.mat.
}
  \item{annot.mat}{
Matrix with annotations obtained with featurecounts function in subread package it has to have a column with gene symbols named "Geneid".
}
  \item{cond}{
Vector with sample conditions used in the contrasts, sample order with same order as colnames(expr.mat).
}
  \item{fit.main}{
Model to use in limma analysis.
}
  \item{contrast}{
List of vectors with each contrast to use.
}
  \item{specie}{
Specie used to make annotations (human or mouse). It is required when using the option GOndKEGG=TRUE.
}
  \item{GOndKEGG}{
TRUE if we want to make KEGG and GO annotations for genes.
}
}
\details{
Build a data frame with the columns required to run the macro data4Tyers. This data frame contains information regarding gene anotation, statistics of the analysis and expression values (normalized or not) of RNAseq data.
}
\value{
data.frame object
}
\author{
Magdalena Arnal Segura <marnal@imim.es>
}
