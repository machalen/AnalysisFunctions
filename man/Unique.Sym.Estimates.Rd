\name{Unique.Sym.Estimates}
\alias{Unique.Sym.Estimates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Resolve duplicated gene IDs by computing the mean or the median
}
\description{
This function removes duplicated gene IDs by calculating the mean or the median of their expression values
}
\usage{
Unique.Sym.Estimates(Annot.Table, operator = "Mean", colName = "Symbol")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Annot.Table}{
Numeric data frame with a unique column of characters corresponding to the gene IDs, the other columns are samples.
}
  \item{operator}{
Mean or median operator to resolve repeated gene IDs 
}
  \item{colName}{
Name of the column that contains gene symbols.
}
}
\details{
A gene expression matrix can have repeated gene IDs due to exon or probe annotation. Some methods need to be unique gene ID based. To summarize expression values to a single geneID the mean and median operators can be applied. 
}
\value{
Numeric data frame where row names are gene IDs 
}
\author{
Magdalena Arnal Segura <marnal@imim.es>
}
