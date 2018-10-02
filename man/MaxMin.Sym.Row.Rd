\name{MaxMin.Sym.Row}
\alias{MaxMin.Sym.Row}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter out duplicated genes getting the maximum or the minimum
}
\description{
When there are repeated Gene IDs in a matrix this function resolves those duplicated elements getting the gene value with the minimum or the maximum value in another column.
}
\usage{
MaxMin.Sym.Row(Annot.Table, ColField, operator = "Max", colName = "Symbol")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Annot.Table}{
Matrix with rows equal to geneId.
}
  \item{ColField}{
Name of the column used to chose the minimum or the maximum value of repeated gene Ids.
}
  \item{operator}{
Maximum or minimum operator to apply.
}
  \item{colName}{
Name of the column with symbol IDs.
}
}
\details{
When there are repeated Gene IDs in the column colName in a matrix this function resolves those duplicated elements getting the gene value with the minimum or the maximum value in the column ColField.
}
\value{
Matrix without repeated GeneIDs in the column colName.
}
\author{
Magdalena Arnal Segura <marnal@imim.es>
}
