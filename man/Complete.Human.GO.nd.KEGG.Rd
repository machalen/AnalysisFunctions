\name{Complete.Human.GO.nd.KEGG}
\alias{Complete.Human.GO.nd.KEGG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functional annotation of RNA-seq data with GO and KEGG databases.
}
\description{
Functional annotation of RNA-seq data with GO and KEGG databases for human and human.
}
\usage{
Complete.Human.GO.nd.KEGG(annot.mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annot.mat}{
Matrix with a column "Geneid" corresponding to gene symbols human.
}
}
\details{
Functional annotation of RNAseq results using the GO and KEGG databases. Annotation is made with the r ‘GO.db’ ‘KEGG.db’ and ‘org.Hs.eg.db’ packages and they need to be installed.
}
\value{
data.frame object
}
\author{
Magdalena Arnal Segura <marnal@imim.es>
}
