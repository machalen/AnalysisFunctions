\name{GSEA_report_HS}
\alias{GSEA_report_HS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to make a GSEA summary table with human gene symbol annotation
}
\description{
In order fo the results of GSEA to be more comprehensive, we make a function to summarize the results
}
\usage{
GSEA_report_HS(GSEA_result, GeneSetAnnot=F ,GeneSets_annot=NULL, TabCols=NULL, 
                             GeneStats_annot=NULL, StatsCols=NULL, GSEA_Dir)
}

\arguments{
  \item{GSEA_result}{
Directory path with results of one GSEA analysis
}
  \item{GeneSetAnnot}{
TRUE When gene sets are described in a table and we want them to be included
}
  \item{GeneSets_annot}{
Table with gene sets annotation and rownames(GeneSets_annot) == gene set name
}
  \item{TabCols}{
Columns of the GeneSets_annot to be included
}
  \item{GeneStats_annot}{
Table with the statistics of genes and rownames == Genesymbol
}
  \item{StatsCols}{
Columns of the table GeneStats_annot to be included
}
  \item{GSEA_Dir}{
Directory to save the results
}
}
\details{
This function makes a summary of the GSEA results with the different sub-folders. The output .csv contains all gene sets in the top (selected by the user) together with the genes in each case.
}
\value{
data.frame ofject
}
\examples{
for (i in 1:length(inDirFiles)) {
  GSEA_result <- file.path(GSEA_Dir, inDirFiles[i])
  GSEA_report_HS(GSEA_result=GSEA_result, GeneSetAnnot=T ,GeneSets_annot=GSEA_info, 
                   TabCols=c("source","descriptionBrief","descriptionFull", "nGenes"))
  
}
}
