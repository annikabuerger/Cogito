\name{MurEpi.RRBS.small}
\docType{data}
\alias{MurEpi.RRBS.small}
\title{Example data set: Murine methylation status data set of GSE77004}
\description{
This murine data from King et al. was downloaded from the NCBI GEO database 
under accession number GSE77004. 

The methylation status, mesured by RRBS, was similarly taken from the published 
files (accession number GSE84103), which contain the fraction of methylated 
cytosine for every CpG context supported by a minimum of 5 reads.

To reduce the storage size and complexity of the murine example data, the data 
only contains data of chr5 and four sample conditions (3x TFX and 1x mut) were 
removed.
}
\usage{MurEpi.RRBS.small}
\format{A GRanges Object containing a lot of Ranges with methylation status.}
\source{
NCBI GEO database accession number GSE84103 
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84103}
}
\references{
King AD, Huang K, Rubbi L, Liu S, Wang CY, Wang Y, Pellegrini M, Fan G. 
Reversible Regulation of Promoter and Enhancer Histone Landscape by DNA 
Methylation in Mouse Embryonic Stem Cells. Cell Rep. 2016 Sep 27;17(1):289-302. 
doi: 10.1016/j.celrep.2016.08.083
}
\keyword{datasets}