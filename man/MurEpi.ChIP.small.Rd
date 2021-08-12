\name{MurEpi.ChIP.small}
\docType{data}
\alias{MurEpi.ChIP.small}
\title{Example data set: Murine ChIP-seq data of GEO GSE77004}
\description{
This murine data from King et al. was downloaded from the NCBI GEO database 
under accession number GSE77004. 

The available ChIP-seq data (GSE77002) was then processed as described in King 
et al.: After alignment with bowtie with parameters selecting for uniquely 
mapped, best-matching reads and a maximum of two mismatches per read, the peak 
calling was done with homer findPeaks algorithm and an input control. 
Subsequently, the raw peaks were filtered with the following parameters: -F 8 
for H3K4me3, -size 1000 -minDist 3000 -F 4 -tagThreshold 32 for H3K27me3, -F 4 
for H3K27ac and -size 1000 -minDist1000 -nfr for H3K4me1.

To reduce the storage size and complexity of the murine example data, the data 
only contains data of chr5 and four sample conditions (3x TFX and 1x mut) were 
removed.
}
\usage{MurEpi.ChIP.small}
\format{A GRanges Object containing a lot of Ranges with scores.}
\source{
NCBI GEO database accession number GSE77002 
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77002}
}
\references{
King AD, Huang K, Rubbi L, Liu S, Wang CY, Wang Y, Pellegrini M, Fan G. 
Reversible Regulation of Promoter and Enhancer Histone Landscape by DNA 
Methylation in Mouse Embryonic Stem Cells. Cell Rep. 2016 Sep 27;17(1):289-302. 
doi: 10.1016/j.celrep.2016.08.083
}
\keyword{datasets}