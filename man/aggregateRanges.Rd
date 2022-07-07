\name{aggregateRanges}
\docType{methods}
\alias{aggregateRanges}
\title{Aggregate GRanges with columns of attached values to genes}
\author{Annika Bürger}
\usage{
    aggregateRanges(ranges, configfile = NULL, organism = NULL,
                    referenceRanges = NULL, name = "", verbose = FALSE)
}
\arguments{
    \item{ranges}{
        list of GRanges, GRangesList or CompressedGRangesList
        with names in "RRBS|DNA|CNV|RNA|CHiP"}
    \item{configfile}{
        character, path to configuration file in json format}
    \item{organism}{TxDb or OrganismDb object, default value NULL}
    \item{referenceRanges}{list of GRanges, GRangesList or
        CompressedGRangesList with length one and name of reference,
        default value NULL}
    \item{name}{character, default value ""}
    \item{verbose}{logical, default value FALSE}
}
\value{
    List object with three members: One GRanges object with one gene
    or range of given reference per line and one column per sample,
    configuration information, and name.
}
\seealso{
    \code{\link{summarizeRanges}}
}
\description{
    Aggregates multiple GRanges objects with present columns of attached 
    values (mcols) to genes or ranges of given reference of given organism.
}
\examples{
mm9 <- TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene

### small artificial example ###
ranges.RNA.control <-
    GRanges(seq = "chr10",
            IRanges(c(41023369, 41211825, 41528287, 41994926, 42301673,
                        43256520, 43618919, 49503584, 51349066, 52099001),
                    c(41023544, 41212385, 41528663, 41995357, 42302290,
                        43257075, 43619492, 49504033, 51349425, 52099521)),
            seqinfo = GenomeInfoDb::seqinfo(mm9),
            expr = runif(5, 0, 1))
ranges.RNA.condition <-
    GRanges(seq = "chr10",
            IRanges(c(41013942, 41208731, 41535166, 41999999, 42292275, 
                        43256194, 43615562, 49497888, 51347046, 52092180),
                    c(41014274, 41209664, 41536039, 42000182, 42292965, 
                        43256430, 43615866, 49498362, 51347969, 52092733)),
            seqinfo = GenomeInfoDb::seqinfo(mm9),
            expr = runif(5, 0, 1))
ranges.ChIP.control <-
    GRanges(seq = "chr10",
            IRanges(c(41022835, 41307587, 42197924, 42302387, 42893825,
                        43259749, 43620352, 43721891, 44248812, 45207572,
                        49508713, 51309978, 51348779, 52101900, 52265513),
                    c(41022954, 41307745, 42198201, 42302555, 42893974,
                        43259889, 43620604, 43722051, 44248920, 45207704,
                        49508859, 51310187, 51348921, 52102030, 52265689)),
            seqinfo = GenomeInfoDb::seqinfo(mm9),
            score = round(runif(15, 5, 90)))

example.dataset <- list(RNA = GRangesList(control = ranges.RNA.control, 
                                            condition = ranges.RNA.condition), 
                        ChIP = ranges.ChIP.control)

aggregated.ranges <- aggregateRanges(ranges = example.dataset,
                                        organism = mm9, 
                                        name = "art.example", 
                                        verbose = TRUE)

names(aggregated.ranges)
head(aggregated.ranges$genes)
}
