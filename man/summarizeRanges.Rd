\name{summarizeRanges}
\docType{methods}
\alias{summarizeRanges}
\title{Summarize Aggregated GRanges}
\author{Annika Bürger}
\usage{
    summarizeRanges(aggregated.ranges, outputFormat = "pdf", verbose = FALSE)
}
\arguments{
    \item{aggregated.ranges}{list of GRanges, cofiguration information, and name
        for example result from function Cogito::summarizeRanges}
    \item{outputFormat}{character, can be pdf or html default value pdf}
    \item{verbose}{logical, default value FALSE}
}
\value{
    No return value, only side effects: creation of a rmd, a pdf or html and a
    data file (RData).
}
\seealso{
    \code{\link{aggregateRanges}}
}
\description{
    Summarize GRanges with present columns of attached values (mcols).
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

# adding information about conditions
aggregated.ranges$config$conditions <- list(condition = c("RNA.condition.expr"),
                                            control = c("RNA.control.expr",
                                                        "ChIP.score"))
summarizeRanges(aggregated.ranges = aggregated.ranges,
                outputFormat = "pdf",
                verbose = TRUE)
}
