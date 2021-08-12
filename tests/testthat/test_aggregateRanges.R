context("function Cogito::aggregateRanges")
library(Cogito)

test_that("function Cogito::aggregateRanges returns error on wrong parameter", 
          {
            expect_error(aggregateRanges("abc"))
})

test_that("function Cogito::aggregateRanges returns correct value", 
          {
            mm9 <- 
              TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
            ranges.RNA.control <-
              GRanges(seq = "chr10",
                      IRanges(c(41023369, 41211825, 41528287, 41994926, 
                                42301673, 43256520, 43618919, 49503584, 
                                51349066, 52099001),
                              c(41023544, 41212385, 41528663, 41995357, 
                                42302290, 43257075, 43619492, 49504033, 
                                51349425, 52099521)),
                      seqinfo = GenomeInfoDb::seqinfo(mm9),
                      expr = c(0.79, 0.11, 0.07, 0.34, 0.54))
            ranges.RNA.condition <-
              GRanges(seq = "chr10",
                      IRanges(c(41013942, 41208731, 41535166, 41999999, 
                                42292275, 43256194, 43615562, 49497888, 
                                51347046, 52092180),
                              c(41014274, 41209664, 41536039, 42000182, 
                                42292965, 43256430, 43615866, 49498362, 
                                51347969, 52092733)),
                      seqinfo = GenomeInfoDb::seqinfo(mm9),
                      expr = c(0.20, 0.65, 0.22, 0.45, 0.11))
            ranges.ChIP.control <-
              GRanges(seq = "chr10",
                      IRanges(c(41022835, 41307587, 42197924, 42302387, 
                                42893825, 43259749, 43620352, 43721891, 
                                44248812, 45207572, 49508713, 51309978, 
                                51348779, 52101900, 52265513),
                              c(41022954, 41307745, 42198201, 42302555, 
                                42893974, 43259889, 43620604, 43722051, 
                                44248920, 45207704, 49508859, 51310187, 
                                51348921, 52102030, 52265689)),
                      seqinfo = GenomeInfoDb::seqinfo(mm9),
                      score = c(24, 59, 17, 12, 29, 7, 45, 34, 28, 14, 58,
                               74, 24, 61, 32))
            
            example.dataset <- 
              list(RNA = GRangesList(control = ranges.RNA.control, 
                                     condition = ranges.RNA.condition), 
                   ChIP = ranges.ChIP.control)
            
            aggregated.ranges <- aggregateRanges(ranges = example.dataset,
                                                 organism = mm9, 
                                                 name = "art.example")
            
            expect_type(aggregated.ranges, "list")
            
            expect_equal(names(aggregated.ranges), 
                             c("genes", "config", "name"))
            
            expect_equal(aggregated.ranges$name, "art.example")
            
            expect_type(aggregated.ranges$config, "list")
            expect_equal(names(aggregated.ranges$config),
                         c("organism", "MaxDistToGene", 
                           "technologies", "conditions"))
            expect_equal(aggregated.ranges$config$organism, "mm9")
            expect_equal(aggregated.ranges$config$MaxDistToGene, 100000)
            expect_type(aggregated.ranges$config$technologies, "list")
            expect_equal(names(aggregated.ranges$config$technologies), 
                         c("RNA", "ChIP"))
            expect_equal(aggregated.ranges$config$technologies$RNA, 
                         c("RNA.control.expr", "RNA.condition.expr"))
            expect_equal(aggregated.ranges$config$technologies$ChIP, 
                         c("ChIP.score"))
            expect_type(aggregated.ranges$config$conditions, "list")
            expect_equal(length(aggregated.ranges$config$conditions), 0)
            
            expect_s4_class(aggregated.ranges$genes, "GRanges")
            expect_equal(length(aggregated.ranges$genes), 20)
            expect_equal(colnames(mcols(aggregated.ranges$genes)),
                         c("gene_id", "RNA.control.expr", 
                           "RNA.condition.expr", "ChIP.score"))
})
