context("function Cogito::summarizeRanges")
library(Cogito)

test_that("function Cogito::summarizeRanges returns error on wrong parameter", 
          {
            expect_error(summarizeRanges("abc"))
})
