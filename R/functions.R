# ranges have to be list of GRanges, GRangesList or CompressedGRangesList
# names of ranges have to be in "RRBS|DNA|CVN|RNA|ChIP"
# organism has to be a TxDb or OrganismDb
aggregateRanges <- function(ranges, configfile = NULL, 
                            organism = NULL,
                            name = "", verbose = FALSE) {
    if (!is.list(ranges)) {
        stop("Cogito::aggregateRanges parameter ranges is no list")
    }
    w <- unlist(lapply(ranges, function(range) {
        return(
            is.list(range) || is(range, "GRanges") || is(range, "GRangesList")
            )
    }))
    ranges <- ranges[w]
    w <- which(unlist(lapply(ranges, is.list)))
    for (i in w) {
        ranges[[i]] <-
            ranges[[i]][unlist(lapply(ranges[[i]], is, "GRanges"))]
        if (length(ranges[[i]]) == 0) {
            ranges[[i]] <- NULL
        }
    }
    ranges <-
        ranges[grep("RRBS|DNA|CNV|RNA|ChIP", names(ranges), ignore.case = TRUE)]
    if (length(ranges) == 0) {
        stop(
            "Cogito::aggregateRanges ",
            "members of parameter ranges has to be GRanges, GRangesList, ",
            "CompressedGRangesList or list of GRanges and the given base ",
            "technologies (names of the parameter ranges) have to be one of ",
            "RRBS, DNA, CNV, RNA, or ChIP"
        )
    }

    if (is.null(organism)) {
        organism <- 
            TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
    }
    if (!is(organism, "TxDb") && !is(organism, "OrganismDb")) {
        stop(
            "Cogito::aggregateRanges ",
            "wrong class of argument organism, only TxDb or ",
            "OrganismDb are allowed"
        )
    }
    genome <- organism
    
    if (is.na(verbose) | !is.logical(verbose)) {
        verbose <- FALSE
    }
    if (verbose) {
        message("start function Cogito::aggregateRanges")
    }
    ## load configlist
    configlistauto <- list(organism = names(summary(as.factor(
        genome(genome)
    )))[1], MaxDistToGene = 100000)
    if (is.null(configfile) || !is.character(configfile)) {
        configfile <- paste0(getwd(), "/", Sys.Date(), "_", name, "_config.txt")
        configlist <- configlistauto
    } else {
        configlist <- jsonlite::fromJSON(paste(readLines(
            configfile,
            warn = FALSE
        ), collapse = ""))
        configlist$organism <- configlistauto$organism
        if (is.null(configlist$MaxDistToGene) ||
            !is.numeric(configlist$MaxDistToGene) ||
            configlist$MaxDistToGene < 0) {
            configlist$MaxDistToGene <- configlistauto$MaxDistToGene
        }
    }

    if ("SYMBOL" %in% AnnotationDbi::columns(genome)) {
        genes <- sort(suppressMessages(
            GenomicFeatures::genes(genome, columns = c("SYMBOL"))
        ))
    } else {
        genes <- sort(suppressMessages(
            GenomicFeatures::genes(genome)
        ))
    }

    ranges <- lapply(ranges, function(range) {
        if (methods::is(range, "GRanges")) {
            return(list(range))
        } else {
            return(range)
        }
    })

    maxdist <- configlist$MaxDistToGene

    ## MEAN: distribute RRBS/RNA/DNA/CNV over genes
    ## (numeric: mean, ordered: median, factor: modus, logical: any)
    l <- ranges[grep("RRBS|RNA|DNA|CNV", names(ranges), ignore.case = TRUE)]
    for (ll in seq_along(l)) {
        if (!is.null(l[[ll]])) {
            if (verbose) {
                message("- annotating with ", names(l)[[ll]])
            }
            for (i in seq_along(l[[ll]])) {
                if (verbose && length(l[[ll]]) > 1) {
                    Cgt.progressbar(i, seq_along(l[[ll]]))
                }
                tmp <- l[[ll]][[i]]
                colnames(mcols(tmp)) <-
                    paste0(
                        names(l)[ll], ".", names(l[[ll]])[i],
                        ifelse(is.null(names(l[[ll]])[i]), "",
                            "."
                        ), colnames(mcols(tmp))
                    )
                if (all(unlist(lapply(mcols(tmp), is.numeric)))) {
                    genes <- Cgt.aggr(genes, tmp,
                        scale = "numeric",
                        mode = "mean", maxdist
                    )
                } else if (all(unlist(lapply(mcols(tmp), is.ordered)))) {
                    genes <- Cgt.aggr(genes, tmp,
                        scale = "ordered",
                        mode = "median", maxdist
                    )
                } else if (all(unlist(lapply(mcols(tmp), is.factor)))) {
                    genes <- Cgt.aggr(genes, tmp,
                        scale = "factor",
                        mode = "modus", maxdist
                    )
                } else if (all(unlist(lapply(mcols(tmp), is.logical)))) {
                    genes <- Cgt.aggr(genes, tmp,
                        scale = "logical",
                        mode = "any", maxdist
                    )
                } else {
                    stop("Cogito::aggregateRanges no uniform scaling")
                }
            }
        }
    }
    ## MAX: distribute ChIP over genes
    ## (numeric: max, ordered: max, factor/logical: modus)
    l <- ranges[grep("ChIP", names(ranges), ignore.case = TRUE)]
    for (ll in seq_along(l)) {
        if (!is.null(l[[ll]])) {
            if (verbose) {
                message("- annotating with ", names(l)[[ll]])
            }
            for (i in seq_along(l[[ll]])) {
                if (verbose && length(l[[ll]]) > 1) {
                    Cgt.progressbar(i, seq_along(l[[ll]]))
                }
                tmp <- l[[ll]][[i]]
                colnames(mcols(tmp)) <-
                    paste0(
                        names(l)[ll], ".", names(l[[ll]])[i],
                        ifelse(is.null(names(l[[ll]])[i]), "",
                            "."
                        ), colnames(mcols(tmp))
                    )
                if (all(unlist(lapply(mcols(tmp), is.numeric)))) {
                    genes <- Cgt.aggr(genes, tmp,
                        scale = "numeric",
                        mode = "max", maxdist
                    )
                } else if (all(unlist(lapply(mcols(tmp), is.ordered)))) {
                    genes <- Cgt.aggr(genes, tmp,
                        scale = "ordered",
                        mode = "max", maxdist
                    )
                } else if (all(unlist(lapply(mcols(tmp), is.factor)))) {
                    genes <- Cgt.aggr(genes, tmp,
                        scale = "factor",
                        mode = "modus", maxdist
                    )
                } else if (all(unlist(lapply(mcols(tmp), is.logical)))) {
                    genes <- Cgt.aggr(genes, tmp,
                        scale = "logical",
                        mode = "modus", maxdist
                    )
                } else {
                    stop("Cogito::aggregateRanges no uniform scaling")
                }
            }
        }
    }

    configlist$technologies <- lapply(names(ranges), function(tech) {
        colnames(mcols(genes))[
            grep(paste0(tech, "."), colnames(mcols(genes)), fixed = TRUE)
        ]
    })
    names(configlist$technologies) <- names(ranges)

    # due to automated typecasts from list to matrix in edge cases
    if (is.matrix(configlist$technologies)) {
        configlist$technologies <- lapply(as.list(data.frame(
            configlist$technologies
        )), as.character)
    }

    # guess conditions from names
    shortcolnames <- unlist(lapply(strsplit(colnames(mcols(genes))[-1],
                                            ".", fixed=TRUE), "[[", 2))
    shorttechs <- lapply(configlist$technologies, function(tech) {
        unlist(lapply(strsplit(tech, ".", fixed=TRUE), "[[", 2))})
    groupsIdxs <- Cgt.guessConditions(shortcolnames, shorttechs)
    configlist$conditions <- lapply(groupsIdxs, function(group) {
        colnames(mcols(genes))[-1][group]})

    write(jsonlite::toJSON(configlist, pretty = TRUE),
        append = FALSE,
        file = configfile
    )

    genes <- genes[!apply(mcols(genes), 1, function(row) all(is.na(row[-1]))), ]
    return(list(genes = genes, config = configlist, name = name))
}

summarizeRanges <- function(aggregated.ranges, verbose = FALSE) {
    if (!is.list(aggregated.ranges) ||
        !("genes" %in% names(aggregated.ranges)) ||
        !methods::is(aggregated.ranges$genes, "GRanges")) {
        stop("Cogito::summarizeRanges wrong format")
    }
    if (length(aggregated.ranges$genes) == 0) {
        stop("Cogito::summarizeRanges parameter ranges is empty")
    }
    if (is.na(verbose) || !is.logical(verbose)) {
        verbose <- FALSE
    }
    if (verbose) {
        message("start function Cogito::summarizeRanges")
    }
    ## preparations ####
    dat <- Cgt.getDataFromGRanges(aggregated.ranges$genes)
    scales <- dat$scales
    dat <- dat$dat
    config <- aggregated.ranges$config

    if (is.null(aggregated.ranges$name) ||
        !is.character(aggregated.ranges$name)) {
        name <- ""
    } else {
        name <- aggregated.ranges$name
    }

    ## 1. summary of single columns of attached values ####
    rmdSummAnn <- vector("list", length(scales))
    names(rmdSummAnn) <- names(dat)
    for (i in seq_along(scales)) {
        if (verbose) {
            Cgt.progressbar(i, seq_along(scales), "summary of attached values")
        }
        switch(scales[i],
            bin = {
                summar <-
                    Cgt.plotAttributeNominal(
                        as.list(dat)[i],
                        relativScale = FALSE,
                        sort = TRUE, na.rm = TRUE,
                        usedLevelsOnly = FALSE,
                        main = paste(name, names(dat)[i])
                    )
            },
            nom = {
                summar <-
                    Cgt.plotAttributeNominal(
                        as.list(dat)[i],
                        relativScale = FALSE,
                        sort = TRUE, na.rm = TRUE,
                        usedLevelsOnly = TRUE,
                        main = paste(name, names(dat)[i])
                    )
                locname <- names(summar$location)
                summar$location <-
                    paste0(
                        summar$location,
                        " (", summar$data[1], " = ",
                        round(summar$data[1] / sum(summar$data) * 100), "%)"
                    )
                if (sum(summar$data[1] == summar$data) > 1) {
                    summar$location <-
                        paste(
                            summar$location, "among",
                            sum(summar$data == summar$data[1]),
                            "levels with equal frequency."
                        )
                }
                names(summar$location) <- locname
            },
            ord = {
                summar <-
                    Cgt.plotAttributeOrdinal(
                        as.list(dat)[i],
                        relativScale = FALSE,
                        usedLevelsOnly = FALSE, na.rm = TRUE,
                        main = paste(name, names(dat)[i])
                    )
            },
            int = ,
            rat = {
                summar <- Cgt.plotAttributeCont(
                    as.list(dat)[i],
                    main = paste(name, names(dat)[i]),
                    scale = scales[i]
                )
                summar$location <- Cgt.roundNicely(summar$location)
                summar$dispersion <- Cgt.roundNicely(summar$dispersion)
            }
        )
        rmdSummAnn[[i]] <- summar
    }
    objectstosave <- "rmdSummAnn"
    
    ## 2. summary of columns of attached values ####
    ## of same base technology and condition (=groups)
    if (!is.null(config$conditions) && length(config$conditions) > 0 &&
        !is.null(config$technologies) && length(config$technologies) > 0) {
        groups <- lapply(
            seq_len(length(config$technologies) * length(config$conditions)),
            function(group) {
                intersect(
                    config$technologies[[
                    (group - 1) %% length(config$technologies) + 1]],
                    config$conditions[[
                    (group - 1) %/% length(config$technologies) + 1]]
                )
            })
        names(groups) <- as.character(outer(
            paste0(names(config$technologies), "."),
            names(config$conditions), paste0
        ))
        groups <- lapply(
            lapply(groups, match, colnames(dat)), function(x) x[!is.na(x)]
        )
        groups <- groups[lapply(groups, length) > 1]
        rmdSummGroups <- vector("list", length(groups))
        names(rmdSummGroups) <- names(groups)
        for (i in seq_along(groups)) {
            if (verbose) {
                Cgt.progressbar(
                    i, seq_along(groups),
                    paste(
                        "summary of attached values of tracks with the",
                        "same base technology and condition"
                    )
                )
            }
            summar <- list(idx = groups[[i]])
            switch(scales[groups[[i]][1]],
                rat = ,
                int = {
                    summar$dat <-
                        Cgt.plotAttributeCont(as.list(dat)[
                            groups[[i]]
                        ],
                        main = names(groups)[i],
                        scale = scales[groups[[i]][1]],
                        groups = config$conditions
                        )
                },
                bin = {
                    summar$dat <-
                        Cgt.plotAttributeNominal(as.list(dat)[
                            groups[[i]]
                        ],
                        relativScale = FALSE,
                        usedLevelsOnly = TRUE,
                        main = names(groups)[i],
                        na.rm = TRUE
                        )
                },
                nom = {
                    summar$dat <-
                        Cgt.plotAttributeNominal(as.list(dat)[
                            groups[[i]]
                        ],
                        relativScale = FALSE,
                        sort = TRUE,
                        usedLevelsOnly = TRUE,
                        main = names(groups)[i],
                        na.rm = TRUE
                        )
                },
                ord = {
                    summar$dat <-
                        Cgt.plotAttributeOrdinal(as.list(dat)[
                            groups[[i]]
                        ],
                        relativScale = FALSE,
                        usedLevelsOnly = FALSE,
                        main = names(groups)[i],
                        na.rm = TRUE
                        )
                }
            )
            rmdSummGroups[[i]] <- summar
        }
        objectstosave <- c(objectstosave, "rmdSummGroups")
    }

    ## 3. summary of columns of attached values of same base technology ####
    if (!is.null(config$technologies) && length(config$technologies) > 0) {
        techs <- config$technologies
        techs <- lapply(
            lapply(techs, match, colnames(dat)),
            function(x) x[!is.na(x)]
        )
        techs <- techs[lapply(techs, length) != 0]
        rmdSummTechs <- vector("list", length(techs))
        names(rmdSummTechs) <- names(techs)
        for (i in seq_along(techs)) {
            if (verbose) {
                Cgt.progressbar(
                    i, seq_along(techs),
                    paste(
                        "summary of attached values of tracks",
                        "with the same technology"
                    )
                )
            }
            summar <- list(idx = techs[[i]])
            if (length(techs[[i]]) == 1) { # fun with type conversions
                argAtts <- list(dat[, techs[[i]]])
            } else {
                argAtts <- as.list(dat[, techs[[i]]])
            }
            switch(scales[techs[[i]][1]],
                "rat" = ,
                "int" = {
                    summar$dat <-
                        Cgt.plotAttributeCont(
                            argAtts,
                            main = names(techs)[i],
                            scale = scales[techs[[i]][1]],
                            groups = config$conditions
                        )
                },
                "bin" = {
                    summar$dat <-
                        Cgt.plotAttributeNominal(
                            argAtts,
                            relativScale = FALSE,
                            usedLevelsOnly = TRUE,
                            main = names(techs)[i], na.rm = TRUE,
                            groups = config$conditions
                        )
                },
                "nom" = {
                    summar$dat <-
                        Cgt.plotAttributeNominal(
                            argAtts,
                            relativScale = FALSE, sort = TRUE, na.rm = TRUE,
                            usedLevelsOnly = TRUE, main = names(techs)[i],
                            groups = config$conditions
                        )
                },
                "ord" = {
                    summar$dat <-
                        Cgt.plotAttributeOrdinal(
                            argAtts,
                            relativScale = FALSE, usedLevelsOnly = FALSE,
                            main = names(techs)[i], na.rm = TRUE,
                            groups = config$conditions
                        )
                }
            )
            rmdSummTechs[[i]] <- summar
        }
        objectstosave <- c(objectstosave, "rmdSummTechs")
    }

    ## 4. comparison of columns of attached values ####
    pairs <- t(combn(seq_along(scales), 2))
    rmdCompAnn <- vector("list", nrow(pairs))
    names(rmdCompAnn) <-
        apply(pairs, 1, function(pair)
            paste0(names(dat)[pair], collapse = "_vs_"))
    for (pair in seq_len(nrow(pairs))) {
        if (verbose) {
            Cgt.progressbar(pair, seq_len(nrow(pairs)), 
                            paste(
                                "pairwise comparisons of",
                                "attached values from different tracks"
                                )
                            )
        }
        rmdCompAnn[[pair]] <-
            do.call(
                paste0(
                    "Cgt.plotComparison",
                    c("CatToCat", "CatToCont", "ContToCont")[
                        ((scales %in% c("int", "rat"))[pairs[, 1]] +
                            (scales %in% c("int", "rat"))[pairs[, 2]])[pair] + 1
                    ]
                ),
                list(
                    datComp = as.list(dat[, pairs[pair, ]]),
                    scale1 = scales[pairs[pair, 1]],
                    scale2 = scales[pairs[pair, 2]]
                )
            )
        # adjust p-value with bonferroni
        rmdCompAnn[[pair]]$p.value <- nrow(pairs) * rmdCompAnn[[pair]]$p.value
    }
    if (verbose) {
        message("save data")
    }
    objectstosave <- c(objectstosave, "rmdCompAnn")

    ## save data ####
    datasetfile <- paste0(getwd(), "/", Sys.Date(), "_", name, "_dataset.RData")
    save(file = datasetfile, list = objectstosave, compress = TRUE)

    ## write header in rmd-file for html or PDF ####
    if (verbose) {
        message("write report")
    }
    outputFormat <- "PDF"
    rmdFile <- paste0(getwd(), "/", Sys.Date(), "_", name, "_custom.rmd")
    if (outputFormat == "PDF") {
        write(paste0(
            "---\ntitle: \"", name, " - Cogito report ",
            Sys.Date(), "\"\noutput: pdf_document\n",
            "fontsize: 10pt\nlinkcolor: blue\n---\n"
        ),
        file = rmdFile
        )
    } else {
        write(paste0(
            "---\ntitle: \"", name, " - Cogito report ",
            Sys.Date(), "\"\noutput:\nhtml_document:\n",
            "css: style.css\n---\n"
        ), file = rmdFile)
    }

    write(Cgt.rcodestring(
        paste0("knitr::opts_chunk$set(echo = FALSE, out.width='60%', ",
                "fig.align=\"center\")\n",
                "library(ggplot2)\nload(\"", datasetfile, "\")"),
        "set chunk options and load libraries and data", echo=FALSE
    ),
    append = TRUE, file = rmdFile
    )

    write(paste0(
        "The given ", 
        paste(unique(GenomeInfoDb::genome(aggregated.ranges$genes)),
                collapse = ", "
        ), 
        " dataset *", name ,"* consists of ", 
        ncol(mcols(aggregated.ranges$genes)) - 1,
        " tracks with different base technologies and different ",
        "conditions as shown in the following table.\n"
    ),
    append = TRUE, file = rmdFile
    )

    if (!is.null(config$conditions) && length(config$conditions) > 0 &&
        !is.null(config$technologies) && length(config$technologies) > 0) {
        overviewMat <- matrix(NA, nrow=length(config$conditions),
                                ncol=length(config$technologies))
        rownames(overviewMat) <- names(config$conditions)
        colnames(overviewMat) <- names(config$technologies)
        for (row in seq_len(nrow(overviewMat)))
            for (col in seq_len(ncol(overviewMat)))
                overviewMat[row, col] <-
            sum(config$conditions[[row]] %in% config$technologies[[col]])
        write(Cgt.tablestring(overviewMat), append = TRUE, file = rmdFile)
    }

    write(paste0(
        "\n### Table of content\n",
        "1. ", Cgt.intlinkstring(
            "Summary of attached values of single tracks",
            "summaryann"
        ), "\n",
        "2. ", Cgt.intlinkstring(
            paste(
                "Summary of attached values of tracks of the same condition",
                "and technology"
            ), "summarygroupsann"
        ), "\n",
        "3. ", Cgt.intlinkstring(
            "Summary of attached values of tracks of the same technology", 
            "summarytechsann"
        ), "\n",
        "4. ", Cgt.intlinkstring(
            "Comparisons of attached values of different tracks", 
            "summarycomp"
        )
    ),
    append = TRUE, file = rmdFile
    )

    ## write 1. in rmd-file for html or PDF ####
    write("\n\\newpage", append = TRUE, file = rmdFile)
    write("\n# 1. Summary of attached values of single tracks{#summaryann}",
        append = TRUE, file = rmdFile
    )
    tmptext <- paste(
        "In the following table each row summarizes one column of attached",
        "values of one track of the given dataset. The scale",
        "is one of binary, nominal, ordinal, interval and rational. To get an",
        "overview about the attached values, a location and a deviation",
        "parameter is shown in the third and fourth column.",
        "Which location and deviation parameter is used depends on the scale",
        "of the attached values and is explained in the particular subsection."
    )
    write(paste("\n ", tmptext), append = TRUE, file = rmdFile)

    write("\n", append = TRUE, file = rmdFile)
    summtable <- data.frame(
        "." = substr(
            names(dat), 1,
            unlist(regexec(".", names(dat), fixed = TRUE)) - 1
        ),
        "name" = Cgt.intlinkstring(
            substr(
                names(dat),
                unlist(regexec(".", names(dat), fixed = TRUE)) + 1,
                max(width(names(dat)))
            ), names(dat)
        ),
        scale = c("rational", "interval", "ordinal", "nominal", "binary")[
            as.integer(factor(scales,
                levels = c("rat", "int", "ord", "nom", "bin")
            ))
        ],
        mean = unlist(lapply(rmdSummAnn, function(sAnn) sAnn$location)),
        dispersion = unlist(lapply(rmdSummAnn, function(sAnn) sAnn$dispersion)),
        row.names = NULL
    )
    summtable$'.'[duplicated(summtable$'.')] <- "."
    write(
        Cgt.tablestring(
            summtable,
            row.names = FALSE,
            col.widths =
                c(
                    max(width(names(dat))) -
                        max(unlist(regexec(".", names(dat), fixed = TRUE))),
                    8, 6, 6
                )
        ),
        append = TRUE, file = rmdFile
    )
    write("\n", append = TRUE, file = rmdFile)
    for (i in seq_along(rmdSummAnn)) {
        if (summtable[i,'.'] != ".") {
            write("\\newpage", append = TRUE, file = rmdFile)
            write(
                paste0("\n## ", summtable[i,'.']),
                append = TRUE, file = rmdFile
            )
        }
        write(paste0("\n### ", names(dat)[i], "{#", names(dat)[i], "}"),
            append = TRUE, file = rmdFile
        )
        write(paste0(
            "\nLocation parameter (",
            names(rmdSummAnn[[i]]$location), "): ",
            rmdSummAnn[[i]]$location
        ),
        append = TRUE, file = rmdFile
        )
        write(paste0(
            "\nDispersion parameter (",
            names(rmdSummAnn[[i]]$dispersion), "): ",
            rmdSummAnn[[i]]$dispersion
        ),
        append = TRUE, file = rmdFile
        )
        write(Cgt.rcodestring(
            paste0(
                "plotdata <- rmdSummAnn[[", i,
                "]]$ggplotdata\nggplot(plotdata",
                rmdSummAnn[[i]]$ggplotstring
            ),
            paste0("Step1-", i, "-", length(rmdSummAnn))
        ),
        append = TRUE, file = rmdFile
        )
        write(paste0("\n", Cgt.intlinkstring("back to overview", "summaryann")),
            append = TRUE, file = rmdFile
        )
    }

    ## write 2. in rmd-file for html or PDF ####
    # Summary of attached values of same condition and 
    # base technology (=groups)
    write("\n\\newpage", append = TRUE, file = rmdFile)
    write(paste(
        "\n# 2. Summary of attached values of tracks of the same condition",
        "and technology{#summarygroupsann}"
    ),
    append = TRUE, file = rmdFile
    )

    if (exists("rmdSummGroups") && length(rmdSummGroups) > 0) {
        write("\n", append = TRUE, file = rmdFile)
        tmptext <-
            paste("*", Cgt.intlinkstring(
                names(rmdSummGroups),
                paste0("group", names(rmdSummGroups))
                ),
                "\n\t+", unlist(lapply(
                    lapply(groups, function(x) names(dat)[x]),
                    function(y) paste(y, collapse = "\n\t+ "))),
                collapse = "\n")
        write(tmptext, append = TRUE, file = rmdFile)

        for (i in seq_along(rmdSummGroups)) {
            write("\n\\newpage", append = TRUE, file = rmdFile)
            write(paste0(
                "\n## ", names(rmdSummGroups)[i],
                "{#group", names(rmdSummGroups)[i], "}"
            ),
            append = TRUE, file = rmdFile
            )
            write(Cgt.rcodestring(
                paste0(
                    "plotdata <- rmdSummGroups[[", i,
                    "]]$dat$ggplotdata\nggplot(plotdata",
                    rmdSummGroups[[i]]$dat$ggplotstring
                ),
                paste0("Step2-", i, "-", length(rmdSummGroups))
            ),
            append = TRUE, file = rmdFile
            )
            write(paste0("\n", Cgt.intlinkstring(
                "back to overview",
                "summarygroupsann"
            )),
            append = TRUE, file = rmdFile
            )
        }
    } else {
        write(paste0(
            "\n", "No such groups of attached values with equal ",
            "scales are present in the given dataset."
            ),
            append = TRUE, file = rmdFile
        )
    }

    ## write 3. in rmd-file for html or PDF ####
    # Summary of attached values of same technology
    write("\n\\newpage", append = TRUE, file = rmdFile)
    write(paste(
        "\n# 3. Summary of attached values of tracks of the same",
        "technology{#summarytechsann}"
    ), append = TRUE, file = rmdFile)

    if (exists("rmdSummTechs") && length(rmdSummTechs) > 0) {
        write("\n", append = TRUE, file = rmdFile)
        # new
        tmptext <-
            paste("*",
                Cgt.intlinkstring(
                    names(rmdSummTechs),
                    paste0("technology", names(rmdSummTechs))
                ),
                "\n\t+", unlist(lapply(lapply(techs, function(x) {
                    names(dat)[x]
                }), function(y) paste(y, collapse = "\n\t+ "))),
                collapse = "\n"
            )
        write(tmptext, append = TRUE, file = rmdFile)

        for (i in seq_along(rmdSummTechs)) {
            write("\n\\newpage", append = TRUE, file = rmdFile)
            write(paste0(
                "\n## ", names(rmdSummTechs)[i],
                "{#technology", names(rmdSummTechs)[i], "}"
            ),
            append = TRUE, file = rmdFile
            )
            write(Cgt.rcodestring(
                paste0(
                    "plotdata <- rmdSummTechs[[", i,
                    "]]$dat$ggplotdata\nggplot(plotdata",
                    rmdSummTechs[[i]]$dat$ggplotstring
                ),
                paste0("Step3-", i, "-", length(rmdSummTechs))
            ),
            append = TRUE, file = rmdFile
            )
            write(paste0("\n", Cgt.intlinkstring(
                "back to overview",
                "summarytechsann"
            )),
            append = TRUE, file = rmdFile
            )
        }
    } else {
        write(paste0("\n", "No technology with equal scale present."),
            append = TRUE, file = rmdFile
        )
    }

    ## write 4. in rmd-file for html or PDF ####
    write("\n\\newpage", append = TRUE, file = rmdFile)
    write(paste(
        "\n# 4. Comparisons of attached values of different tracks",
        "{#summarycomp}"
    ), append = TRUE, file = rmdFile)
    write(paste(
        "Each pair of attached values has been tested",
        "for correlation. \n"
    ), append = TRUE, file = rmdFile)

    # matrix with significance-values
    # (and internal links and hover text)
    # adjusted p-values of correlation and significance
    corMat <- matrix("", ncol(dat) - 1, ncol(dat) - 1)
    corMatPvals <- matrix(0, ncol(dat) - 1, ncol(dat) - 1)
    for (pair in seq_len(nrow(pairs))) {
        corMat[pairs[pair, 1], ncol(dat) - pairs[pair, 2] + 1] <-
            Cgt.tooltipstring(Cgt.intlinkstring(
                ifelse(is.na(rmdCompAnn[[pair]]$p.value) ||
                            rmdCompAnn[[pair]]$p.value > 0.05 ||
                    (!is.na(rmdCompAnn[[pair]]$correlation) &&
                            rmdCompAnn[[pair]]$correlation < 0.5),
                "ns",
                ifelse(rmdCompAnn[[pair]]$p.value <= 0.001,
                    "***",
                    ifelse(rmdCompAnn[[pair]]$p.value <= 0.01, "**", "*")
                )
                ),
                paste0("comp", pair)
            ),
            paste(
                ifelse(is.na(rmdCompAnn[[pair]]$p.value) ||
                            rmdCompAnn[[pair]]$p.value > 0.05,
                    "not significant",
                    Cgt.roundNicely(rmdCompAnn[[pair]]$p.value)
                ),
                names(rmdCompAnn[[pair]]$p.value)
            ),
            format = outputFormat
            )
        corMatPvals[pairs[pair, 1], ncol(dat) - pairs[pair, 2] + 1] <-
            rmdCompAnn[[pair]]$p.value
    }
    colnames(corMatPvals) <- rev(names(dat)[-1])
    rownames(corMatPvals) <- names(dat)[-ncol(dat)]

    ## heatmap
    write(Cgt.rcodestring(
        paste0(
            "ggplotdata <- data.frame(\n",
            "  val1=c(unlist(lapply(strsplit(",
            "names(rmdCompAnn), \"_vs_\"), \"[\",",
            " 1)), unlist(lapply(strsplit(names(",
            "rmdCompAnn), \"_vs_\"), \"[\", 2))),\n",
            "  val2=c(unlist(lapply(strsplit(names(",
            "rmdCompAnn), \"_vs_\"), \"[\", 2)), ",
            "unlist(lapply(strsplit(names(",
            "rmdCompAnn), \"_vs_\"), \"[\", 1))),\n",
            "  p.value=rep(unlist(lapply(rmdCompAnn,",
            " function(comp) comp$p.value)), 2),\n",
            "  cor=rep(unlist(lapply(rmdCompAnn, ",
            "function(comp) comp$correlation)), ",
            "2))\n",
            "ggplotdata$correlation <- factor",
            "(rep(\"ns\", nrow(ggplotdata)), ",
            "levels=c(\"ns\", \"*\", \"**\", ",
            "\"***\"), ordered=TRUE)\n",
            "w <- ggplotdata$p.value <= 0.05 & ",
            "(is.na(ggplotdata$cor) | ",
            "ggplotdata$cor >= 0.5)\n",
            "ggplotdata$correlation[w] <- \"*\"\n",
            "ggplotdata$correlation[w &",
            " ggplotdata$p.value <= 0.01] <- ",
            "\"**\"\n",
            "ggplotdata$correlation[w & ",
            "ggplotdata$p.value <= 0.001] <- ",
            "\"***\"\n",
            "ggplot(ggplotdata, aes(x=val1, ",
            "y=val2, fill=correlation)) + ",
            "geom_tile() + labs(x=\"\", y=\"\") +\n",
            "  theme(axis.text.x=element_text(",
            "angle = 45, vjust = 1, hjust=1), ",
            "aspect.ratio=1) +\n",
            "  scale_fill_manual(values=c(",
            "\"lightgray\", \"gold\", \"orange\", ",
            "\"darkorange\")[summary(",
            "ggplotdata$correlation) > 0])"
        ),
        "Heatmap comparisons between attached values"
    ),
    append = TRUE, file = rmdFile
    )

    # table with significant comparisons only
    tabledata <- data.frame(
        Attached.Values.1 = unlist(lapply(strsplit(
            names(rmdCompAnn), "_vs_"
        ), "[", 1)),
        Attached.Values.2 = unlist(lapply(strsplit(
            names(rmdCompAnn), "_vs_"
        ), "[", 2)),
        p.value = unlist(lapply(
            rmdCompAnn,
            function(comp) comp$p.value
        )),
        cor = unlist(lapply(
            rmdCompAnn,
            function(comp) comp$correlation
        ))
    )

    tabledata$correlation <- factor(rep("ns", nrow(tabledata)),
        levels = c("ns", "*", "**", "***"),
        ordered = TRUE
    )
    w <- tabledata$p.value <= 0.05 & (is.na(tabledata$cor) |
        tabledata$cor >= 0.5)
    tabledata$correlation[w] <- "*"
    tabledata$correlation[w & tabledata$p.value <= 0.01] <- "**"
    tabledata$correlation[w & tabledata$p.value <= 0.001] <- "***"
    w <- which(tabledata$correlation != "ns")
    tabledata <- tabledata[w, ]
    if (nrow(tabledata) > 0) {
        tabledata$correlation <-
            Cgt.intlinkstring(tabledata$correlation, paste0("comp", w))
        write(Cgt.tablestring(tabledata[, c(
            "Attached.Values.1", "Attached.Values.2",
            "correlation"
        )], row.names = FALSE), append = TRUE, file = rmdFile)

        write(paste0(
            "\n* ns = not significant (p-value > 0.05 or cor.",
            "coef. < 0.5)\n* <span>*</span> = adjusted p-value ",
            "$\\leq$ 0.05\n* <span>**</span> = adjusted p-value",
            " $\\leq$ 0.01\n* <span>***</span> = adjusted",
            " p-value $\\leq$ 0.001\n"
        ), append = TRUE, file = rmdFile)
    } else {
        write("No significant correlations were found. \n",
            append = TRUE, file = rmdFile
        )
    }

    # only significant values (and !is.na)
    pvals <- unlist(lapply(rmdCompAnn, function(tmp) tmp$p.value))
    cors <- unlist(lapply(rmdCompAnn, function(tmp) tmp$correlation))
    w <- which(pvals <= 0.05 & (is.na(cors) | cors >= 0.5))

    if (length(w) > 0) {
        for (pair in w) {
            write(paste0(
                "\n## Compare ", names(dat)[pairs[pair, 1]],
                " with ", names(dat)[pairs[pair, 2]],
                "{#comp", pair, "}"
            ),
            append = TRUE, file = rmdFile
            )
            if (!is.na(rmdCompAnn[[pair]]$correlation)) {
                write(paste0(
                    names(
                        rmdCompAnn[[pair]]$correlation
                    ), " is: ",
                    ifelse(is.numeric(rmdCompAnn[[pair]]$correlation),
                        Cgt.roundNicely(rmdCompAnn[[pair]]$correlation),
                        rmdCompAnn[[pair]]$correlation
                    ), "."
                ),
                append = TRUE, file = rmdFile
                )
            }
            write(Cgt.rcodestring(
                paste0(
                    "plotdata <- rmdCompAnn[[", pair,
                    "]]$ggplotdata\nggplot(plotdata",
                    rmdCompAnn[[pair]]$ggplotstring
                ),
                paste0("Step4-", pair, "-", nrow(pairs))
            ),
            append = TRUE, file = rmdFile
            )
            write(paste0("\n", Cgt.intlinkstring(
                "back to overview",
                "summarycomp"
            )),
            append = TRUE, file = rmdFile
            )
        }
    }

    ## save rmd-file and generate html or PDF ####
    outfile <- paste0(
        getwd(), "/", Sys.Date(), "_", name, "_report.",
        ifelse(outputFormat == "PDF", "pdf", "html")
    )

    rmarkdown::render(rmdFile, output_file = outfile, quiet = !verbose)
    if (verbose) {
        message("created: ", rmdFile, "\ncreated: ", outfile)
    }
}


################################################################################
############# internally used functions ########################################

## prints a progressbar on console, use for loops in verbose == TRUE
Cgt.progressbar <- function(i, all, name = NULL) {
    if (!is.numeric(i) || !is.numeric(all)) {
        stop("Cogito::Cgt.progressbar")
    }
    if (all[1] == i) {
        if (!is.null(name) && is.character(name))
            cat(name, "\n")
        cat("0%", paste0(rep("-", 21), collapse = ""), "25%",
            paste0(rep("-", 21), collapse = ""), "50%",
            paste0(rep("-", 22), collapse = ""), "75%",
            paste0(rep("-", 21), collapse = ""), "100%\n",
            sep = ""
        )
    }
    w <- seq(from = 0, to = 100, length.out = length(all) + 1)
    ww <- which(all == i)
    cat(paste0(rep("x", round(w[ww + 1]) - round(w[ww])), collapse = ""))
    if (all[length(all)] == i) {
        cat("\n")
    }
}

## aggregate attached values (ann) on genes
## case 1: scale "numeric" (mode "min", "max", "mean", "median")
## case 2: scale "ordered" (mode "median")
## case 3: scale "factor" or "logical" (mode "modus")
## case 4: scale "logical" (mode "any")
Cgt.aggr <- function(genes, ann, scale, mode, maxdist){
    # check if arguments genes and ann are GRanges and ann has mcols
    if (!methods::is(genes, "GRanges") || !methods::is(ann, "GRanges") ||
        ncol(mcols(ann)) == 0) {
        stop(
            "Cogito::Cgt.aggr ",
            "parameter genes and parameter ann ",
            "have to be GRanges, ann have to have mcols"
        )
    }
    # check argument scale
    if (!is.character(scale) || length(scale) != 1 || !scale %in%
        c("numeric", "ordered", "factor", "logical")) {
        stop(
            "Cogito::Cgt.aggr ",
            "parameter scale has to be character with length 1, ",
            "and \"numeric\", \"ordered\", \"factor\" or \"logical\""
        )
    }
    # check argument mode
    if (!is.character(mode) || length(scale) != 1 ||
        !mode %in% c("min", "max", "mean", "median", "modus", "any")) {
        stop(
            "Cogito::Cgt.aggr ",
            "parameter mode has to be character with length 1 ",
            "and \"min\", \"max\", \"mean\", \"median\", \"any\" or \"modus\""
        )
    }
    # check argument maxdist
    if (!is.numeric(maxdist) || maxdist < 0) {
        stop(
            "Cogito::Cgt.aggr ",
            "parameter maxdist has to be numeric and >= 0"
        )
    }
    # check combination of scale and mode and combination with mcols of ann
    ## case 1: scale "numeric" (mode "min", "max", "mean", "median")
    if (!(scale == "numeric" && mode %in% c("min", "max", "mean", "median") &&
        all(unlist(lapply(mcols(ann), is.numeric)))) &&
        ## case 2: scale "ordered" (mode "median")
        !(scale == "ordered" && mode %in% c("median") &&
            all(unlist(lapply(mcols(ann), is.ordered)))) &&
        ## case 3: scale "factor" or "logical" (mode "modus")
        !(scale %in% c("factor", "logical") && mode == "modus" &&
            (all(unlist(lapply(mcols(ann), is.factor))) ||
                all(unlist(lapply(mcols(ann), is.logical))))) &&
        ## case 4: scale "logical" (modus "any")
        !(scale == "logical" && mode == "any")) {
        stop(
            "Cogito::Cgt.aggr ",
            "parameter combination scale and mode have to be: ",
            "case 1 numeric with min, max, mean or median, ",
            "case 2 ordered with median, ",
            "case 3 factor or logical with modus ",
            "or case 4 logical with modus any"
        )
    }

    dtn <- distanceToNearest(ann, resize(genes, 1))
    dtn <- dtn[as.data.frame(dtn)[, 3] <= maxdist]
    dtns <- split(S4Vectors::queryHits(dtn), S4Vectors::subjectHits(dtn))
    sgl <- unlist(lapply(dtns, length)) == 1
    for (i in seq_along(mcols(ann))) {
        genes$blub <- NA
        if (scale %in% c("factor", "ordered")) {
            genes$blub <-
                factor(genes$blub,
                    levels = levels(mcols(ann)[, i]),
                    ordered = is.ordered(mcols(ann)[, 1])
                )
        }
        # calculate values, which have only one overlap
        genes$blub[as.integer(names(dtns)[sgl])] <-
            mcols(ann)[unlist(dtns[sgl]), i]
        # calculate values, which have more than one overlap
        mydf <- mcols(ann)[, i]
        # scale="numeric" and mode %in% "max", "min", "median", "mean" case 1
        if (scale == "numeric") {
            genes$blub[as.integer(names(dtns)[!sgl])] <-
                unlist(lapply(
                    dtns[!sgl],
                    function(idx) match.fun(mode)(mydf[idx])
                ))
        } else if (scale == "ordered") { # and mode == "median" case 2
            genes$blub[as.integer(names(dtns)[!sgl])] <-
                lapply(dtns[!sgl], function(idx) Cgt.medianOrdered(mydf[idx]))
        } else if (scale %in% c("factor", "logical") && mode == "modus") {
            # case 3
            genes$blub[as.integer(names(dtns)[!sgl])] <-
                unlist(lapply(dtns[!sgl], function(idx) Cgt.modus(mydf[idx])))
        } else { # case 4
            genes$blub[as.integer(names(dtns)[!sgl])] <-
                unlist(lapply(dtns[!sgl], function(idx) Cgt.any(mydf[idx])))
        }
        colnames(mcols(genes))[
            ncol(mcols(genes))
        ] <- colnames(mcols(ann))[i]
    }
    return(genes)
}

Cgt.medianOrdered <- function(x, na.rm = TRUE) {
    if (length(x) == 1) {
        return(x)
    }
    if (!is.factor(x) || !is.ordered(x)) {
        stop("Cogito::Cgt.medianOrdered")
    }
    if (is.na(na.rm) || !is.logical(na.rm)) {
        na.rm <- TRUE
    }

    if (all(is.na(x)) || (!na.rm && sum(is.na(x)) >= length(x) / 2)) {
        return(NA)
    }
    return(levels(x)[min(which(cumsum(summary(x)) >=
        (length(x) - sum(is.na(x))) / 2))])
}

Cgt.any <- function(x) {
    if (!is.logical(x) && !is.factor(x)) {
        stop("Cogito::Cgt.any")
    }
    if (all(is.na(x))) {
        return(NA)
    }
    return(any(x, na.rm = TRUE))
}

Cgt.modus <- function(x, na.rm = TRUE) {
    if (!is.logical(x) && !is.factor(x)) {
        stop("Cogito::Cgt.modus")
    }
    if (is.na(na.rm) || !is.logical(na.rm)) {
        na.rm <- TRUE
    }
    if (all(is.na(x))) {
        return(NA)
    }
    if (is.logical(x) && !na.rm) {
        return(c(TRUE, FALSE, NA)[which.max(c(
            sum(x, na.rm = TRUE),
            sum(!x, na.rm = TRUE),
            sum(is.na(x))
        ))])
    }
    if (is.logical(x) && na.rm) {
        return(c(TRUE, FALSE)[which.max(c(
            sum(x, na.rm = TRUE),
            sum(!x, na.rm = TRUE)
        ))])
    }
    if (!na.rm || !any(is.na(x))) {
        return(levels(x)[which.max(summary(x))])
    }
    return(levels(x)[which.max(summary(x)[-length(summary(x))])])
}

Cgt.plotAttributeNominal <- function(attributes, usedLevelsOnly = TRUE,
                                        relativScale = TRUE, sort = FALSE,
                                        maxNrLevels = 10, na.rm = FALSE,
                                        groups = NULL, main = "") {
    if (!is.list(attributes)) {
        attributes <- list(attributes = attributes)
    }
    if (is.null(names(attributes))) {
        names(attributes) <- paste0("attributes ", seq_along(attributes))
    } else {
        names(attributes)[names(attributes) == ""] <-
            paste0("attributes ", seq_along(attributes))[
                names(attributes) == ""
            ]
    }
    attributes <- lapply(attributes, function(x) {
        if (!is.factor(x)) {
            return(factor(x))
        } else {
            return(x)
        }
    })
    if (!is.null(groups) &&
        (!is.list(groups) || !all(unlist(lapply(groups, is.character))))) {
        groups <- NULL
    }
    if (!is.null(groups)) {
        groups <- lapply(
            groups,
            function(group) group[group %in% names(attributes)]
        )
        groups$Other <- names(attributes)[
            !names(attributes) %in% unlist(groups)
        ]
        groups <- groups[unlist(lapply(groups, length)) > 0]
        attributes <- attributes[unlist(groups)]
    }

    # all same levels
    for (i in seq_along(attributes)) {
        levels(attributes[[i]]) <- unique(unlist(lapply(attributes, levels)))
    }

    if (na.rm) {
        attributes <- lapply(attributes, function(x) x[!is.na(x)])
        attNA <- rep(FALSE, length(attributes))
        names(attNA) <- names(attributes)
    } else {
        attNA <- unlist(lapply(attributes, function(x) sum(is.na(x))))
    }

    if (maxNrLevels < 2) {
        warning(
            "Cogito::Cgt.plotAttributeNominal ",
            "parameter maxNrLevels < 2, set to default 10"
        )
        maxNrLevels <- 10
    }

    if (length(attributes) == 1) {
        plotdata <- summary(attributes[[1]], maxsum = maxNrLevels)
        if (sort) {
            plotdata <- sort(plotdata, decreasing = TRUE)
        }
        if (usedLevelsOnly) {
            plotdata <- plotdata[plotdata > 0]
        }
        n <- names(plotdata)
        plotdata <- matrix(plotdata, 1, length(plotdata))
        colnames(plotdata) <- n
    } else {
        plotdata <- t(matrix(unlist(lapply(attributes, function(att) {
            if (any(attNA) & all(!is.na(att))) {
                return(c(summary(att, maxsum = maxNrLevels), 0))
            } else {
                return(summary(att, maxsum = maxNrLevels))
            }
        })), ncol = length(attributes)))
        if (any(attNA)) {
            colnames(plotdata) <- names(summary(attributes[[
            which(attNA > 0)[1]]], maxNrLevels))
        } else {
            colnames(plotdata) <- names(summary(attributes[[1]], maxNrLevels))
        }
        if (sort) {
            plotdata <- plotdata[, order(plotdata[
                which.max(apply(plotdata, 1, sum)),
            ], decreasing = TRUE)]
        }
        if (usedLevelsOnly) {
            plotdata <- plotdata[, colSums(plotdata) > 0]
        }
    }

    # shorten names, if all prefixes are equal to main title
    if (!is.null(main) && is.character(main) && main != "" &&
        all(substr(names(attributes), 1, width(main)) == main)) {
        names(attributes) <-
            sub("^[[:punct:]]*", "", substr(
                names(attributes), width(main) + 1,
                max(width(names(attributes)))
            ))
        groups <- lapply(
            groups, function(group) {
                sub("^[[:punct:]]*", "", substr(
                    group, width(main) + 1,
                    max(width(group))
                ))
            }
        )
    }
    # set y-axis label to postfix, if it is equal in all samples
    # and shorten names
    xlab <- "value"
    if (length(attributes) > 1) {
        postfix <-
            unique(unlist(lapply(
                strsplit(names(attributes), "[[:punct:]]"), function(sstr) {
                    sstr[length(sstr)]
                }
            )))
        if (length(postfix) == 1) {
            xlab <- postfix
            names(attributes) <-
                sub(
                    "[[:punct:]]*$", "",
                    substr(
                        names(attributes), 1,
                        width(names(attributes)) - width(postfix)
                    )
                )
            groups <- lapply(groups, function(group) {
                sub("[[:punct:]]*$", "", substr(
                    group, 1,
                    width(group) - width(postfix)
                ))
            })
        }
    }
    rownames(plotdata) <- names(attributes)

    # location parameter modus
    loc <- colnames(plotdata)[apply(plotdata, 1, which.max)]
    # dispersion parameter nr of used levels (also NA's)
    dis <- apply(as.matrix(apply(plotdata, 1, "!=", 0)), 2, sum)
    if (length(attributes) == 1) {
        names(loc) <- "modus"
        names(dis) <- "number of used levels"
    } else {
        names(loc) <- rownames(plotdata)
    }

    # ggplot2
    ggplotdata <- data.frame(
        quantity = as.vector(plotdata),
        value = factor(rep(colnames(plotdata),
            each = nrow(plotdata)
        ),
        levels = colnames(plotdata)
        ),
        sample = rep(rownames(plotdata), ncol(plotdata))
    )
    if (max(m <- (apply(plotdata / apply(plotdata, 1, sum), 2, min))) > 0.9) {
        if (relativScale) {
            m <- max(m)
        } else {
            m <- min(plotdata[, which.max(m)])
        }
        ggplotstring <- paste0(
            ", aes(x=value, y=quantity, fill=sample)) + \n",
            "geom_bar(data=plotdata %>% mutate(ss=\"all\"),",
            " stat=\"identity\", position=position_dodge())",
            " + \ngeom_bar(data=plotdata %>% filter(",
            "quantity <", m, ") %>% mutate(ss=\"subset\"), ",
            "stat=\"identity\", ",
            "position=position_dodge()) +",
            " \nfacet_wrap(~ ss, scales=\"free_y\") + xlab(\"", xlab, "\")"
        )
    } else {
        ggplotstring <- paste0(
            ", aes(x=value, y=quantity, fill=sample)) + ",
            "geom_bar(stat=\"identity\", ",
            "position=position_dodge()) + xlab(\"", xlab, "\")"
        )
    }

    if (!is.null(main) && is.character(main) && main != "") {
        ggplotstring <- paste0(ggplotstring, " + labs(title=\"", main, "\")")
    }

    if (!is.null(groups)) {
        ggplotdata$group <- NA
        for (group in seq_along(groups)) {
            ggplotdata$group[ggplotdata$sample %in% groups[[group]]] <-
                names(groups)[group]
        }
        ggplotstring <- paste(ggplotstring, "+ facet_wrap(~group)")
    }

    if (relativScale) {
        ggplotdata$quantity <- as.vector(plotdata / apply(plotdata, 1, sum))
    }

    result <- list(
        data = plotdata, n = unlist(lapply(attributes, length)),
        location = loc, dispersion = dis,
        ggplotdata = ggplotdata, ggplotstring = ggplotstring
    )
    return(result)
}

Cgt.plotAttributeOrdinal <- function(attributes, usedLevelsOnly = FALSE,
                                    relativScale = TRUE, maxLevelsWindow = 5,
                                    na.rm = FALSE, groups = NULL, main = "") {
    if (!is.list(attributes)) {
        attributes <- list(attributes = attributes)
    }
    if (is.null(names(attributes))) {
        names(attributes) <- paste0("attribute ", seq_along(attributes))
    }
    attributes <- lapply(attributes, function(x) {
        if (!is.factor(x)) {
            return(factor(x))
        } else {
            return(x)
        }
    })
    if (na.rm) {
        attributes <- lapply(attributes, function(x) x[!is.na(x)])
        attNA <- rep(FALSE, length(attributes))
        names(attNA) <- names(attributes)
    } else {
        attNA <- unlist(lapply(attributes, function(x) sum(is.na(x)) > 0))
    }

    # all levels the same and ordered
    if (length(attributes) == 1) {
        attributes[[1]] <- factor(attributes[[1]], ordered = TRUE)
    } else {
        # if some levels(attribute) contains all other and
        # is ordered take it, else make new order
        w <- which(unlist(lapply(attributes, is.ordered)) & !
            vapply(seq_along(attributes), function(atti) {
                any(unlist(
                    lapply(lapply(attributes, levels), function(x) {
                        any(!x %in% levels(attributes[[atti]]))
                    })
                ))
            }, c(TRUE)))
        if (length(w) > 0) {
            attributes <- lapply(attributes, factor,
                levels = levels(attributes[[w[1]]]),
                ordered = TRUE
            )
        } else {
            attributes <-
                lapply(attributes, factor,
                        levels = unique(unlist(lapply(attributes, levels))),
                        ordered = TRUE)
        }
    }

    if (maxLevelsWindow < 1) {
        warning(
            "Cogito::Cgt.plotAttributeOrdinal ",
            "parameter maxLevelsWindow < 1, set to default 5"
        )
        maxLevelsWindow <- 5
    }

    if (length(attributes) == 1) {
        plotdata <- summary(attributes[[1]], maxsum = nlevels(attributes[[1]]) +
            ifelse(attNA, 1, 0))[seq_along(levels(attributes[[1]]))]
        if (usedLevelsOnly) {
            plotdata <- plotdata[plotdata > 0]
        }
        median <- S4Vectors::runValue(
            S4Vectors::Rle(names(plotdata), plotdata)[sum(plotdata) / 2 + 0.5]
        )
        medianIx <- which(names(plotdata) == median)
        plotdataResult <- plotdata
        # window with width <= maxLevelsWindow centered at median
        if (medianIx - maxLevelsWindow > 1) {
            plotdata[medianIx - maxLevelsWindow] <-
                sum(plotdata[seq_len(medianIx - maxLevelsWindow)])
            plotdata <- plotdata[(medianIx - maxLevelsWindow):length(plotdata)]
            names(plotdata)[1] <- "other (<)"
            medianIx <- medianIx - (medianIx - maxLevelsWindow - 1)
        }
        if (medianIx + maxLevelsWindow < length(plotdata)) {
            plotdata[medianIx + maxLevelsWindow] <-
                sum(plotdata[(medianIx + maxLevelsWindow):length(plotdata)])
            plotdata <- plotdata[seq_len(medianIx + maxLevelsWindow)]
            names(plotdata)[length(plotdata)] <- "other (>)"
        }
        if (attNA) {
            plotdata <- c(plotdata, sum(is.na(attributes[[1]])))
            names(plotdata)[length(plotdata)] <- "NA's"
            plotdataResult <- c(plotdataResult, sum(is.na(attributes[[1]])))
            names(plotdataResult)[length(plotdataResult)] <- "NA's"
        }
        n <- names(plotdata)
        plotdata <- matrix(plotdata, 1, length(plotdata))
        colnames(plotdata) <- n
        n <- names(plotdataResult)
        plotdataResult <- matrix(plotdataResult, 1, length(plotdataResult))
        colnames(plotdataResult) <- n
        rownames(plotdataResult) <- names(attributes)
    } else {
        if (any(attNA)) {
            plotdata <- t(vapply(seq_along(attributes), function(atti) {
                if (attNA[atti]) {
                    return(summary(attributes[[atti]]))
                } else {
                    return(c(summary(attributes[[atti]]), 0))
                }
            }, rep(0, nlevels(attributes[[1]]) + 1)))
            colnames(plotdata) <- names(summary(attributes[[which(attNA)[1]]]))
        } else {
            plotdata <- t(vapply(
                seq_along(attributes),
                function(atti) summary(attributes[[atti]]),
                rep(0, nlevels(attributes[[1]]))
            ))
        }
        if (usedLevelsOnly) {
            plotdata <- plotdata[, colSums(plotdata) > 0]
        }
    }

    # ggplot
    # shorten names, if all prefixes are equal to main title
    if (!is.null(main) && is.character(main) && main != "" &&
        all(substr(names(attributes), 1, width(main)) == main)) {
        names(attributes) <-
            sub("^[[:punct:]]*", "", substr(
                names(attributes), width(main) + 1,
                max(width(names(attributes)))
            ))
        groups <- lapply(
            groups, function(group) {
                sub("^[[:punct:]]*", "", substr(
                    group, width(main) + 1,
                    max(width(group))
                ))
            }
        )
    }
    # set y-axis label to postfix, if it is equal in all samples
    # and shorten names
    xlab <- "value"
    if (length(attributes) > 1) {
        postfix <-
            unique(unlist(lapply(
                strsplit(names(attributes), "[[:punct:]]"), function(sstr) {
                    sstr[length(sstr)]
                }
            )))
        if (length(postfix) == 1) {
            xlab <- postfix
            names(attributes) <-
                sub(
                    "[[:punct:]]*$", "",
                    substr(
                        names(attributes), 1,
                        width(names(attributes)) - width(postfix)
                    )
                )
            groups <- lapply(groups, function(group) {
                sub("[[:punct:]]*$", "", substr(
                    group, 1,
                    width(group) - width(postfix)
                ))
            })
        }
    }
    rownames(plotdata) <- names(attributes)

    ggplotdata <- data.frame(
        quantity = as.vector(plotdata),
        value = factor(
            rep(colnames(plotdata),
                each = nrow(plotdata)
            ),
            colnames(plotdata)
        ),
        sample = rep(rownames(plotdata), ncol(plotdata))
    )
    ggplotstring <- paste0(
        ", aes(x=value, y=quantity, fill=sample)) + ",
        "geom_bar(stat=\"identity\", ",
        "position=position_dodge()) + xlab(\"", xlab, "\")"
    )

    if (!is.null(main) && is.character(main) && main != "") {
        ggplotstring <- paste0(ggplotstring, " + labs(title=\"", main, "\")")
    }

    if (!is.null(groups)) {
        ggplotdata$group <- NA
        for (group in seq_along(groups)) {
            ggplotdata$group[ggplotdata$sample %in% groups[[group]]] <-
                names(groups)[group]
        }
        ggplotstring <- paste(ggplotstring, "+ facet_wrap(~group)")
    }

    if (relativScale) {
        ggplotdata$quantity <- as.vector(plotdata / apply(plotdata, 1, sum))
    }

    if (length(attributes) == 1) {
        loc <- median
        names(loc) <- "median"
        plotdata <- plotdataResult
    } else {
        if (any(attNA)) {
            loc <- apply(plotdata[, -ncol(plotdata)], 1, function(x) {
                S4Vectors::runValue(S4Vectors::Rle(
                    colnames(plotdata[, -ncol(plotdata)]), x
                )[sum(x) / 2])
            })
        } else {
            loc <- apply(plotdata, 1, function(x) {
                S4Vectors::runValue(S4Vectors::Rle(
                    colnames(plotdata), x
                )[sum(x) / 2])
            })
        }
    }
    if (any(attNA)) {
        dis <- apply(matrix(plotdata[, -ncol(plotdata)],
            nrow = nrow(plotdata)
        ), 1, function(x) {
            if (sum(is.na(x)) > 0) {
                max(which(x[-length(x)] != 0)) -
                    min(which(x[-length(x)] != 0)) + 1
            } else {
                max(which(x != 0)) - min(which(x != 0)) + 1
            }
        })
    } else {
        dis <- apply(plotdata, 1, function(x) {
            if (sum(is.na(x)) > 0) {
                max(which(x[-length(x)] != 0)) -
                    min(which(x[-length(x)] != 0)) + 1
            } else {
                max(which(x != 0)) - min(which(x != 0)) + 1
            }
        })
    }
    if (length(attributes) == 1) {
        names(dis) <- "number of values between lowest and highest used value"
    }

    result <- list(
        data = plotdata, n = unlist(lapply(attributes, length)),
        location = loc, dispersion = dis, ggplotdata = ggplotdata,
        ggplotstring = ggplotstring
    )
    return(result)
}

Cgt.plotAttributeCont <- function(attributes, outline = FALSE, scale = "int",
                                    groups = NULL, main = "") {
    if (!is.list(attributes)) {
        attributes <- list(attributes = attributes)
    }
    if (is.null(names(attributes))) {
        names(attributes) <- paste0("attributes ", seq_along(attributes))
    }
    attributes <- lapply(attributes, function(x) {
        if (!is.numeric(x)) {
            return(as.numeric(x))
        } else {
            return(x)
        }
    })

    if (!is.null(groups) &&
        (!is.list(groups) || !all(unlist(lapply(groups, is.character))))) {
        groups <- NULL
    }
    if (!is.null(groups)) {
        groups <- lapply(groups, function(group) {
            group[group %in% names(attributes)]
        })
        groups$Other <-
            names(attributes)[!names(attributes) %in% unlist(groups)]
        groups <- groups[unlist(lapply(groups, length)) > 0]
        attributes <- attributes[unlist(groups)]
    }

    if (scale == "rat") {
        # location parameter rat geom mean
        loc <- unlist(lapply(attributes, function(x) {
            exp(sum(log(x[x > 0]), na.rm = TRUE) / sum(!is.na(x)))
        }))
        # dispersion parameter rat attribute (norm.)
        # (coefficient of variation standardVariation/mean)
        dis <- unlist(lapply(
            attributes,
            function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
        ))
        if (length(attributes) == 1) {
            names(loc) <- "geometric mean"
            names(dis) <- "coefficient of variation"
        }
    } else {
        # location parameter int arith mean
        loc <- unlist(lapply(
            attributes,
            function(x) mean(x, na.rm = TRUE) / sum(!is.na(x))
        ))
        # dispersion parameter int attribute standard deviation
        dis <- unlist(lapply(attributes, function(x) sd(x, na.rm = TRUE)))
        if (length(attributes) == 1) {
            names(loc) <- "arithmetic mean"
            names(dis) <- "standard deviation"
        }
    }

    # how long whiskers depends on range
    range <- 1.5
    plotdata <- boxplot(attributes, plot = FALSE, range = range)
    tmp <- which(vapply(seq_len(max(width(plotdata$names))), function(prefix) {
        all(unlist(gregexpr(
            substr(plotdata$names[1], 1, prefix),
            plotdata$names
        )) == 1)
    }, c(TRUE)))
    if (length(tmp > 0)) {
        plotdata$names <- substr(
            plotdata$names, max(tmp) + 1,
            max(width(plotdata$names))
        )
    }
    tempDat <- t(plotdata$stats)
    rownames(tempDat) <- names(attributes)
    colnames(tempDat) <- c(
        paste0(">-range", range), "25%", "50%", "75%", paste0("<+range", range)
    )
    result <- list(
        data = tempDat, n = unlist(lapply(attributes, length)),
        location = loc, dispersion = dis
    )

    # ggplot
    # shorten names, if all prefixes are equal to main title
    if (!is.null(main) && is.character(main) && main != "" &&
        all(substr(names(attributes), 1, width(main)) == main)) {
        names(attributes) <-
            sub("^[[:punct:]]*", "", substr(
                names(attributes), width(main) + 1,
                max(width(names(attributes)))
            ))
        groups <- lapply(
            groups, function(group) {
                sub("^[[:punct:]]*", "", substr(
                    group, width(main) + 1,
                    max(width(group))
                ))
            }
        )
    }
    # set y-axis label to postfix, if it is equal in all samples
    # and shorten names
    ylab <- "value"
    if (length(attributes) > 1) {
        postfix <-
            unique(unlist(lapply(
                strsplit(names(attributes), "[[:punct:]]"), function(sstr) {
                    sstr[length(sstr)]
                }
            )))
        if (length(postfix) == 1) {
            ylab <- postfix
            names(attributes) <-
                sub(
                    "[[:punct:]]*$", "",
                    substr(
                        names(attributes), 1,
                        width(names(attributes)) - width(postfix)
                    )
                )
            groups <- lapply(groups, function(group) {
                sub("[[:punct:]]*$", "", substr(
                    group, 1,
                    width(group) - width(postfix)
                ))
            })
        }
    }
    ggplotdata <- data.frame(
        value = unlist(attributes),
        sample = factor(rep(
            names(attributes),
            lapply(attributes, length)
        ))
    )
    ggplotdata <- ggplotdata[!is.na(ggplotdata$value), ]
    if (!is.null(groups)) {
        ggplotdata$group <- NA
        for (group in seq_along(groups)) {
            ggplotdata$group[ggplotdata$sample %in% groups[[group]]] <-
                names(groups)[group]
        }
    }
    ggplotstring <- paste0(
        ", aes(x=sample, y=value",
        ifelse(!is.null(groups), ", fill=group", ""),
        ")) + geom_boxplot(", ifelse(outline, ")",
            paste0(
                "outlier.shape=NA) + coord_cartesian(",
                "ylim=c(", min(tempDat[, 1]), ", ",
                max(tempDat[, 5]), "))"
            )
        ), " + xlab(\"sample\") + ylab(\"", ylab, "\")"
    )
    if (length(attributes) > 1){
        ggplotstring <- 
            paste(ggplotstring, "+ theme(axis.text.x = ",
                    "element_text(angle = 45, vjust = 1, hjust = 1))")
    }

    if (!is.null(main) && is.character(main) && main != "") {
        ggplotstring <- paste0(ggplotstring, " + labs(title=\"", main, "\")")
    }

    result$ggplotdata <- ggplotdata
    result$ggplotstring <- ggplotstring
    return(result)
}

Cgt.plotComparisonContToCont <- function(datComp, scale1 = "int",
                                            scale2 = "int") {
    if (!is.list(datComp) || !length(datComp) >= 2 ||
        length(datComp[[1]]) != length(datComp[[2]]) ||
        sum(!is.na(datComp[[1]])) == 0 || sum(!is.na(datComp[[2]])) == 0) {
        stop(
            "Cogito::Cgt.plotComparisonContToCont ",
            "parameter datComp was no list or has length < 2 ",
            "or different lengths or one value was only NA"
        )
    }
    if (!is.character(scale1) || !scale1 %in% c("int", "rat")) {
        warning(
            "Cogito::Cgt.plotComparisonContToCont ",
            "parameter scale1 was no character ",
            "or not in \"int\", \"rat\" and set to \"int\"."
        )
        scale1 <- "int"
    }
    if (!is.character(scale2) || !scale2 %in% c("int", "rat")) {
        warning(
            "Cogito::Cgt.plotComparisonContToCont ",
            "parameter scale2 was no character ",
            "or not in \"int\", \"rat\" and set to \"int\"."
        )
        scale2 <- "int"
    }
    
    if (sum(!is.na(datComp[[1]]) & !is.na(datComp[[2]])) <= 2) {
        return(list(
            data = NULL, correlation = NA, significant = FALSE,
            mutinf = NA, normmutinf = NA, p.value = NA, ggplotstring = ")"
        ))
    }
    
    cort <- stats::cor.test(datComp[[1]], datComp[[2]],
        method = "spearman",
        use = "na.or.complete", 
        exact = FALSE
    )
    cor <- cort$estimate
    names(cor) <- "Spearman's correlation coefficient"
    pval <- cort$p.value
    names(pval) <- "correlation test p-value"
    
    ## mutual information
    nrBins <- 11
    discDat1 <- datComp[[1]]
    discDat1[order(datComp[[1]])] <-
        (seq_along(datComp[[1]]) %/%
            round(length(datComp[[1]]) / nrBins + 0.5) + 1)
    discDat2 <- datComp[[2]]
    discDat2[order(datComp[[2]])] <-
        (seq_along(datComp[[2]]) %/%
            round(length(datComp[[2]]) / nrBins + 0.5) + 1)
    tab <- table(discDat1, discDat2)
    mutinf <- entropy::mi.plugin(tab)
    # normalized mutual information
    normmutinf <-
        mutinf / min(
            entropy::entropy.plugin(rowSums(tab)),
            entropy::entropy.plugin(colSums(tab))
        )

    mylm <- lm(datComp[[2]] ~ datComp[[1]])
    sign <- (!is.na(pval) && pval <= 0.05)

    plottitle <- paste("compare", names(datComp)[1], "with", names(datComp)[2])
    plotsubtitle <- paste(names(cor), round(cor, 3))

    ggplotdata <- as.data.frame(datComp)
    ggplotdata <-
        ggplotdata[!is.na(ggplotdata[, 1]) & !is.na(ggplotdata[, 2]), ]
    ggplotstring <- paste0(
        ", aes(x=", names(datComp)[1], ", y=",
        names(datComp)[2], ")) + geom_point() + ",
        "labs(title=\"", plottitle, "\", subtitle=\"",
        plotsubtitle, "\")"
    )

    return(list(
        data = mylm, correlation = cor, significant = sign, mutinf = mutinf,
        normmutinf = normmutinf, p.value = pval,
        ggplotdata = ggplotdata, ggplotstring = ggplotstring
    ))
}

Cgt.plotComparisonCatToCont <- function(datComp, scale1 = "nom",
                                        scale2 = "int", outline = FALSE) {
    if (!is.list(datComp) || !length(datComp) >= 2 ||
        length(datComp[[1]]) != length(datComp[[2]]) ||
        sum(!is.na(datComp[[1]])) == 0 || sum(!is.na(datComp[[2]])) == 0) {
        stop(
            "Cogito::Cgt.plotComparisonCatToCont ",
            "parameter datComp was no list ",
            "or has length < 2 or differnt lengths or one value was only NA"
        )
    }
    if (!is.character(scale1) ||
        !scale1 %in% c("ord", "nom", "bin", "int", "rat")) {
        warning(
            "Cogito::Cgt.plotComparisonCatToCont ",
            "parameter scale1 was no character or not in \"nom\", \"ord\", ",
            "\"bin\", \"int\", \"rat\" and set to \"nom\"."
        )
        scale1 <- "nom"
    }
    if (!is.character(scale2) ||
        !scale2 %in% c("ord", "nom", "bin", "int", "rat")) {
        warning(
            "Cogito::Cgt.plotComparisonCatToCont ",
            "parameter scale2 was no character or not in \"nom\", \"ord\", ",
            "\"bin\", \"int\", \"rat\" and set to \"int\"."
        )
        scale2 <- "int"
    }

    if (scale1 %in% c("int", "rat") && scale2 %in% c("nom", "ord", "bin")) {
        temp <- scale1
        scale1 <- scale2
        scale2 <- temp
        datComp <- datComp[2:1]
    }
    plotdata <- split(datComp[[2]], datComp[[1]], drop = TRUE)
    # remove groups with only NA's
    plotdata <- plotdata[lapply(plotdata, function(g) sum(!is.na(g))) != 0]

    ## mutual information
    if (length(plotdata) != 1) {
        nrBins <- 11
        discDat <- datComp[[2]]
        discDat[order(datComp[[2]])] <-
            (seq_along(datComp[[2]]) %/% round(length(datComp[[2]])
            / nrBins + 0.5) + 1)
        tab <- table(datComp[[1]], discDat)
        tab <- tab[apply(tab, 1, sum) > 0, apply(tab, 2, sum) > 0]
        mutinf <- entropy::mi.plugin(tab)
        # normalized mutual information
        normmutinf <- mutinf / min(
            entropy::entropy.plugin(rowSums(tab)),
            entropy::entropy.plugin(colSums(tab))
        )
    } else {
        mutinf <- NA
        normmutinf <- NA
    }

    ### "cat-->cont" ob cont in den verschiedenen cats verschieden ist
    if (length(plotdata) == 1) { # no test possible
        pval <- NA
        names(pval) <- "no test"
    } else if (length(plotdata) == 2) { # for example if bin
        pval <- wilcox.test(plotdata[[1]], plotdata[[2]],
            alternative = "two.sided",
            exact = FALSE
        )$p.value
        names(pval) <- "Wilcox rank sum test p-value"
        # if normally distributed (test with shapiro.test),
        # t-test also possible
    } else {
        pval <- kruskal.test(plotdata)$p.value
        names(pval) <- "Kruskal-Wallis rank sum test p-value"
    }

    # how long whiskers depends on range
    range <- 1.5
    plottitle <- paste("compare", names(datComp)[1], "with", names(datComp)[2])
    plotsubtitle <- ifelse(is.na(pval), "", paste(names(pval), round(pval, 3)))
    bxpdata <- boxplot(plotdata, range = range, plot = FALSE)

    bxpdata <- t(bxpdata$stats)
    rownames(bxpdata) <- names(plotdata)
    colnames(bxpdata) <- c(
        paste0(">-range", range), "25%", "50%",
        "75%", paste0("<+range", range)
    )

    ggplotdata <- as.data.frame(datComp)
    ggplotdata <- ggplotdata[!is.na(ggplotdata[, 2]), ] # not displayable
    ggplotstring <-
        paste0(
            ", aes(x=", names(datComp)[1], ", y=",
            names(datComp)[2], ")) + ", "labs(title=\"", plottitle,
            "\", subtitle=\"", plotsubtitle, "\") + geom_boxplot(",
            ifelse(outline, ")", paste0(
                "outlier.shape=NA) + ",
                "coord_cartesian(ylim=c(",
                min(bxpdata[, 1]), ", ",
                max(bxpdata[, 5]), "))"
            ))
        )

    return(list(
        data = bxpdata, correlation = NA,
        significant = (!is.na(pval) && pval <= 0.05), mutinf = mutinf,
        normmutinf = normmutinf, p.value = pval,
        ggplotdata = ggplotdata, ggplotstring = ggplotstring
    ))
}

Cgt.plotComparisonCatToCat <- function(datComp, scale1 = "nom",
                                        scale2 = "nom") {
    if (!is.list(datComp) || !length(datComp) >= 2 ||
        length(datComp[[1]]) != length(datComp[[2]]) ||
        sum(!is.na(datComp[[1]])) == 0 || sum(!is.na(datComp[[2]])) == 0) {
        stop(
            "Cogito::Cgt.plotComparisonCatToCat ",
            "parameter datComp was no list or ",
            "has length < 2 or different lengths or one value was only NA"
        )
    }
    if (!is.character(scale1) || !scale1 %in% c("ord", "nom", "bin")) {
        warning(
            "Cogito::Cgt.plotComparisonCatToCat ",
            "parameter scale1 was no character or not in \"nom\", ",
            "\"ord\", \"bin\" and set to \"nom\"."
        )
        scale1 <- "nom"
    }
    if (!is.character(scale2) || !scale2 %in% c("ord", "nom", "bin")) {
        warning(
            "Cogito::Cgt.plotComparisonCatToCat ",
            "parameter scale2 was no character or not in \"nom\", ",
            "\"ord\", \"bin\" and set to \"nom\"."
        )
        scale2 <- "nom"
    }

    if (scale1 == "ord" && scale2 != "ord") {
        scale1 <- scale2
        scale2 <- "ord"
        datComp <- datComp[2:1]
    }
    plotdata <- table(datComp[[1]], datComp[[2]])

    if (nrow(plotdata) <= 1 || ncol(plotdata) <= 1) {
        return(list(
            data = plotdata, correlation = NA, significant = FALSE,
            mutinf = NA, normmutinf = NA, p.value = NA, ggplotstring = ")"
        ))
    }
    if (any(w <- apply(plotdata, 1, sum) == 0)) {
        plotdata <- plotdata[-which(apply(plotdata, 1, sum) == 0), ]
        # none or only one row left
        if (any(dim(plotdata) == 0) || is.vector(plotdata)) {
            return(list(
                data = plotdata, correlation = NA, significant = FALSE,
                mutinf = NA, normmutinf = NA, p.value = NA, ggplotstring = ")"
            ))
        }
    }
    if (any(apply(plotdata, 2, sum) == 0)) {
        plotdata <- plotdata[, -which(apply(plotdata, 2, sum) == 0)]
        if (is.vector(plotdata)) { # only one col left
            return(list(
                data = plotdata, correlation = NA, significant = FALSE,
                mutinf = 0, normmutinf = 0, p.value = NA, ggplotstring = ")"
            ))
        }
    }
    cor <- NA
    if (scale1 %in% c("bin", "nom") && scale2 %in% c("bin", "nom")) {
        plotdata <- plotdata[
            order(apply(plotdata, 1, sum), decreasing = TRUE),
            order(apply(plotdata, 2, sum), decreasing = TRUE)
        ]
        if (sum(plotdata < 5) == 0) {
            pval <- chisq.test(plotdata)$p.value
            names(pval) <- "Pearsons's Chi-squared test p-value"
        } else {
            pval <- fisher.test(plotdata, alternative = "two.sided",
                                simulate.p.value = TRUE, B = 1e5)$p.value
            names(pval) <- "Fisher-exact test p-value"
        }
    } else if (scale1 %in% c("bin", "nom") && scale2 == "ord") {
        plotdata <-
            plotdata[order(apply(plotdata, 1, sum), decreasing = TRUE), ]
        pval <- kruskal.test(datComp[[2]], datComp[[1]])$p.value
        names(pval) <- "Kruskal-Wallis rank sum test p-value"
    } else {
        if (sum(!is.na(datComp[[1]]) & !is.na(datComp[[2]])) <= 2)
            cort <- list(estimate = cor(as.numeric(datComp[[1]]), 
                                        as.numeric(datComp[[2]]),
                                        method = "spearman", 
                                        use = "na.or.complete"),
                            p.value = NA)
        else
            cort <- stats::cor.test(as.numeric(datComp[[1]]), 
                                    as.numeric(datComp[[2]]),
                                    method = "spearman", 
                                    use = "na.or.complete", 
                                    exact = FALSE
                                    )
        cor <- cort$estimate
        names(cor) <- "Spearman's rank correlation coefficient"
        pval <- cort$p.value
        names(pval) <- "correlation test p-value"
    }
    ## mutual information
    mutinf <- entropy::mi.plugin(plotdata)
    # normalized mutual information
    normmutinf <- mutinf / min(
        entropy::entropy.plugin(rowSums(plotdata)),
        entropy::entropy.plugin(colSums(plotdata))
    )

    plottitle <- paste("compare", names(datComp)[1], "with", names(datComp)[2])
    plotsubtitle <- paste(names(pval), round(pval, 3))

    ggplotdata <- data.frame(
        val1 = rep(rownames(plotdata), ncol(plotdata)),
        val2 = rep(colnames(plotdata), each = nrow(plotdata)),
        number = as.vector(plotdata)
    )
    colnames(ggplotdata)[seq_len(2)] <- names(datComp)
    if (scale1 == "ord") {
        ggplotdata[, 1] <- factor(ggplotdata[, 1],
            levels = levels(datComp[[1]]),
            ordered = is.ordered(datComp[[1]])
        )
    }
    if (scale2 == "ord") {
        ggplotdata[, 2] <- factor(ggplotdata[, 2],
            levels = levels(datComp[[2]]),
            ordered = is.ordered(datComp[[2]])
        )
    }
    ggplotstring <- paste0(
        ", aes(x=", colnames(ggplotdata)[1], ", y=",
        colnames(ggplotdata)[2], ", fill=number)) + ",
        "geom_tile() + labs(title=\"", plottitle,
        "\", subtitle=\"", plotsubtitle, "\")"
    )

    return(list(
        data = plotdata, correlation = cor,
        significant = (!is.na(pval) && pval <= 0.05), mutinf = mutinf,
        normmutinf = normmutinf, p.value = pval, ggplotdata = ggplotdata,
        ggplotstring = ggplotstring
    ))
}

Cgt.getDataFromGRanges <- function(ranges, scales=NULL){
    if (!methods::is(ranges, "GRanges")) {
        stop("Cogito::Cgt.getDataFromGRanges")
    }
    dat <- mcols(ranges)
    # fit given scales to data
    if (ncol(dat) != 0) {
        scalesAlt <- Cgt.fitscale(dat)
        if (!is.character(scales)) {
            scales <- scalesAlt
        } else {
            if (ncol(dat) < length(scales)) {
                warning(
                    "Cogito::Cgt.getDataFromGRanges ",
                    "smaller number of columns of attached values than ",
                    "scale specifications: only first ones used"
                )
                scales <- scales[seq_along(dat)]
            } else if (ncol(dat) > length(scales)) {
                warning(
                    "Cogito::Cgt.getDataFromGRanges ",
                    "larger number of columns of attached values than ",
                    "scale specifications: other ones added"
                )
                scales <- c(scales, scalesAlt[-seq_along(scales)])
            }
            if (sum(w <-
                !scales %in% c("", "nom", "ord", "int", "rat", "bin"))) {
                warning(
                    "Cogito::Cgt.getDataFromGRanges ",
                    "bad scale identification in scales: that ones changed"
                )
                scales[w] <- scalesAlt[w]
            }
            if (sum(w <- factor(scales,
                levels = c("", "nom", "ord", "int", "rat", "bin"),
                ordered = TRUE
            ) >
                factor(scalesAlt,
                    levels = c("", "nom", "ord", "int", "rat", "bin"),
                    ordered = TRUE
                ))) {
                warning(
                    "Cogito::Cgt.getDataFromGRanges ",
                    "bad scale identification in scales: that ones changed"
                )
                scales[w] <- scalesAlt[w]
            }
        }
        dat <- dat[, scales != ""]
        scales <- scales[scales != ""]
    }

    # further analysis only makes sense if at least two different values appear
    if (length(w <- which(unlist(lapply(dat, function(x) {
        if (length(unique(x[seq_len(min(10, length(x)))])) > 1) {
            return(2)
        } else {
            return(length(unique(x)))
        }
    })) == 1)) > 0) {
        warning(
            "Cogito::Cgt.getDataFromGRanges ",
            "attribute(s) ", paste0(names(dat)[w], collapse = ", "),
            " has/have only one value: discarded for further analysis"
        )
        dat <- dat[, -w]
        scales <- scales[-w]
    }

    # if an attribute is not numeric and have all different values,
    # then it is a names column
    if (length(w <- which(unlist(lapply(dat, function(x) {
        if (is.numeric(x) ||
            length(unique(x[seq_len(min(10, length(x)))])) <
                min(10, length(x))) {
            return(2)
        } else {
            return(length(unique(x)))
        }
    })) == nrow(dat))) > 0) {
        warning(
            "Cogito::Cgt.getDataFromGRanges ",
            "attribute(s) ", paste0(names(dat)[w], collapse = ", "),
            " is/are not numeric and has/have all different values and ",
            "assumed to be names: discarded for further analysis"
        )
        dat <- dat[, -w]
        scales <- scales[-w]
    }

    # correct to right data format
    for (col in seq_along(dat)) {
        switch(scales[col],
            "bin" = dat[, col] <- as.logical(dat[, col]),
            "rat" = ,
            "int" = dat[, col] <- as.numeric(dat[, col]),
            "ord" = if (!is.factor(dat[, col])) {
                dat[, col] <- factor(dat[, col], ordered = TRUE)
            } else if (!is.ordered(dat[, col])) {
                dat[, col] <- ordered(dat[, col], levels = levels(dat[, col]))
            },
            "nom" = dat[, col] <- as.factor(dat[, col])
        )
    }

    names(scales) <- colnames(dat)
    values <- lapply(seq_along(dat), function(col) {
        switch(scales[col],
            "bin" = summary(dat[, col])[-1],
            "rat" = ,
            "int" = {
                if (sum(is.na(dat[, col])) > 0) {
                    result <-
                        c(
                            quantile(dat[, col], seq(0, 1, 0.1), na.rm = TRUE),
                            sum(is.na(dat[, col]))
                        )
                    names(result)[12] <- "NA's"
                } else {
                    result <- quantile(dat[, col], seq(0, 1, 0.1), na.rm = TRUE)
                }
                return(result)
            },
            "ord" = ,
            "nom" = {
                result <- summary(dat[, col])
                return(result[result != 0])
            }
        )
    })
    names(values) <- names(scales)
    return(list(dat = dat, scales = scales, values = values))
}

Cgt.fitscale <- function(df) {
    if (!methods::is(df, "DataFrame") && !is.data.frame(df)) {
        stop("Cogito::Cgt.fitscale")
    }
    return(unlist(lapply(df, function(dfc) {
        if (is.numeric(dfc)) {
            if (sum(!is.na(dfc) & dfc < 0) == 0) {
                return("rat")
            } else {
                return("int")
            }
        } else if (is.factor(dfc)) {
            if (is.ordered(dfc)) {
                return("ord")
            } else {
                return("nom")
            }
        } else if (is.character(dfc) &&
            length(unique(dfc)) <= max(10, 0.1 * nrow(df))) {
            return("nom")
        } else if (is.logical(dfc)) {
            return("bin")
        } else {
            return("")
        }
    })))
}

Cgt.guessConditions <- function(names, techs=list(), all=FALSE){
    if (!is.character(names)){
        stop("Cogito::Cgt.guessConditions")
    }
    names <- names[width(names) > 0]
    if (length(names) == 0) {
        stop("Cogito::Cgt.guessConditions")
    }
    posgroups <- unique(unlist(strsplit(names, "[[:punct:]]")))
    if (is.logical(all) && all)
        posgroups <- unique(unlist(lapply(names, function(x) {
            apply(combn(seq_len(nchar(x)), 2), 2, function(pair) {
                substr(x, pair[1], pair[2])
            })
        })))
    wgroups <- lapply(posgroups, function(x) {
        grep(paste0(x, "([[:punct:]]|$)"), names, fixed=FALSE)})
    posgroups <- posgroups[w <- unlist(lapply(wgroups, length)) > 1]
    wgroups <- wgroups[w]
    groups <- unique(wgroups)
    names(groups) <- lapply(groups, function(x) {
        n <- posgroups[unlist(lapply(wgroups, identical, y=x))]
        return(n[which.max(nchar(n))])
    })
    w <- grep("^[[:punct:]]", names(groups))
    todel <- integer()
    if (length(w) > 0) {
        for (g in w) {
            newname <- substr(names(groups)[g], 2, max(nchar(names(groups)[g])))
            if (newname %in% names(groups)) {
                todel <- c(todel, g)
            } else {
                names(groups)[g] <- newname
            }
        }
    }
    if (length(todel) > 0)
        groups <- groups[-todel]
    todel <- integer()
    w <- grep("[[:punct:]]$", names(groups))
    if (length(w) > 0) {
        for (g in w) {
            newname <- substr(names(groups)[g], 1, nchar(names(groups)[g])-1)
            if (newname %in% names(groups)) {
                todel <- c(todel, g)
            } else {
                names(groups)[g] <- newname
            }
        }
    }
    if (length(todel) > 0)
        groups <- groups[-todel]
    groups <- groups[nchar(names(groups)) > 1]
    groups <- groups[unlist(lapply(names(groups), function(name) {
        any(width(strsplit(name, "[[:punct:]]|[[:digit:]]")[[1]]) != 0)
    }))]
    if (length(groups) == 0)
        return(groups)
    if (is.list(techs) && all(unlist(lapply(techs, is.character))))
        groupsn <- groups[!unlist(lapply(groups, function(group) {
            any(unlist(lapply(techs, function(tech) {
                all(names[group] %in% tech)
            })))
        }))]
    if (length(groupsn) > 0)
        groups <- groupsn
    return(groups)
}

Cgt.roundNicely <- function(x, digits = 2) {
    if (!is.numeric(x) || !is.numeric(digits) || digits < 0) {
        stop("Cogito::Cgt.roundNicely")
    }
    if (abs(x) >= 10^(-digits)) {
        result <- paste0(round(x, digits))
    } else {
        result <- formatC(x, format = "e", digits = 1)
    }
    if (result == "0.0e+00") {
        result <- "0"
    }
    names(result) <- names(x)
    return(result)
}

Cgt.rcodestring <- function(code, label = "", echo = NA) {
    if (!is.character(code)) {
        stop("Cogito::Cgt.rcodestring")
    }
    if (!is.character(label)) {
        label <- ""
    }
    if (!is.na(echo) && is.logical(echo)) {
        return(paste0(
            "\n```{r ", label, ", echo=", echo, "}\n",
            code, "\n```\n"
        ))
    } else {
        return(paste0("\n```{r ", label, "}\n", code, "\n```\n"))
    }
}

Cgt.tooltipstring <- function(tooltip, tooltiptext, format = "html") {
    if (!is.character(tooltip) || !is.character(tooltiptext)) {
        stop("Cogito::Cgt::tooltipstring")
    }
    if (is.character(format) && format == "PDF") {
        return(tooltip)
    } else {
        return(paste0(
            "<span class=\"tt\">", tooltip, "<span class=\"ttt\">",
            tooltiptext, "</span></span>"
        ))
    }
}

Cgt.intlinkstring <- function(text, linkto) {
    if ((!is.character(text) && !is.factor(text)) || !is.character(linkto)) {
        stop("Cogito::Cgt.intlinkstring")
    }
    return(paste0("[", text, "](#", linkto, ")"))
}

Cgt.tablestring <- function(data, dotsAfter = integer(), row.names = TRUE,
                            col.widths = NULL) {
    if ((!is.vector(data) && !is.matrix(data) &&
        !methods::is(data, "DataFrame") && !is.data.frame(data)) ||
        !is.integer(dotsAfter)) {
        stop("Cogito::Cgt.tablestring")
    }
    if (is.vector(data)) {
        tmp <- matrix(data, nrow = 1)
        colnames(tmp) <- names(data)
        data <- tmp
    }
    for (i in seq_len(ncol(data))) {
        data[, i] <- as.character(data[, i])
    }
    if (length(dotsAfter) > 0) {
        tmp <- c(
            apply(data, 1, paste, collapse = " | "),
            rep("", length(dotsAfter))
        )
        for (i in seq_along(dotsAfter)) {
            tmp[(dotsAfter[i] + 2):length(tmp)] <-
                tmp[(dotsAfter[i] + 1):(length(tmp) - 1)]
            tmp[dotsAfter[i] + 1] <-
                paste(rep("...", ncol(data)), collapse = " | ")
            dotsAfter <- dotsAfter + 1
        }
    } else {
        tmp <- apply(data, 1, paste, collapse = " | ")
    }
    if (is.null(cnames <- colnames(data))) {
        cnames <- rep(" ", ncol(data))
    }

    if (!is.null(col.widths) && is.vector(col.widths) &&
        is.numeric(col.widths) && length(col.widths) == ncol(data)) {
        tmp2 <- col.widths
    } else {
        tmp2 <- apply(apply(data, 2, width), 2, max)
    }
    tmp2 <- vapply(tmp2, function(x) {
        return(paste0(rep("-", x), collapse = ""))
    }, c(""))
    if (is.null(rownames(data)) || (is.logical(row.names) && !row.names)) {
        result <- paste0(
            paste(cnames, collapse = " | "), "\n",
            paste(tmp2, collapse = " | "), "\n",
            paste(tmp, collapse = "\n")
        )
    } else {
        result <- paste0(
            ". | ", paste(cnames, collapse = " | "), "\n",
            " ", paste(rep("-", max(width(rownames(data)))), collapse = ""),
            " | ", paste(tmp2, collapse = " | "), "\n",
            paste(paste0(rownames(data), " | ", tmp), collapse = "\n")
        )
    }
    return(result)
}



