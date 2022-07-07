# ranges have to be list of GRanges, GRangesList or CompressedGRangesList
# names of ranges have to be in "RRBS|DNA|CVN|RNA|ChIP"
# organism has to be a TxDb or OrganismDb
aggregateRanges <- function(ranges, configfile = NULL, 
                            organism = NULL, referenceRanges = NULL,
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

    if ((is.list(referenceRanges) || is(referenceRanges, "GRangesList")) &&
        (length(referenceRanges) >= 1 && is(referenceRanges[[1]], "GRanges"))) {
        reference <- referenceRanges[[1]]
    } else if ("SYMBOL" %in% AnnotationDbi::columns(genome)) {
        reference <- sort(
            GenomicFeatures::genes(genome, columns = c("SYMBOL"))
        )
    } else {
        reference <- sort(
            GenomicFeatures::genes(genome)
        )
    }

    ranges <- lapply(ranges, function(range) {
        if (methods::is(range, "GRanges")) {
            return(list(range))
        } else {
            return(range)
        }
    })

    maxdist <- configlist$MaxDistToGene

    ## MEAN: distribute RRBS/RNA/DNA/CNV over reference
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
                    reference <- Cgt.aggr(reference, tmp,
                        scale = "numeric",
                        mode = "mean", maxdist
                    )
                } else if (all(unlist(lapply(mcols(tmp), is.ordered)))) {
                    reference <- Cgt.aggr(reference, tmp,
                        scale = "ordered",
                        mode = "median", maxdist
                    )
                } else if (all(unlist(lapply(mcols(tmp), is.factor)))) {
                    reference <- Cgt.aggr(reference, tmp,
                        scale = "factor",
                        mode = "modus", maxdist
                    )
                } else if (all(unlist(lapply(mcols(tmp), is.logical)))) {
                    reference <- Cgt.aggr(reference, tmp,
                        scale = "logical",
                        mode = "any", maxdist
                    )
                } else {
                    stop("Cogito::aggregateRanges no uniform scaling")
                }
            }
        }
    }
    ## MAX: distribute ChIP over reference
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
                    reference <- Cgt.aggr(reference, tmp,
                        scale = "numeric",
                        mode = "max", maxdist
                    )
                } else if (all(unlist(lapply(mcols(tmp), is.ordered)))) {
                    reference <- Cgt.aggr(reference, tmp,
                        scale = "ordered",
                        mode = "max", maxdist
                    )
                } else if (all(unlist(lapply(mcols(tmp), is.factor)))) {
                    reference <- Cgt.aggr(reference, tmp,
                        scale = "factor",
                        mode = "modus", maxdist
                    )
                } else if (all(unlist(lapply(mcols(tmp), is.logical)))) {
                    reference <- Cgt.aggr(reference, tmp,
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
        colnames(mcols(reference))[
            grep(paste0(tech, "."), colnames(mcols(reference)), fixed = TRUE)
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
    shortcolnames <- unlist(lapply(strsplit(colnames(mcols(reference))[-1],
                                            ".", fixed=TRUE), "[[", 2))
    shorttechs <- lapply(configlist$technologies, function(tech) {
        unlist(lapply(strsplit(tech, ".", fixed=TRUE), "[[", 2))})
    groupsIdxs <- Cgt.guessConditions(shortcolnames, shorttechs)
    configlist$conditions <- lapply(groupsIdxs, function(group) {
        colnames(mcols(reference))[-1][group]})

    write(jsonlite::toJSON(configlist, pretty = TRUE),
        append = FALSE,
        file = configfile
    )

    reference <- reference[!apply(mcols(reference), 1,
                                  function(row) all(is.na(row[-1]))), ]
    return(list(genes = reference, config = configlist, name = name))
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

Cgt.guessConditions <- function(names, techs = list(), all = FALSE){
    if (!is.character(names)){
        stop("Cogito::Cgt.guessConditions")
    }
    names <- names[width(names) > 0]
    if (length(names) == 0) {
        stop("Cogito::Cgt.guessConditions")
    }
    posgroups <- unique(unlist(strsplit(names, "[[:punct:]]| ")))
    if (is.logical(all) && all)
        posgroups <- unique(unlist(lapply(names, function(x) {
            apply(combn(seq_len(nchar(x)), 2), 2, function(pair) {
                substr(x, pair[1], pair[2])
            })
        })))
    wgroups <- lapply(posgroups, function(x) {
        grep(paste0(x, "([[:punct:]]| |$)"), names, fixed=FALSE)})
    posgroups <- posgroups[w <- unlist(lapply(wgroups, length)) > 1]
    wgroups <- wgroups[w]
    groups <- unique(wgroups)
    names(groups) <- lapply(groups, function(x) {
        n <- posgroups[unlist(lapply(wgroups, identical, y=x))]
        return(n[which.max(nchar(n))])
    })
    w <- grep("^([[:punct:]]| )", names(groups))
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
    w <- grep("([[:punct:]]| )$", names(groups))
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
        any(width(strsplit(name, "[[:punct:]]|[[:digit:]]| ")[[1]]) != 0)
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
