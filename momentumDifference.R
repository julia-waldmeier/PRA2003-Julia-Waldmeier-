# PRA2003 - Week 5: WT - mutant difference as a function of momentum
#
# For each of the 6 WT/mutant strain pairs, checks whether the difference
# between the WT/normal form and the mutant/resistant form depends on the
# momentum of the bacteria (p = sqrt(px^2 + py^2 + pz^2)).
#
# Builds on momentumAsymmetry.R (Week 4): the counting is the same, but
# instead of a ratio it now uses the DIFFERENCE per event, the same quantity
# as in Q2 (asymmetryTest.R), and adds a fit of a constant and a plot.
#
# Method:
#   1. Every occurrence of one of the 12 target strains has its momentum
#      computed and is sorted into one of nBins momentum bins. The bin
#      boundaries are the quantiles of Set 1 (so each bin holds about 20% of
#      the bacteria) and are reused for every file. The first and last bins
#      are open-ended, so every bacterium is counted.
#   2. For each file, pair and bin: count the WT and the mutant.
#   3. For each file, pair and bin: the difference per event (like Q1/Q2),
#        D_file = countWT / nEvents - countMutant / nEvents
#      where nEvents = number of non-empty events in that file.
#   4. Per pair and bin, the SUB-SAMPLING method (same as Q1/Q2):
#        point       = mean of the per-file values of D
#        uncertainty = standard deviation of the per-file values of D
#        significance = point / uncertainty   (significant if >= 3 sigma)
#   5. Per pair, a weighted fit of a CONSTANT through the nBins points,
#      i.e. a flat line D = c, to test for an asymmetry that is the same at
#      every momentum. With weight w = 1 / uncertainty^2 the least-squares
#      result is the weighted mean of the points:
#        c = sum(w * y) / sum(w),   uncertainty = 1 / sqrt(sum(w))
#      c / uncertainty = significance of a momentum-independent asymmetry.
#      chi2 = sum(((y - c) / uncertainty)^2) with ndf = nBins - 1 tells
#      whether a constant describes the points (chi2/ndf around 1 = good).
#
# Set 10 has a very different momentum distribution from Sets 1-9 (almost
# no bacteria above the highest bin edge, 58% in the 2.68-5.87 bin, while
# Sets 1-9 have 20% in every bin), so it is not a repeat of the same
# measurement. The main result therefore uses Sets 1-9 (setting
# excludedSets). The same analysis with all 10 sets is also written, for
# comparison.
#
# HOW TO RUN:
#   Run analyzeAllSubsamples.R first (for the number of events per file in
#   subsample_results.csv). Then: Rscript momentumDifference.R
#   The counts per file are cached in momentum_perfile_results.csv and the
#   bin boundaries in momentum_bin_edges.csv. If both exist, the raw data is
#   not read again (a few seconds). Delete them to redo the counting
#   (about 20 minutes).
#
# Writes (main result, Sets 1-9):
#   momentum_difference_results.csv - point +- uncertainty per pair and bin
#   momentum_difference_fit.csv     - fitted constant +- uncertainty per pair
#   momentum_difference.png         - the points, error bars and fitted constant
# and the same three files with "_all_sets" in the name for all 10 sets.

# ---- Settings ----
nBins <- 5
excludedSets <- c(10)   # left out of the main result (see above); use c() to keep all

filenames <- paste0("output-Set", 1:10, ".txt")
nFiles <- length(filenames)

pairPositiveID <- c(211, 321, 2212, 3122, 3312, 3334)
pairNegativeID <- c(-211, -321, -2212, -3122, -3312, -3334)
pairName <- c("E. coli", "Bacillus subtilis", "Pseudomonas aeruginosa",
              "Streptococcus pneumoniae", "Mycobacterium tuberculosis",
              "Salmonella")
targetIDs <- c(pairPositiveID, pairNegativeID)
nPairs <- length(pairPositiveID)

chunkSize <- 1000000

countsFile <- "momentum_perfile_results.csv"
edgesFile <- "momentum_bin_edges.csv"
eventsFile <- "subsample_results.csv"

if (!file.exists(eventsFile)) {
  stop("Cannot find '", eventsFile, "'. Run analyzeAllSubsamples.R first.")
}

# ---- Reads one file, returning the momentum and ID of every occurrence
# of a target strain in it (same function as in momentumAsymmetry.R) ----
extractMomentumAndID <- function(filepath) {
  con <- file(filepath, "r")
  on.exit(close(con))

  momentumList <- list()
  idList <- list()
  chunkIndex <- 0

  repeat {
    chunk <- readLines(con, n = chunkSize)
    if (length(chunk) == 0) {
      break
    }
    chunkIndex <- chunkIndex + 1

    nFields <- lengths(gregexpr(" ", chunk, fixed = TRUE, useBytes = TRUE)) + 1
    isData <- (nFields == 4)
    dataLines <- chunk[isData]

    if (length(dataLines) > 0) {
      ids <- as.integer(sub("^.* ", "", dataLines, perl = TRUE, useBytes = TRUE))
      keep <- ids %in% targetIDs

      if (any(keep)) {
        keptLines <- dataLines[keep]
        parsed <- read.table(text = keptLines, colClasses = "numeric")
        p <- sqrt(parsed[[1]]^2 + parsed[[2]]^2 + parsed[[3]]^2)

        momentumList[[chunkIndex]] <- p
        idList[[chunkIndex]] <- ids[keep]
      }
    }
  }

  list(momentum = unlist(momentumList), id = unlist(idList))
}

# ---- Steps 1-2: bin boundaries and counts per file ----
for (f in filenames) {
  if (!file.exists(countsFile) && !file.exists(f)) {
    stop("Cannot find '", f, "' in the current working directory (", getwd(), ").")
  }
}

if (file.exists(countsFile) && file.exists(edgesFile)) {
  cat("Using the saved counts in", countsFile, "\n\n")
  breaks <- read.csv(edgesFile)$Edge

} else if (file.exists(countsFile)) {
  # Counts are saved, only the bin boundaries are missing: they come from
  # Set 1 alone, so only Set 1 has to be read again
  cat("Recomputing the bin boundaries from", filenames[1], "...\n")
  set1 <- extractMomentumAndID(filenames[1])
  breaks <- quantile(set1$momentum, probs = seq(0, 1, length.out = nBins + 1), names = FALSE)
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf
  write.csv(data.frame(Edge = breaks), edgesFile, row.names = FALSE)

} else {
  # countPerFile[file, pair, bin, 1] = WT count, [,,,2] = mutant count
  countPerFile <- array(0, dim = c(nFiles, nPairs, nBins, 2))
  breaks <- NULL

  for (fileIndex in 1:nFiles) {
    cat("Processing", filenames[fileIndex], "...\n")
    fileData <- extractMomentumAndID(filenames[fileIndex])

    # Set 1 defines the bin boundaries (quantiles), reused for every file
    if (fileIndex == 1) {
      breaks <- quantile(fileData$momentum, probs = seq(0, 1, length.out = nBins + 1),
                         names = FALSE)
      breaks[1] <- -Inf
      breaks[length(breaks)] <- Inf
    }

    bin <- as.integer(cut(fileData$momentum, breaks = breaks, include.lowest = TRUE))

    for (i in 1:nPairs) {
      for (b in 1:nBins) {
        inBin <- (bin == b)
        countPerFile[fileIndex, i, b, 1] <- sum(fileData$id == pairPositiveID[i] & inBin)
        countPerFile[fileIndex, i, b, 2] <- sum(fileData$id == pairNegativeID[i] & inBin)
      }
    }
  }

  perFileRows <- data.frame()
  for (fileIndex in 1:nFiles) {
    for (i in 1:nPairs) {
      for (b in 1:nBins) {
        perFileRows <- rbind(perFileRows, data.frame(
          Set = fileIndex,
          Strain = pairName[i],
          MomentumBin = b,
          CountWT = countPerFile[fileIndex, i, b, 1],
          CountMutant = countPerFile[fileIndex, i, b, 2]
        ))
      }
    }
  }
  options(scipen = 999)
  write.csv(perFileRows, countsFile, row.names = FALSE)
  write.csv(data.frame(Edge = breaks), edgesFile, row.names = FALSE)
  cat("\nWrote", countsFile, "and", edgesFile, "\n")
}

counts <- read.csv(countsFile)

# Labels for the momentum bins, e.g. "0.65-1.32" and ">5.87"
binLabel <- character(nBins)
for (b in 1:nBins) {
  if (b == 1) {
    binLabel[b] <- paste0("<", round(breaks[2], 2))
  } else if (b == nBins) {
    binLabel[b] <- paste0(">", round(breaks[nBins], 2))
  } else {
    binLabel[b] <- paste0(round(breaks[b], 2), "-", round(breaks[b + 1], 2))
  }
}

# ---- Number of non-empty events in each file (from the Week 4 analysis) ----
subsampleResults <- read.csv(eventsFile)
nEventsPerFile <- numeric(nFiles)
for (fileIndex in 1:nFiles) {
  nEventsPerFile[fileIndex] <- subsampleResults$nEvents[subsampleResults$Set == fileIndex][1]
}

# ---- Steps 3-5 and the plot, for a chosen group of sets ----
analyseSets <- function(setsUsed, nameSuffix) {
  nUsed <- length(setsUsed)

  # Steps 3-4: difference per event per file, then mean and SD over the
  # files, per pair and bin
  results <- data.frame()

  for (i in 1:nPairs) {
    for (b in 1:nBins) {
      perFileDiff <- numeric(nUsed)

      for (k in 1:nUsed) {
        fileIndex <- setsUsed[k]
        row <- counts[counts$Set == fileIndex & counts$Strain == pairName[i] &
                        counts$MomentumBin == b, ]
        averageWT <- row$CountWT / nEventsPerFile[fileIndex]
        averageMutant <- row$CountMutant / nEventsPerFile[fileIndex]
        perFileDiff[k] <- averageWT - averageMutant
      }

      difference <- mean(perFileDiff)
      uncertainty <- sd(perFileDiff)

      results <- rbind(results, data.frame(
        Strain = pairName[i],
        MomentumBin = b,
        MomentumRange = binLabel[b],
        Difference = signif(difference, 4),
        Uncertainty = signif(uncertainty, 4),
        Significance = round(difference / uncertainty, 2)
      ))
    }
  }

  # Step 5: weighted fit of a constant per pair
  fitResults <- data.frame()

  for (i in 1:nPairs) {
    pairRows <- results[results$Strain == pairName[i], ]
    y <- pairRows$Difference
    w <- 1 / pairRows$Uncertainty^2   # precise points count more

    constant <- sum(w * y) / sum(w)
    constantUncertainty <- 1 / sqrt(sum(w))
    chi2 <- sum(((y - constant) / pairRows$Uncertainty)^2)
    ndf <- nBins - 1

    fitResults <- rbind(fitResults, data.frame(
      Strain = pairName[i],
      Constant = signif(constant, 4),
      ConstantUncertainty = signif(constantUncertainty, 4),
      Significance = round(constant / constantUncertainty, 2),
      Chi2 = round(chi2, 2),
      ndf = ndf,
      pValue = round(pchisq(chi2, ndf, lower.tail = FALSE), 2)
    ))
  }

  setsText <- paste0("sets ", paste(setsUsed, collapse = ","))
  options(width = 200)
  cat("\n==== Using", setsText, "====\n\n")
  cat("Difference per event (WT - mutant) per momentum bin, mean +- SD of the files:\n\n")
  print(results, row.names = FALSE)
  cat("\nFit of a constant (difference = c, the same in every momentum bin):\n\n")
  print(fitResults, row.names = FALSE)

  write.csv(results, paste0("momentum_difference_results", nameSuffix, ".csv"), row.names = FALSE)
  write.csv(fitResults, paste0("momentum_difference_fit", nameSuffix, ".csv"), row.names = FALSE)

  # Plot: one panel per pair, points +- SD and the fitted constant
  png(paste0("momentum_difference", nameSuffix, ".png"), width = 1800, height = 1100, res = 140)
  par(mfrow = c(2, 3), mar = c(6, 5.5, 4.5, 1), las = 1, mgp = c(4, 0.7, 0),
      oma = c(0, 0, 2, 0))

  for (i in 1:nPairs) {
    pairRows <- results[results$Strain == pairName[i], ]
    fit <- fitResults[fitResults$Strain == pairName[i], ]
    x <- pairRows$MomentumBin
    y <- pairRows$Difference
    low <- y - pairRows$Uncertainty
    high <- y + pairRows$Uncertainty

    plot(x, y, ylim = range(c(low, high, 0)), xlim = c(0.7, nBins + 0.3),
         xaxt = "n", pch = 19, col = "darkblue",
         xlab = "", ylab = "Difference per event (WT - mutant)",
         main = pairName[i], font.main = 3)
    mtext(sprintf("constant = %.3g +/- %.3g (%.1f sigma), chi2/ndf = %.1f/%d",
                  fit$Constant, fit$ConstantUncertainty, fit$Significance, fit$Chi2, fit$ndf),
          side = 3, line = 0.4, cex = 0.75)
    axis(1, at = x, labels = binLabel, cex.axis = 0.8)
    mtext("Momentum bin |p| (10^-20 kg m/s)", side = 1, line = 3, cex = 0.75)

    # fitted constant with its +- 1 sigma band
    rect(0, fit$Constant - fit$ConstantUncertainty, nBins + 1,
         fit$Constant + fit$ConstantUncertainty, col = rgb(1, 0, 0, 0.15), border = NA)
    abline(h = fit$Constant, col = "red", lwd = 2)
    abline(h = 0, lty = 2, col = "grey")                       # no difference
    arrows(x, low, x, high, angle = 90, code = 3, length = 0.04) # +- SD error bars
    points(x, y, pch = 19, col = "darkblue")
  }
  mtext(paste0("Difference per event vs momentum (", setsText, ")"),
        outer = TRUE, side = 3, line = 0.5, font = 2)

  dev.off()
  cat("\nWrote momentum_difference_results", nameSuffix, ".csv, momentum_difference_fit",
      nameSuffix, ".csv and momentum_difference", nameSuffix, ".png\n", sep = "")
}

# Main result: without the excluded set(s)
analyseSets(setdiff(1:nFiles, excludedSets), "")

# Comparison: all 10 sets
analyseSets(1:nFiles, "_all_sets")
