# PRA2003 - Analyze all 10 sub-sample files
#
# Replaces having 10 near-identical scripts (one per data file) with a
# single reusable function, looped over all 10 filenames. Each file is
# read in memory-safe chunks (never loading the whole ~800MB file at
# once), and for each of the 12 bacterial strains we compute:
#   average per event = total count / number of valid events
#   uncertainty        = sqrt(total count) / number of valid events
#     (Poisson statistics - see week 3)
# Events with 0 particles (e.g. header "77 0") are excluded from the
# event count used for the average, since nothing was observed in them.
#
# Writes subsample_results.csv (one row per file x strain), which
# Week4deliverable.R reads in to compute the combined result.

# ---- The 12 bacterial strains of interest, in a fixed column order ----
strainID <- c(211, -211, 321, -321, 2212, -2212,
              3122, -3122, 3312, -3312, 3334, -3334)
strainName <- c("E. coli WT", "E. coli mutant",
                "Bacillus subtilis WT", "Bacillus subtilis mutant",
                "Pseudomonas aeruginosa WT", "Pseudomonas aeruginosa antibiotic-resistant",
                "Streptococcus pneumoniae", "Capsule-deficient S. pneumoniae",
                "Mycobacterium tuberculosis", "Drug-resistant M. tuberculosis",
                "Salmonella enterica", "Salmonella mutant")
nStrains <- length(strainID)

# ---- Analyze a single data file, returning per-strain average + uncertainty ----
analyseFile <- function(filepath, chunkSize = 1000000) {
  if (!file.exists(filepath)) {
    stop("Cannot find '", filepath, "' in the current working directory (",
         getwd(), "). Set your working directory to the folder containing ",
         "the data files, or edit the filenames used below.")
  }

  totalCount <- rep(0, nStrains)
  nEvents <- 0

  con <- file(filepath, "r")
  on.exit(close(con))

  repeat {
    chunk <- readLines(con, n = chunkSize)
    if (length(chunk) == 0) {
      break  # reached end of file
    }

    # Count fields per line via number of spaces (+1). fixed=TRUE/useBytes=TRUE
    # skip regex-engine and locale overhead, which matters a lot at 25M lines.
    nFields <- lengths(gregexpr(" ", chunk, fixed = TRUE, useBytes = TRUE)) + 1

    isHeader <- (nFields == 2)
    isData   <- (nFields == 4)

    # Some events have 0 particles (e.g. "77 0") - excluded from nEvents
    headerLines <- chunk[isHeader]
    if (length(headerLines) > 0) {
      nParticlesInEvent <- as.integer(
        sub("^\\S+ ", "", headerLines, perl = TRUE, useBytes = TRUE)
      )
      nEvents <- nEvents + sum(nParticlesInEvent > 0)
    }

    dataLines <- chunk[isData]
    if (length(dataLines) > 0) {
      # bacterial ID is always the last whitespace-separated field
      ids <- as.integer(
        sub("^.* ", "", dataLines, perl = TRUE, useBytes = TRUE)
      )

      for (i in 1:nStrains) {
        totalCount[i] <- totalCount[i] + sum(ids == strainID[i])
      }
    }
  }

  if (nEvents == 0) {
    stop("No events found in '", filepath, "' - check the file format.")
  }

  list(
    nEvents = nEvents,
    totalCount = totalCount,
    avgPerEvent = totalCount / nEvents,
    uncertainty = sqrt(totalCount) / nEvents
  )
}

# ---- Loop over all 10 sub-sample files ----
nFiles <- 10
filenames <- paste0("output-Set", 1:nFiles, ".txt")

subsampleRows <- vector("list", nFiles * nStrains)
row <- 1

for (setNum in 1:nFiles) {
  cat("Analyzing", filenames[setNum], "...\n")
  result <- analyseFile(filenames[setNum])

  for (i in 1:nStrains) {
    subsampleRows[[row]] <- data.frame(
      Set = setNum,
      ID = strainID[i],
      Strain = strainName[i],
      AveragePerEvent = round(result$avgPerEvent[i], 5),
      Uncertainty = round(result$uncertainty[i], 5),
      nEvents = result$nEvents
    )
    row <- row + 1
  }
}

subsampleResults <- do.call(rbind, subsampleRows)

options(scipen = 999)  # avoid scientific notation (e.g. 5e-05) in the CSV
write.csv(subsampleResults, "subsample_results.csv", row.names = FALSE)

cat("\nWrote subsample_results.csv (", nrow(subsampleResults), "rows).\n")
