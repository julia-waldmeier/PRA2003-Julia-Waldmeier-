# PRA2003 - Week 4 deliverable
#
# Combines the 10 independent sub-sample results (from subsample_results.csv,
# produced by analyzeAllSubsamples.R running the shared analysis function
# across all 10 ~500K-event data files) into a single best estimate for
# each of the 12 bacterial strains, using the sub-sampling method from the
# lecture:
#   1. Combine the 10 sub-sample averages into one weighted average
#      (weight = 1 / uncertainty^2, so more precise sub-samples count more)
#   2. Use the spread (standard deviation) of the 10 sub-sample averages
#      directly as the statistical uncertainty on the combined result:
#          uncertainty_combined = sd(sub-sample averages)

subsampleFile <- "subsample_results.csv"
if (!file.exists(subsampleFile)) {
  stop("Cannot find '", subsampleFile, "'. Run analyzeAllSubsamples.R first ",
       "to generate it from the 10 raw data files.")
}

subsampleResults <- read.csv(subsampleFile, stringsAsFactors = FALSE)

# ---- The 12 bacterial strains, in the order they first appear in the CSV ----
strainOrder <- unique(subsampleResults$ID)
strainID <- strainOrder
strainName <- subsampleResults$Strain[match(strainID, subsampleResults$ID)]
nStrains <- length(strainID)

nFiles <- length(unique(subsampleResults$Set))
nEventsPerFile <- subsampleResults$nEvents[match(unique(subsampleResults$Set),
                                                  subsampleResults$Set)]

# ---- Combine the 10 sub-samples, per strain ----
weightedAverage <- numeric(nStrains)
subSamplingUncertainty <- numeric(nStrains)

for (i in 1:nStrains) {
  rows <- subsampleResults[subsampleResults$ID == strainID[i], ]
  rows <- rows[order(rows$Set), ]

  x <- rows$AveragePerEvent      # the 10 sub-sample averages for this strain
  sigma <- rows$Uncertainty      # their individual (Poisson) uncertainties

  weights <- 1 / sigma^2
  weightedAverage[i] <- sum(x * weights) / sum(weights)

  # Spread (standard deviation) of the 10 sub-sample averages, used
  # directly as the statistical uncertainty of the combined measurement
  subSamplingUncertainty[i] <- sd(x)
}

totalEvents <- sum(nEventsPerFile)

results <- data.frame(
  ID = strainID,
  Strain = strainName,
  WeightedAverage = round(weightedAverage, 5),
  Uncertainty = round(subSamplingUncertainty, 5)
)

cat("Combined", totalEvents, "events across", nFiles, "sub-samples.\n\n")
print(results)

# ---- Write out the final combined results for reproducibility ----
options(scipen = 999)  # avoid scientific notation (e.g. 5e-05) in the CSV
write.csv(results, "final_results.csv", row.names = FALSE)

cat("\nWrote final_results.csv (", nrow(results), "rows).\n")
