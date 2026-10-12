# PRA2003 - Week 4: asymmetry test (WT vs mutant/resistant strains)
#
# For each of the 6 strain pairs, computes the difference WITHIN each of
# the 10 files first (rather than subtracting the two already-combined
# final averages), since WT and mutant counts are correlated within a
# file (a file with more bacteria overall has more of both forms
# together) - subtracting per-file cancels that shared noise before
# combining.
#
#   For each file i:  diff_i = average_i(X) - average_i(-X)
#   Reported difference = mean of the 10 diff_i values
#   Reported uncertainty = standard deviation of the 10 diff_i values
#     (same sub-sampling convention as Week4deliverable.R, no division
#     by sqrt(10))
#   Significance = |difference| / uncertainty, in units of sigma
#
# Reads subsample_results.csv (produced by analyzeAllSubsamples.R).
# Writes asymmetry_results.csv.

subsampleFile <- "subsample_results.csv"
if (!file.exists(subsampleFile)) {
  stop("Cannot find '", subsampleFile, "'. Run analyzeAllSubsamples.R first ",
       "to generate it from the 10 raw data files.")
}

subsampleResults <- read.csv(subsampleFile, stringsAsFactors = FALSE)

# ---- The 6 wild-type/mutant strain pairs of interest ----
pairPositiveID <- c(211, 321, 2212, 3122, 3312, 3334)
pairNegativeID <- c(-211, -321, -2212, -3122, -3312, -3334)
pairName <- c("E. coli (WT vs mutant)",
              "Bacillus subtilis (WT vs mutant)",
              "Pseudomonas aeruginosa (WT vs resistant)",
              "Streptococcus pneumoniae (normal vs capsule-def.)",
              "Mycobacterium tuberculosis (WT vs drug-resistant)",
              "Salmonella (WT vs mutant)")
nPairs <- length(pairPositiveID)

difference <- numeric(nPairs)
uncertainty <- numeric(nPairs)
significance <- numeric(nPairs)

for (i in 1:nPairs) {
  posRows <- subsampleResults[subsampleResults$ID == pairPositiveID[i], ]
  negRows <- subsampleResults[subsampleResults$ID == pairNegativeID[i], ]
  posRows <- posRows[order(posRows$Set), ]
  negRows <- negRows[order(negRows$Set), ]

  # difference computed WITHIN each file first (one value per file)
  perFileDiff <- posRows$AveragePerEvent - negRows$AveragePerEvent

  difference[i] <- mean(perFileDiff)
  uncertainty[i] <- sd(perFileDiff)
  significance[i] <- abs(difference[i]) / uncertainty[i]
}

results <- data.frame(
  Pair = pairName,
  Difference = round(difference, 5),
  Uncertainty = round(uncertainty, 5),
  Significance = round(significance, 2),
  Significant = significance > 3
)

print(results)

options(scipen = 999)
write.csv(results, "asymmetry_results.csv", row.names = FALSE)

cat("\nWrote asymmetry_results.csv (", nrow(results), "rows).\n")
