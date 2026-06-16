# SDRF_FIXED creation 
# *due to differences in naming of technical replicates between our generated protein groups and the original SDRF
# **the changes in the SDRF were subsequently validated

library(dplyr)

sdrf_design <- read.delim("sdrf.tsv", check.names = FALSE)

# 1 Extract correct design columns
sdrf_design <- sdrf_design %>%
  mutate(
    genotype = `characteristics[genotype]`,
    polstate = `characteristics[phase]`
  )

# 2 Safety check
stopifnot(all(c("genotype", "polstate") %in% colnames(sdrf_design)))

# 3 Build sample names
sdrf_design <- sdrf_design %>%
  group_by(genotype, polstate) %>%
  mutate(rep = row_number()) %>%
  ungroup() %>%
  mutate(sample = paste(genotype, polstate, rep, sep = "_"))

# 4 Clean output
sdrf_design <- dplyr::select(
  sdrf_design,
  sample,
  genotype,
  polstate,
  everything()
)

# 5 Save new SDRF
write.table(
  sdrf_design,
  file = "sdrf_FIXED.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("SDRF fixed and saved: sdrf_FIXED.tsv")

