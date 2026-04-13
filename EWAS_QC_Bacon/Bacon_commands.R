# from https://www.bioconductor.org/packages/release/bioc/vignettes/bacon/inst/doc/bacon.html#1_Introduction

#Running Bacon
bc <- bacon(teststatistics = df$tstat,
              effectsizes = df$beta,
              standarderrors = df$SE)

# Bias estimated by Bacon
bias(bc)

# Inflation estimated by Bacon
inflation(bc)

# See original/corrected p-values
pval(bc,object = )

# Bacon empirical null distribution fit
fit(bc)

# Bacon estimates for empirical null distribution
estimates(bc)

# Bacon corrected vs uncorrected p-values QQ-plot
plot(bc, type="qq")
