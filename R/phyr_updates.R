# Updates for phyr variance components
library(dplyr)

# sim iteration
iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

# set up 
res.path <- paste0("results/a/res_", iter, ".RDATA")
dat.path <- paste0("data/a/simdat_", iter, ".RDATA")

# load result and data of simulation iteration
load(res.path)
load(dat.path)

# for checks
message("iter = ", iter)
message("wd   = ", getwd())
message("res  = ", res.path)
message("dat  = ", dat.path)


# get s2r for species and phylo effects
s2r.species <- result$s2_sp[which(result$model=="phyr")] * result$s2_resid[which(result$model=="phyr")]
s2r.phylo <- result$s2_phylo[which(result$model=="phyr")] * result$s2_resid[which(result$model=="phyr")]

# scale factor
V <- phylo.mat
V_scaled <- V / max(V)
scale_factor <- max(V) * exp(determinant(V_scaled)$modulus[1] / nrow(V))

# get standardise s2 phylo
s2_phylo_std <- s2r.phylo / scale_factor

# update s2_sp and s2_phylo in result only for phyr and corresponding bias
s2_sp_bias = s2r.species - result$sigma2.s[1]
s2_phylo_bias = s2_phylo_std - result$sigma2.p[1]

# new phyr row with correct variance components for sp and phylo
phyr_result <- result[which(result$model=="phyr"), ]
phyr_result$s2_sp <- s2r.species
phyr_result$s2_phylo <- s2_phylo_std
phyr_result$s2_sp_bias <- s2_sp_bias
phyr_result$s2_phylo_bias <- s2_phylo_bias

# #checks
# phyr_old <- result[which(result$model=="phyr"), ]
# testthat::expect_equal(phyr_old, phyr_result)

# remove current row "phyr" model and add new row with updated variance components
result <- result |> 
  filter(model != "phyr") |> 
  bind_rows(phyr_result)


# Save results as RDATA
save(list = "result", file = paste0("results/a/res_", iter, ".RDATA"))
