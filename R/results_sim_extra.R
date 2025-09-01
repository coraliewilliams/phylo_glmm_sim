# This script has the results for the simulated models with one measure per species (n_rep=1)
# using the model in equation 2 in the main text

################### Load results #############################

library(ggplot2);library(devtools);library(cowplot);library(ggdark);theme_set(dark_theme_bw())
library(ggdist);library(tidyverse);library(dplyr);library(patchwork); library(scales);
library(latex2exp); library(xtable)

#---
# load results set 1b (one measure per species)
load("~/PhD/3_PhyloTMB/sims/sim_set1b_results_pglmm.RDATA")
res1.dat <- dat

# load results set 2b (one measure per species)
load("~/PhD/3_PhyloTMB/sims/sim_set2b_results_pglmm.RDATA")
res2.dat <- dat



################# 1. Convergence #########################


# Check convergence %
#table(res1.dat$model, res1.dat$convergence)/24000*100
#table(res2.dat$model, res2.dat$convergence)/3000*100

# Check percentage ESS>400
#table(res1.dat$model, res1.dat$ESS>400)/24000*100
#table(res2.dat$model, res2.dat$ESS>400)/3000*100


### Summary tables ---
tab.conv <- as.data.frame(table(res1.dat$convergence, res1.dat$model, res1.dat$species_size)/8000*100)
res1.dat_conv_summary <- tab.conv |> 
  filter(Var1 == TRUE) |> 
  dplyr::select(Var2, Var3, Freq) |> 
  pivot_wider(names_from = Var3, values_from = Freq) |> 
  mutate(model = factor(Var2, levels=c("brms", "MCMCglmm", "INLA", "glmmTMB", "phyr"))) |> 
  arrange(model)

tab.conv <- as.data.frame(table(res2.dat$convergence, res2.dat$model, res2.dat$species_size)/1000*100)
res2.dat_conv_summary <- tab.conv |> 
  filter(Var1 == TRUE) |> 
  dplyr::select(Var2, Var3, Freq) |> 
  pivot_wider(names_from = Var3, values_from = Freq) |> 
  mutate(model = factor(Var2, levels=c("brms", "MCMCglmm", "INLA", "glmmTMB", "phyr"))) |> 
  arrange(model)

# # merge convergence info
# conv.tab <- full_join(res1.dat_conv_summary, res2.dat_conv_summary, by = "model") %>%
#   rename(Package = model) %>%
#   select(Package, `25`, `50`, `100`, `200`, `400`, `800`)


# # Create summary tables of convergence % for Supp. information
# print(xtable(conv.tab, digits = 1),
#       include.rownames = FALSE)

# All MCMCglmm models have ESS>400
# check ESS
table(res1.dat$ESS>400, res1.dat$model)
table(res2.dat$ESS>400, res2.dat$model)


# Remove all models of seed if one didn't converge
res1 <- res1.dat |> 
  group_by(seed, scenario, sigma2.p) |> 
  filter(all(convergence == TRUE, na.rm = FALSE)) |> 
  ungroup()

res2 <- res2.dat |> 
  group_by(seed, scenario, sigma2.p) |> 
  filter(all(convergence == TRUE, na.rm = FALSE)) |> 
  ungroup()



table(res1$species_size)/5/64000
table(res2$species_size)/3/8000

table(res1$species_size)
table(res2$species_size)

# assess model errors
table(is.na(res1$model_error), res1$model)
# assess model warnings
table(is.na(res1.dat$model_warning), res1.dat$model)
table(is.na(res2.dat$model_warning), res2.dat$model)


# all glmmTMB models have the same warning
table(res1$model_warning[which(res1$model=="glmmTMB")])
table(res2$model_warning[which(res1$model=="glmmTMB")])


table(res1$s2_phylo > 4, res1$model)
table(res1$s2_resid > 4, res1$model)

table(res1$s2_phylo_bias > 4, res1$model)
table(res1$s2_resid_bias > 4, res1$model)


hist(res1$s2_phylo_bias[which(res1$model=="INLA")], breaks=100)
hist(res1$s2_phylo_bias[which(res1$model=="INLA" & res1$s2_phylo_bias < 100)], breaks=100)




