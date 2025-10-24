# This script has the results for the simulated models with repeated measures per species (n_rep>1)
# using the model in equation 3 in the main text

################### Load results #############################

library(ggplot2);library(devtools);library(cowplot);library(ggdark);theme_set(dark_theme_bw())
library(ggdist);library(tidyverse);library(dplyr);library(patchwork); library(scales);
library(latex2exp); library(xtable); library(grid)

#---
# load results set 1a (repeated measures)
load("~/PhD/3_PhyloTMB/sims/sim_set1a_results_pglmm.RDATA")
res1.dat <- dat

# load results set 2a (repeated measures)
load("~/PhD/3_PhyloTMB/sims/sim_set2a_results_pglmm.RDATA")
res2.dat <- dat




################# 1. Convergence #########################


# Check convergence %
#table(res1.dat$convergence, res1.dat$model)
#table(res2.dat$convergence, res2.dat$model)
#table(res1.dat$model, res1.dat$convergence)/192000*100
#table(res2.dat$model, res2.dat$convergence)/24000*100

# Check percentage ESS>400
#table(res1.dat$model, res1.dat$ESS>400)/192000*100
#table(res2.dat$model, res2.dat$ESS>400)/24000*100


#### Set specific convergence check for INLA models if s2_sp and s2_phylo are > 4 --> alternatively change to 3
res1.dat <- res1.dat |> 
  mutate(convergence = ifelse(model == "INLA" & (s2_sp > 4 | s2_phylo > 4), FALSE, convergence))

res2.dat <- res2.dat |>
  mutate(convergence = ifelse(model == "INLA" & (s2_sp > 4 | s2_phylo > 4), FALSE, convergence))


### Summary tables ---
tab.conv <- as.data.frame(table(res1.dat$convergence, res1.dat$model, res1.dat$species_size)/64000*100)
res1.dat_conv_summary <- tab.conv |> 
  filter(Var1 == TRUE) |> 
  dplyr::select(Var2, Var3, Freq) |> 
  pivot_wider(names_from = Var3, values_from = Freq) |> 
  mutate(model = factor(Var2, levels=c("brms", "MCMCglmm", "INLA", "glmmTMB", "phyr"))) |> 
  arrange(model)

tab.conv <- as.data.frame(table(res2.dat$convergence, res2.dat$model, res2.dat$species_size)/8000*100)
res2.dat_conv_summary <- tab.conv |> 
  filter(Var1 == TRUE) |> 
  dplyr::select(Var2, Var3, Freq) |> 
  pivot_wider(names_from = Var3, values_from = Freq) |> 
  mutate(model = factor(Var2, levels=c("brms", "MCMCglmm", "INLA", "glmmTMB", "phyr"))) |> 
  arrange(model)

# merge convergence info
conv.tab <- full_join(res1.dat_conv_summary, res2.dat_conv_summary, by = "model") %>%
  rename(Package = model) %>%
  select(Package, `25`, `50`, `100`, `200`, `400`, `800`)


# Create summary tables of convergence % for Supp. information
print(xtable(conv.tab, digits = 1),
      include.rownames = FALSE)
#check ESS
table(res1.dat$ESS>400, res1.dat$model)
table(res2.dat$ESS>400, res2.dat$model)




### Filter out all sim iterations with one failed model ---

# Update convergence variable if ESS <400 only for MCMCglmm and brms
res1.dat <- res1.dat |> 
  mutate(convergence = ifelse(model %in% c("MCMCglmm", "brms") & ESS < 400, FALSE, convergence))
res2.dat <- res2.dat |> 
  mutate(convergence = ifelse(model %in% c("MCMCglmm", "brms") & ESS < 400, FALSE, convergence))


# Remove all models of seed if one didn't converge
res1 <- res1.dat |> 
  group_by(seed, scenario, sigma2.p, sigma2.s) |> 
  filter(all(convergence == TRUE, na.rm = FALSE)) |> 
  ungroup()

res2 <- res2.dat |> 
  group_by(seed, scenario, sigma2.p, sigma2.s) |> 
  filter(all(convergence == TRUE, na.rm = FALSE)) |> 
  ungroup()


# Now check ESS
#table(res1$ESS<400, res1$model)
#table(res2$ESS<400, res2$model)

#### Get number of simulation per package (counts + percentages)
table(res1$species_size)
table(res2$species_size)

table(res1$species_size)/5/64000
table(res2$species_size)/3/8000


#table(res1$species_size, res1$nrep)/5/16000

# # assess model errors
# table(is.na(res1$model_error), res1$model)
# # assess model warnings
# table(is.na(res1.dat$model_warning), res1.dat$model)
# table(is.na(res2.dat$model_warning), res2.dat$model)
# 
# # all INLA models have the same warning
# table(res1$model_warning[which(res1$model=="INLA")])
# table(res2$model_warning[which(res1$model=="INLA")])
# 
# table(res1$s2_sp > 4, res1$model)
# table(res1$s2_phylo > 4, res1$model)
# table(res1$s2_resid > 4, res1$model)
# 
# table(res1$s2_sp_bias > 4, res1$model)
# table(res1$s2_phylo_bias > 4, res1$model)
# table(res1$s2_resid_bias > 4, res1$model)
# 
# 
# hist(res1$s2_phylo_bias[which(res1$model=="INLA")], breaks=100)
# hist(res1$s2_phylo_bias[which(res1$model=="INLA" & res1$s2_phylo_bias < 100)], breaks=100)



################## 2. Formatting ############################

#------ labels for plots

res1$species_lab <- factor(res1$species_size, 
                           labels=c(
                             `25`=parse(text=TeX("$species=25$")),
                             `50`=parse(text=TeX("$species=50$")),
                             `100`=parse(text=TeX("$species=100$"))
                           ))

res2$species_lab <- factor(res2$species_size, 
                           labels=c(
                             `200`=parse(text=TeX("$species=200$")),
                             `400`=parse(text=TeX("$species=400$")),
                             `800`=parse(text=TeX("$species=800$"))
                           ))


#----------- Reorder factors for figures
res1$model <- factor(res1$model, 
                     levels=c("brms", "MCMCglmm", "INLA", "glmmTMB",  "phyr"))

res2$model <- factor(res2$model, 
                     levels=c("MCMCglmm", "INLA", "glmmTMB"))



# create new label nreps based on sample_size
res1 <- res1 |> 
  mutate(nrep = case_when(
    sample_size == 125 & species_size == 25 ~ "5",
    sample_size == 250 & species_size == 25 ~ "10",
    sample_size == 750 & species_size == 25 ~ "30",
    sample_size == 250 & species_size == 50 ~ "5",
    sample_size == 500 & species_size == 50 ~ "10",
    sample_size == 1500 & species_size == 50 ~ "30",
    sample_size == 500 & species_size == 100 ~ "5",
    sample_size == 1000 & species_size == 100 ~ "10",
    sample_size == 3000 & species_size == 100 ~ "30",
    .default = "unbalanced"
  ))

#table(res1$nrep)

res2 <- res2 |> 
  mutate(nrep = case_when(
    sample_size == 1000 & species_size == 200 ~ "5",
    sample_size == 2000 & species_size == 200 ~ "10",
    sample_size == 6000 & species_size == 200 ~ "30",
    sample_size == 2000 & species_size == 400 ~ "5",
    sample_size == 4000 & species_size == 400 ~ "10",
    sample_size == 12000 & species_size == 400 ~ "30",
    sample_size == 4000 & species_size == 800 ~ "5",
    sample_size == 8000 & species_size == 800 ~ "10",
    sample_size == 24000 & species_size == 800 ~ "30",
    .default = "unbalanced"
  ))
#table(res2$nrep)


# re-order factor levels of nrep
res1$nrep <- factor(res1$nrep, 
                    levels=c("5", "10", "30", "unbalanced"),
                    labels=c("5 reps", "10 reps", "30 reps", "unbalanced"))

res2$nrep <- factor(res2$nrep,
                    levels=c("5", "10", "30", "unbalanced"),
                    labels=c("5 reps", "10 reps", "30 reps", "unbalanced"))








################### 3. Derive sim estimates #############################

### all sim studies had true b1=1.5

# bias datasets for b1 fixed effect coeff.
bias.dat.b1 <- res1 |> 
  group_by(model, species_size, species_lab, nrep, sigma2.s, sigma2.p, sigma2.e) |> 
  summarise(bias = mean(mu) - 1.5) |> 
  ungroup()

bias.dat2.b1 <- res2 |> 
  group_by(model, species_size, species_lab, nrep, sigma2.s, sigma2.p, sigma2.e) |> 
  summarise(bias = mean(mu) - 1.5) |> 
  ungroup()


# RMSE datasets for b1 fixed effect coeff.
rmse.dat.b1 <- res1 |> 
  group_by(model, species_size, species_lab, nrep, sigma2.s, sigma2.p, sigma2.e) |> 
  summarise(rmse = sqrt(mean((mu - b1)^2))) |> 
  ungroup()

rmse.dat2.b1 <- res2 |> 
  group_by(model, species_size, species_lab, nrep, sigma2.s, sigma2.p, sigma2.e) |> 
  summarise(rmse = sqrt(mean((mu - b1)^2))) |> 
  ungroup()


# CI width for b1 fixed effect coeff.
ciw.dat.b1 <- res1 |> 
  group_by(model, species_size, species_lab, nrep, sigma2.s, sigma2.p, sigma2.e) |> 
  summarise(ci_width = mean(mu_ci_width)) |> 
  ungroup()

ciw.dat2.b1 <- res2 |> 
  group_by(model, species_size, species_lab, nrep, sigma2.s, sigma2.p, sigma2.e) |> 
  summarise(ci_width = mean(mu_ci_width)) |> 
  ungroup()


# RMSE sigma2.p estimate
rmse.dat.s2p <- res1 |> 
  group_by(model, species_size, species_lab, nrep, sigma2.s, sigma2.p, sigma2.e) |> 
  summarise(rmse = sqrt(mean((s2_phylo - sigma2.p)^2))) |> 
  ungroup()

rmse.dat2.s2p <- res2 |> 
  group_by(model, species_size, species_lab, nrep, sigma2.s, sigma2.p, sigma2.e) |> 
  summarise(rmse = sqrt(mean((s2_phylo - sigma2.p)^2))) |> 
  ungroup()



# RMSE sigma2.s estimate
rmse.dat.s2sp <- res1 |> 
  group_by(model, species_size, species_lab, nrep, sigma2.s, sigma2.p, sigma2.e) |> 
  summarise(rmse = sqrt(mean((s2_sp - sigma2.s)^2))) |> 
  ungroup()

rmse.dat2.s2sp <- res2 |> 
  group_by(model, species_size, species_lab, nrep, sigma2.s, sigma2.p, sigma2.e) |> 
  summarise(rmse = sqrt(mean((s2_sp - sigma2.s)^2))) |> 
  ungroup()




# RMSE sigma2.e estimate
rmse.dat.s2e <- res1 |> 
  group_by(model, species_size, species_lab, nrep, sigma2.s, sigma2.p, sigma2.e) |> 
  summarise(rmse = sqrt(mean((s2_resid - sigma2.e)^2))) |> 
  ungroup()

rmse.dat2.s2e <- res2 |> 
  group_by(model, species_size, species_lab, nrep, sigma2.s, sigma2.p, sigma2.e) |> 
  summarise(rmse = sqrt(mean((s2_resid - sigma2.e)^2))) |> 
  ungroup()


# coverage b1 estimate
res1$cov_b1 <- res1$b1 > res1$mu_ci_low & res1$b1 < res1$mu_ci_high
res2$cov_b1 <- res2$b1 > res2$mu_ci_low & res2$b1 < res2$mu_ci_high

# coverage sigma2.p estimate
res1$cov_sigma2.p <- res1$sigma2.p > res1$s2_phylo_ci_low & res1$sigma2.p < res1$s2_phylo_ci_high
res2$cov_sigma2.p <- res2$sigma2.p > res2$s2_phylo_ci_low & res2$sigma2.p < res2$s2_phylo_ci_high

# coverage sigma2.s estimate
res1$cov_sigma2.s <- res1$sigma2.s > res1$s2_sp_ci_low & res1$sigma2.s < res1$s2_sp_ci_high
res2$cov_sigma2.s <- res2$sigma2.s > res2$s2_sp_ci_low & res2$sigma2.s < res2$s2_sp_ci_high

# coverage sigma2.e estimate
res1$cov_sigma2.e <- res1$sigma2.e > res1$s2_resid_ci_low & res1$sigma2.e < res1$s2_resid_ci_high
res2$cov_sigma2.e <- res2$sigma2.e > res2$s2_resid_ci_low & res2$sigma2.e < res2$s2_resid_ci_high


# derive the coverage dataset for b1 fixed effect coeff.
cov.dat.b1 <- res1 |> 
  group_by(model, species_size, species_lab, nrep, sigma2.s, sigma2.p, sigma2.e) |> 
  summarise(cov_prop = mean(cov_b1, na.rm = TRUE)) |> 
  ungroup()

cov.dat2.b1 <- res2 |> 
  group_by(model, species_size, species_lab, nrep, sigma2.s, sigma2.p, sigma2.e) |> 
  summarise(cov_prop = mean(cov_b1, na.rm = TRUE)) |> 
  ungroup()



## derive sample variance of estimates (for MCSE)
res1.sample_var <- res1 |> 
  group_by(model, species_size) |>
  summarise(bias_b1 = mean(mu) - 1.5,
            rmse_b1 = sqrt(mean((mu - 1.5)^2)),
            b1_S2 = sum((mu - mean(mu))^2) / (n() - 1),
            s2.s_rmse = sqrt(mean((s2_sp - sigma2.s)^2)), 
            s2.p_rmse = sqrt(mean((s2_phylo - sigma2.p)^2)), 
            s2.e_rmse = sqrt(mean((s2_resid - sigma2.e)^2)), 
            s2.s_S2 = sum((s2_sp - mean(s2_sp))^2) / (n() - 1),
            s2.p_S2 = sum((s2_phylo - mean(s2_phylo))^2) / (n() - 1),
            s2.e_S2 = sum((s2_resid - mean(s2_resid))^2) / (n() - 1)) |> 
  ungroup()


res2.sample_var <- res2 |> 
  group_by(model, species_size) |>
  summarise(bias_b1 = mean(mu) - 1.5,
            rmse_b1 = sqrt(mean((mu - 1.5)^2)),
            b1_S2 = sum((mu - mean(mu))^2) / (n() - 1),
            s2.s_rmse = sqrt(mean((s2_sp - sigma2.s)^2)), 
            s2.p_rmse = sqrt(mean((s2_phylo - sigma2.p)^2)), 
            s2.e_rmse = sqrt(mean((s2_resid - sigma2.e)^2)), 
            s2.s_S2 = sum((s2_sp - mean(s2_sp))^2) / (n() - 1),
            s2.p_S2 = sum((s2_phylo - mean(s2_phylo))^2) / (n() - 1),
            s2.e_S2 = sum((s2_resid - mean(s2_resid))^2) / (n() - 1)) |> 
  ungroup()







################### 4. Plots/tables - runtime #############################



####### Fig. 1: run time of all models per species size

options(digits = 1, scipen = 5)

# Plot run times vs packages (y-axis on log10 scale) - set1
runtime_plot_set1 <- res1 |> 
  #filter(model %in% c("glmmTMB", "INLA", "brms", "MCMCglmm")) |> 
  ggplot(aes(x=factor(model), y=run_time, fill=model, color=model)) + 
  geom_violin()+
  scale_y_log10(breaks=c(0, 0.1, 1, 10, 100, 1000, 10000), 
                labels=c("0", "0.1", "1", "10", "100", "1000", "10000")) +
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title="",x="Package", y = "Run time (seconds)")+
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")

#ggsave(filename = "output/Figure_log10runtime_set1.png", width = 10, height = 6)



# Plot run times vs packages (y-axis on log10 scale) - set2
runtime_plot_set2 <- res2 |> 
  filter(model %in% c("MCMCglmm", "INLA", "glmmTMB")) |> 
  ggplot(aes(x=factor(model), y=run_time,
             fill=factor(model), color=factor(model))) + 
  geom_violin()+
  scale_y_log10(breaks=c(0, 0.1, 1, 10, 100, 1000, 10000, 100000), 
                labels=c("0", "0.1", "1", "10", "100", "1000", "10000", "100000")) +
  scale_color_manual(values=c("#7AB47C", "#FBBF24", "#B47AA5")) +
  scale_fill_manual(values=alpha(c("#7AB47C", "#FBBF24", "#B47AA5"), 0.4)) + 
  labs(title="",x="Package", y = "Run time (seconds)")+
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  #theme(legend.position="bottom")+
  theme(legend.position="none")

#ggsave(filename = "output/Figure_log10runtime_set1.png", width = 10, height = 6)



# Combine plots
fig_runtime <- runtime_plot_set1 + 
  runtime_plot_set2 + 
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(tag_levels='A',
                  theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(filename = "output/Figures/Figure_runtime_all.png", width = 10, height = 8)





###----- For supp. information ---

##### Tables 
# Get mean values per species size and package
mean_runtime <- rbind(res1, res2) |> 
  group_by(species_size, model) |> 
  summarise(mean_run_time = mean(run_time, na.rm = TRUE)) |> 
  ungroup()

# Create table of mean run times
tab_runtime <- mean_runtime |> 
  pivot_wider(names_from = species_size, values_from = mean_run_time) |> 
  mutate(model = factor(model, levels=c("brms", "MCMCglmm", "INLA", "glmmTMB", "phyr"))) |> 
  arrange(model)

xtable(tab_runtime, digits=2, 
       caption="Mean run times of models per species size (in seconds).")



###### Plot Log-log runtime vs species number 
# To see which packages take exponential time to compute with larger species number
runtime_log_plot <- bind_rows(res1, res2) |>
  mutate(model = factor(trimws(model))) |>
  ggplot(aes(x = species_size, y = run_time,
             fill = model, colour = model)) +
  geom_boxplot(aes(group = interaction(model, species_size)),
               width = 0.1, position = "identity", outlier.shape = NA) +
  scale_x_log10(breaks = c(25, 50, 100, 200, 400, 800),
                labels = c("25","50","100","200","400","800")) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000),
                labels = c("0.1","1","10","100","1000","10000","100000")) +
  geom_smooth(
    aes(group = model, colour = model),
    method = "lm", formula = y ~ x,
    se = FALSE, linetype = "dashed", linewidth = 0.9
  ) +
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"))) +
  labs(x = "Species number", y = "Run time (seconds)") +
  theme_bw()  
  # # linear (grey dashed, shifted down a bit)
  # geom_smooth(
  #   mapping = aes(x = species_size, y = run_time/170),
  #   method = "lm", formula = y ~ x,
  #   se = FALSE, colour = "grey40", linetype = "dashed", linewidth = 0.9,
  #   inherit.aes = FALSE
  # ) +
  # # quadratic (grey solid, shifted further down)
  # geom_smooth(
  #   mapping = aes(x = species_size, y = run_time*), 
  #   method = "lm", formula = y ~ poly(x, 2),
  #   se = FALSE, colour = "grey40", linewidth = 0.9,
  #   inherit.aes = FALSE
  # )
#+ theme(legend.position = "none")

ggsave(filename = "output/Figure_log_runtime_perspecies.png", width = 10, height = 6)


log10_slope_summary <- bind_rows(res1, res2) |>
  mutate(model = factor(trimws(model))) |>
  group_by(model) |>
  do(broom::tidy(lm(log10(run_time) ~ log10(species_size), data = .), conf.int = TRUE)) |>
  filter(term == "log10(species_size)") |>
  select(model, slope = estimate, conf.low, conf.high, p.value)

log10_slope_summary




################# 5. Plots/tables - fixed effects ##########################



options(digits = 4, scipen = 5) # switch back to display multiple digits

# First: use N_species as facet --


# distribution of b1 ---
b1_plot <- res1 |> 
  ggplot(aes(x=factor(model), y=mu, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("$\\hat{\\beta_1}$ distribution"),x="Package", y = TeX("$\\hat{\\beta_1}$"))+
  #scale_y_continuous(breaks=seq(0, 3, 0.5), limits=c(0, 3)) +
  geom_hline(yintercept=1.5, colour="darkgray", linewidth=0.5)+
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")


# bias of b1 ---
b1_bias_plot <- bias.dat.b1 |> 
  ggplot(aes(x=factor(model), y=bias, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("Bias $\\hat{\\beta_1}$"),x="Package", y = TeX("Bias $\\hat{\\beta_1}$"))+
  #scale_y_continuous(breaks=seq(0, 3, 0.5), limits=c(0, 3)) +
  geom_hline(yintercept=0, colour="darkgray", linewidth=0.5)+
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")


# RMSE ---
b1_rmse_plot <- rmse.dat.b1 |> 
  ggplot(aes(x=factor(model), y=rmse, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("RMSE $\\hat{\\beta_1}$"),x="Package", y = TeX("RMSE $\\hat{\\beta_1}$"))+
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  #facet_wrap(~nrep, ncol=4)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")



# scaled RMSE ---
b1_rmse_scaled_plot <- res1 |> 
  ggplot(aes(x=factor(model), y=rmse_b1_scaled, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("Scaled RMSE $\\hat{\\beta_1}$"),x="Package", y = TeX("Scaled RMSE $\\hat{\\beta_1}$"))+
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")


# CI width ---
b1_ci_width_plot <- ciw.dat.b1 |> 
  ggplot(aes(x=factor(model), y=ci_width, fill=model, color=model)) + 
  geom_boxplot(width=0.4)+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("95% confidence interval width $\\hat{\\beta_1}$"),x="Package", y = TeX("$\\hat{\\beta_1}$ CI width"))+
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+ ## alternatively I can use nrep here --- discuss!
  #facet_wrap(~nrep, ncol=4)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")


# coverage ---
b1_cov_plot <- cov.dat.b1 |> 
  ggplot(aes(x=factor(model), y=cov_prop, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  scale_y_continuous(breaks=seq(0.9, 1, 0.01), limits=c(0.9, 1)) +
  labs(title=TeX("Coverage $\\hat{\\beta_1}$"),x="Package", y = TeX("Coverage $\\hat{\\beta_1}$"))+
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  #facet_wrap(~nrep, ncol=4)+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  theme_bw() +
  theme(legend.position="none")



# combine plots
fig_b1 <- b1_rmse_plot + 
  b1_ci_width_plot + 
  b1_cov_plot + 
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(tag_levels='A',
                  theme = theme(plot.background = element_rect(fill = "white", colour = NA)))


ggsave(filename = "output/Figure_b1_estimate_Nspecies.png", width = 9, height = 11)




# Now use n_rep as the facet --


# distribution of b1 ---
b1_plot <- res1 |> 
  ggplot(aes(x=factor(model), y=mu, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("$\\hat{\\beta_1}$ distribution"),x="Package", y = TeX("$\\hat{\\beta_1}$"))+
  #scale_y_continuous(breaks=seq(0, 3, 0.5), limits=c(0, 3)) +
  geom_hline(yintercept=1.5, colour="darkgray", linewidth=0.5)+
  facet_wrap(~nrep, ncol=4)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")


# bias of b1 ---
b1_bias_plot <- bias.dat.b1 |> 
  ggplot(aes(x=factor(model), y=bias, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("Bias $\\hat{\\beta_1}$"),x="Package", y = TeX("Bias $\\hat{\\beta_1}$"))+
  #scale_y_continuous(breaks=seq(0, 3, 0.5), limits=c(0, 3)) +
  geom_hline(yintercept=0, colour="darkgray", linewidth=0.5)+
  facet_wrap(~nrep, ncol=4)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")


# RMSE ---
b1_rmse_plot <- rmse.dat.b1 |> 
  ggplot(aes(x=factor(model), y=rmse, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("RMSE $\\hat{\\beta_1}$"),x="Package", y = TeX("RMSE $\\hat{\\beta_1}$"))+
  facet_wrap(~nrep, ncol=4)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")


# CI width ---
b1_ci_width_plot <- ciw.dat.b1 |> 
  ggplot(aes(x=factor(model), y=ci_width, fill=model, color=model)) + 
  geom_boxplot(width=0.4)+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("95% confidence interval width $\\hat{\\beta_1}$"),x="Package", y = TeX("$\\hat{\\beta_1}$ CI width"))+
  facet_wrap(~nrep, ncol=4)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")


# coverage ---
b1_cov_plot <- cov.dat.b1 |> 
  ggplot(aes(x=factor(model), y=cov_prop, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  scale_y_continuous(breaks=seq(0.9, 1, 0.01), limits=c(0.9, 1)) +
  labs(title=TeX("Coverage $\\hat{\\beta_1}$"),x="Package", y = TeX("Coverage $\\hat{\\beta_1}$"))+
  facet_wrap(~nrep, ncol=4)+
  geom_boxplot(width=0.4)+
  geom_hline(yintercept=0.95, colour="darkgray", linewidth=0.6)+ 
  theme_bw() +
  theme(legend.position="none")


# combine plots
fig_b1 <- b1_rmse_plot + 
  b1_ci_width_plot + 
  b1_cov_plot + 
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(tag_levels='A',
                  theme = theme(plot.background = element_rect(fill = "white", colour = NA)))


ggsave(filename = "output/Figure_b1_estimate_nreps.png", width = 11, height = 11)






################# 6. Plots/tables - variance components ##########################


### 1. Non-phylo variance estimate ---

options(digits = 3, scipen = 5)

# RMSE of variance estimate 
rmse.s2sp <- rbind(rmse.dat.s2sp, rmse.dat2.s2sp)

s2.sp_plot1_rmse <- rmse.s2sp |> 
  ggplot(aes(x=factor(model), y=rmse, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX(""),x="Package", y = TeX("RMSE $\\hat{\\sigma^2_s}"))+
  #scale_y_continuous(breaks=seq(-0.5, 3, 0.5), limits=c(-0.5, 3)) +
  #facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")




# for 0.25 settings
s2.sp_plot1 <- res1 |> 
  filter(sigma2.s==0.25, sigma2.p==0.25, nrep=="10 reps") |>
  ggplot(aes(x=factor(model), y=s2_sp, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(x="Package", y = TeX("$\\hat{\\sigma^2_s}"))+
  scale_y_continuous(breaks=seq(0, 4, 1), limits=c(0, 4)) +
  geom_hline(yintercept=0.25, colour="darkgray", linewidth=0.5)+
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")


# for 0.05 settings
s2.sp_plot2 <- res1 |> 
  filter(sigma2.s==0.05, sigma2.p==0.05, nrep=="10 reps") |>
  ggplot(aes(x=factor(model), y=s2_sp, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("$\\sigma^2_s=0.05$, $\\sigma^2_p=0.05$"),x="Package", y = TeX("$\\hat{\\sigma^2_s}"))+
  scale_y_continuous(breaks=seq(0, 0.9, 0.05), limits=c(0, 0.9)) +
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  #facet_wrap(~nrep, ncol=4)+
  geom_hline(yintercept=0.05, colour="darkgray", linewidth=0.5)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")





### 2. Phylo variance estimate --

options(digits = 3, scipen = 5)

# RMSE of variance estimate 
rmse.s2p <- rbind(rmse.dat.s2p, rmse.dat2.s2p)

s2.p_plot_rmse <- rmse.s2p |> 
  ggplot(aes(x=factor(model), y=rmse, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX(""),x="Package", y = TeX("RMSE $\\hat{\\sigma^2_p}"))+
  #scale_y_continuous(breaks=seq(-0.5, 3, 0.5), limits=c(-0.5, 3)) +
  #facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")


# for 0.25 settings
s2.p_plot1 <- res1 |> 
  filter(sigma2.s==0.25, sigma2.p==0.25, nrep=="10 reps") |>
  ggplot(aes(x=factor(model), y=s2_phylo, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX(""),x="Package", y = TeX("$\\hat{\\sigma^2_p}"))+
  scale_y_continuous(breaks=seq(0, 4, 1), limits=c(0, 4)) +
  geom_hline(yintercept=0.25, colour="darkgray", linewidth=0.5)+
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")


# for 0.05 settings
s2.p_plot2 <- res1 |> 
  filter(sigma2.s==0.05, sigma2.p==0.05, nrep=="10 reps") |>
  ggplot(aes(x=factor(model), y=s2_phylo, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("$\\sigma^2_s=0.05$, $\\sigma^2_p=0.05$"),x="Package", y = TeX("$\\hat{\\sigma^2_p}"))+
  #scale_y_continuous(breaks=seq(-0.5, 2, 0.5), limits=c(-0.5, 2)) +
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  #facet_wrap(~nrep, ncol=4)+
  geom_hline(yintercept=0.05, colour="darkgray", linewidth=0.5)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")






### 3. Residual variance estimate ---


options(digits = 3, scipen = 5)

# RMSE of variance estimate 
rmse.s2e <- rbind(rmse.dat.s2e, rmse.dat2.s2e)

s2.e_plot_rmse <- rmse.s2e |> 
  ggplot(aes(x=factor(model), y=rmse, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX(""),x="Package", y = TeX("RMSE $\\hat{\\sigma^2_e}"))+
  #scale_y_continuous(breaks=seq(-0.5, 3, 0.5), limits=c(-0.5, 3)) +
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")



# for 0.25 settings
s2.e_plot1 <- res1 |> 
  filter(sigma2.s==0.25, sigma2.p==0.25, nrep=="10 reps") |>
  ggplot(aes(x=factor(model), y=s2_resid, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX(""),x="Package", y = TeX("$\\hat{\\sigma^2_e}"))+
  scale_y_continuous(breaks=seq(0, 0.4, 0.05), limits=c(0.1, 0.3)) +
  geom_hline(yintercept=0.2, colour="darkgray", linewidth=0.5)+
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")


# for 0.05 settings
s2.e_plot2 <- res1 |> 
  filter(sigma2.s==0.05, sigma2.p==0.05, nrep=="10 reps") |>
  ggplot(aes(x=factor(model), y=s2_resid, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("$\\sigma^2_s=0.05$, $\\sigma^2_p=0.05$"),x="Package", y = TeX("$\\hat{\\sigma^2_e}"))+
  scale_y_continuous(breaks=seq(0.1, 0.4, 0.05), limits=c(0.1, 0.4)) +
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  #facet_wrap(~nrep, ncol=4)+
  geom_hline(yintercept=0.2, colour="darkgray", linewidth=0.5)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")




# create labels for plots
#TeX("$\\sigma^2_s=0.25$, $\\sigma^2_p=0.25$, $n_{rep}=10$")


# Combine plots of variance estimates
fig_variance <- s2.sp_plot1 / s2.p_plot1 / s2.e_plot1 +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = TeX("$\\sigma^2_s=0.25, \\sigma^2_p=0.25, n_{rep}=10$"),
    tag_levels = 'A',
    theme = theme(
      plot.title = element_text(size = 14, color="black", face = "bold", hjust = 0.5),
      plot.background = element_rect(fill = "white", colour = NA))
  )



ggsave(filename = "output/Figure_variance_estimates_suppinfo.png", width = 9, height = 10)




# Combine plots of variance estimates (Supp. Info)
fig_variance <- s2.sp_plot1_rmse / s2.p_plot_rmse / s2.e_plot1_rmse +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(plot.background = element_rect(fill = "white", colour = NA))
  )

ggsave(filename = "output/Figure_variance_estimates.png", width = 4, height = 7)








######### 7. Derive Monte Carlo Standard Errors (MCSE) of simulations ############


## 
sample_var <- bind_rows(res1.sample_var, res2.sample_var)


## derive the bias Monte Carlo SE (per model, method and condition) for overall mean
b1_mcse <- sample_var |> 
  group_by(model, species_size) |> 
  summarise(#rmse_b1 = rmse_b1,
            mcse_b1 = round(sqrt(b1_S2/n()),5)) |> 
  arrange(species_size) 
print(xtable(b1_mcse, digits=c(0,2,2,5,5)), include.rownames=FALSE) ##save table for supporting information



## derive the coverage Monte Carlo SE (per model, method and condition)
cov.dat <- bind_rows(cov.dat.b1, cov.dat2.b1)
cov_mcse <- cov.dat |> 
  group_by(model, species_size) |> 
  summarise(mean_cov = mean(cov_prop),  # Compute the mean coverage
            cov_mcse = round(sqrt((mean_cov * (1 - mean_cov)) / n()), 5)) |>  # Apply MCSE formula
  arrange(species_size) 
print(xtable(cov_mcse, digits=c(0,2,2,3,4)), include.rownames=FALSE) ##save table for supporting information


## derive the bias Monte Carlo SE (per model, method and condition) for variance components
s2_mcse <- sample_var |> 
  group_by(model, species_size) |> 
  summarise(s2.p_mcse = round(sqrt(s2.p_S2/n()),5),
            s2.e_mcse = round(sqrt(s2.e_S2/n()),5)) |> 
  arrange(species_size) 
print(xtable(s2_mcse, digits=c(0,3,3,3,3)), include.rownames=FALSE) ##save table for supporting information

