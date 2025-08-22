################### Load results #############################

library(ggplot2);library(devtools);library(cowplot);library(ggdark);theme_set(dark_theme_bw())
library(ggdist);library(tidyverse);library(dplyr);library(patchwork); library(scales);
library(latex2exp); library(xtable)

#---
# load results set 1a (repeated measures)
load("~/PhD/3_PhyloTMB/sims/sim_set1a_results_pglmm.RDATA")
res1.dat <- dat

# load results set 2a (repeated measures)
load("~/PhD/3_PhyloTMB/sims/sim_set2a_results_pglmm.RDATA")
res2.dat <- dat

# #---
# # load results set 1b (one measure per species)
# load("~/PhD/3_PhyloTMB/sims/sim_results_set1b.RDATA")
# res1b <- dat
# 
# # load results set 2b (one measure per species)
# load("~/PhD/3_PhyloTMB/sims/sim_results_set2b.RDATA")
# res2b <- dat




################# 1. Convergence #########################


# Check convergence %
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

# # merge convergence info
# conv.tab <- full_join(res1.dat_conv_summary, res2.dat_conv_summary, by = "model") %>%
#   rename(Package = model) %>%
#   select(Package, `25`, `50`, `100`, `200`, `400`, `800`)


# # Create summary tables of convergence % for Supp. information
# print(xtable(conv.tab, digits = 1),
#       include.rownames = FALSE)
# check ESS
#table(res1.dat$ESS>400, res1.dat$model)
#table(res2.dat$ESS>400, res2.dat$model)




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

table(res1$species_size)/5/64000
table(res2$species_size)/3/8000

table(res1$species_size)
table(res2$species_size)

# assess model errors
table(is.na(res1$model_error), res1$model)
# assess model warnings
table(is.na(res1.dat$model_warning), res1.dat$model)
table(is.na(res2.dat$model_warning), res2.dat$model)


# all INLA models have the same warning
table(res1$model_warning[which(res1$model=="INLA")])
table(res2$model_warning[which(res1$model=="INLA")])


table(res1$s2_sp > 4, res1$model)
table(res1$s2_phylo > 4, res1$model)
table(res1$s2_resid > 4, res1$model)

table(res1$s2_sp_bias > 4, res1$model)
table(res1$s2_phylo_bias > 4, res1$model)
table(res1$s2_resid_bias > 4, res1$model)


hist(res1$s2_phylo_bias[which(res1$model=="INLA")], breaks=100)
hist(res1$s2_phylo_bias[which(res1$model=="INLA" & res1$s2_phylo_bias < 100)], breaks=100)







################### 2. Derive variables #############################

# bias of sigma2.p estimate
res1$s2_phylo_bias <- res1$s2_phylo - res1$sigma2.p
res2$s2_phylo_bias <- res2$s2_phylo - res2$sigma2.p

# bias of sigma2.s estimate
res1$s2_sp_bias <- res1$s2_sp - res1$sigma2.s
res2$s2_sp_bias <- res2$s2_sp - res2$sigma2.s

# bias of sigma2.e estimate
res1$s2_resid_bias <- res1$s2_resid - res1$sigma2.e
res2$s2_resid_bias <- res2$s2_resid - res2$sigma2.e



# RMSE b1 estimate
res1$rmse_b1 <- sqrt((res1$mu - res1$b1)^2)
res2$rmse_b1 <- sqrt((res2$mu - res2$b1)^2)

# RMSE sigma2.p estimate
res1$rmse_s2_phylo <- sqrt((res1$sigma2.p - res1$s2_phylo)^2)
res2$rmse_s2_phylo <- sqrt((res2$sigma2.p - res2$s2_phylo)^2)

# RMSE sigma2.s estimate
res1$rmse_s2_sp <- sqrt((res1$sigma2.s - res1$s2_sp)^2)
res2$rmse_s2_sp <- sqrt((res2$sigma2.s - res2$s2_sp)^2)

# RMSE sigma2.e estimate
res1$rmse_s2_resid <- sqrt((res1$sigma2.e - res1$s2_resid)^2)
res2$rmse_s2_resid <- sqrt((res2$sigma2.e - res2$s2_resid)^2)


# scaled RMSE b1 estimate
res1$rmse_b1_scaled <- sqrt((res1$mu/res1$b1 - 1)^2)
res2$rmse_b1_scaled <- sqrt((res2$mu/res2$b1 - 1)^2)

# scaled RMSE sigma2.p estimate
res1$rmse_s2_phylo_scaled <- sqrt((res1$sigma2.p/res1$s2_phylo - 1)^2)
res2$rmse_s2_phylo_scaled <- sqrt((res2$sigma2.p/res2$s2_phylo - 1)^2)

# scaled RMSE sigma2.s estimate
res1$rmse_s2_sp_scaled <- sqrt((res1$sigma2.s/res1$s2_sp - 1)^2)
res2$rmse_s2_sp_scaled <- sqrt((res2$sigma2.s/res2$s2_sp - 1)^2)

# scaled RMSE sigma2.e estimate
res1$rmse_s2_resid_scaled <- sqrt((res1$sigma2.e/res1$s2_resid - 1)^2)
res2$rmse_s2_resid_scaled <- sqrt((res2$sigma2.e/res2$s2_resid - 1)^2)



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


# create new label nreps based on sample_size
res1 <- res1 |> 
  mutate(nrep = case_when(
    sample_size == 125 & species_size == 25 ~ "5",
    sample_size == 250 & species_size == 25 ~ "10",
    sample_size == 100 & species_size == 25 ~ "30",
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
    sample_size == 1200 & species_size == 400 ~ "30",
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





################## 3. Formatting ############################

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






################### 4. Plots for main text #############################



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

ggsave(filename = "output/Figure_log10runtime_set1.png", width = 10, height = 6)



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

ggsave(filename = "output/Figure_log10runtime_set1.png", width = 10, height = 6)



# Combine plots
fig_runtime <- runtime_plot_set1 + 
  runtime_plot_set2 + 
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(tag_levels='A',
                  theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(filename = "output/Figures/Figure_runtime_all.png", width = 10, height = 8)





#### Fig. 2: variance components







################# 5. Plots/tables for supp. info ##########################


### 1. Run time -------------------------------------

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




### 2. Fixed effect estimate -----------------------------------

# derive the coverage dataset
cov.dat.b1 <- res1 |> 
  group_by(model, species_size, species_lab, nrep, sigma2.s, sigma2.p, sigma2.e) |> 
  summarise(cov_prop = mean(cov_b1, na.rm = TRUE)) |> 
  ungroup()


options(digits = 4, scipen = 5) # switch back to display multiple digits

# bias of b1 ---
b1_bias_plot <- res1 |> 
  ggplot(aes(x=factor(model), y=mu, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("Bias $\\hat{\\beta_1}$"),x="Package", y = TeX("Bias $\\hat{\\beta_1}$"))+
  #scale_y_continuous(breaks=seq(0, 3, 0.5), limits=c(0, 3)) +
  geom_hline(yintercept=1.5, colour="darkgray", linewidth=0.5)+
  #facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")



# RMSE ---
b1_rmse_plot <- res1 |> 
  ggplot(aes(x=factor(model), y=rmse_b1, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("RMSE $\\hat{\\beta_1}$"),x="Package", y = TeX("RMSE $\\hat{\\beta_1}$"))+
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  facet_wrap(~nrep, ncol=4)+
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
b1_ci_width_plot <- res1 |> 
  ggplot(aes(x=factor(model), y=mu_ci_width, fill=model, color=model)) + 
  geom_boxplot(width=0.4)+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("95% confidence interval width $\\hat{\\beta_1}$"),x="Package", y = TeX("$\\hat{\\beta_1}$ CI width"))+
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+ ## alternatively I can use nrep here --- discuss!
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
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
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


ggsave(filename = "output/Figure_b1_estimate_nreps.png", width = 11.5, height = 13)



##


# RMSE ---
b1_rmse_plot <- res1 |> 
  ggplot(aes(x=factor(model), y=rmse_b1, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("RMSE $\\hat{\\beta_1}$"),x="Package", y = TeX("RMSE $\\hat{\\beta_1}$"))+
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")


# CI width ---
b1_ci_width_plot <- res1 |> 
  ggplot(aes(x=factor(model), y=mu_ci_width, fill=model, color=model)) + 
  geom_boxplot(width=0.4)+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("95% confidence interval width $\\hat{\\beta_1}$"),x="Package", y = TeX("$\\hat{\\beta_1}$ CI width"))+
  facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+ ## alternatively I can use nrep here --- discuss!
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


ggsave(filename = "output/Figure_b1_estimate_nspecies.png", width = 11.5, height = 13)






### 3. Non-phylo variance estimate -----------------------------------

options(digits = 1, scipen = 5)
# log10 raw variance estimate values 
s2.sp_plot <- res1 |> 
  filter(sigma2.s==0.25) |>
  #filter(model != "INLA") |>
  ggplot(aes(x=factor(model), y=s2_sp, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("Log10 $\\sigma^2_s$ estimate"),x="Package", y = TeX("Log10 $\\sigma^2_s$ estimate"))+
  scale_y_log10() +
  geom_hline(yintercept=0.05, colour="darkgray", linewidth=0.5)+
  #facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")

ggsave(filename = "output/Fig.test_sigma2.sp_log10.png", width = 10, height = 7)


options(digits = 3, scipen = 5)
# bias of variance estimate
s2.sp_plot_bias <- res1 |> 
  #filter(model != "INLA") |>
  ggplot(aes(x=factor(model), y=s2_sp_bias, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("Bias $\\hat{\\sigma^2_s}$"),x="Package", y = TeX("Bias $\\hat{\\sigma^2_s}"))+
  scale_y_continuous(breaks=seq(-0.5, 3, 0.5), limits=c(-0.5, 3)) +
  geom_hline(yintercept=0, colour="darkgray", linewidth=0.5)+
  #facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")

ggsave(filename = "output/Figure_sigma2_sp_bias.png", width = 10, height = 7)


# RMSE
s2.p_plot_rmse <- res1 |> 
  #filter(model != "INLA") |>
  ggplot(aes(x=factor(model), y=rmse_s2_sp, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title="Plot of RMSE in non-phylo variance estimates",x="Package", y = "RMSE")+
  scale_y_continuous(breaks=seq(-0.5, 2, 0.5), limits=c(-0.5, 2)) +
  #facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  #facet_wrap(~nrep, ncol=4)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")

ggsave(filename = "output/Figure_sigma2_sp_rmse.png", width = 10, height = 7)








### 4. Phylo variance estimate -----------------------------------

options(digits = 1, scipen = 5)
# log10 values 
s2.p_plot <- res1 |> 
  #filter(model != "INLA") |>
  ggplot(aes(x=factor(model), y=s2_phylo, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("Log10 $\\sigma^2_p$ estimate"),x="Package", y = TeX("Log10 $\\sigma^2_p$ estimate"))+
  scale_y_log10() +
  geom_hline(yintercept=0, colour="darkgray", linewidth=0.5)+
  #facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")

ggsave(filename = "output/Fig.test_sigma2.phylo_log10.png", width = 10, height = 7)


options(digits = 3, scipen = 5)
# bias 
s2.p_plot_bias <- res1 |> 
  #filter(model != "INLA") |>
  ggplot(aes(x=factor(model), y=s2_phylo_bias, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("Bias $\\sigma^2_p$ estimate"),x="Package", y = TeX("Bias $\\sigma^2_p$ estimate"))+
  scale_y_continuous(breaks=seq(-0.5, 5, 0.5), limits=c(-0.5, 5)) + 
  geom_hline(yintercept=0, colour="darkgray", linewidth=0.5)+
  #facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")

ggsave(filename = "output/Figure_sigma2_phylo_bias.png", width = 10, height = 7)


# RMSE
s2.p_plot_rmse <- res1 |> 
  #filter(model != "INLA") |>
  ggplot(aes(x=factor(model), y=rmse_s2_phylo, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title="Plot of RMSE in phylo variance estimates",x="Package", y = "RMSE")+
  scale_y_continuous(breaks=seq(-0.5, 5, 0.5), limits=c(-0.5, 5)) +
  #facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  #facet_wrap(~nrep, ncol=4)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")

ggsave(filename = "output/Figure_sigma2_phylo_rmse.png", width = 10, height = 7)





### 5. Residual variance estimate -----------------------------------


options(digits = 1, scipen = 5)
# log10 values 
s2.e_plot <- res1 |> 
  #filter(model != "INLA") |>
  ggplot(aes(x=factor(model), y=s2_resid, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) +  
  labs(title=TeX("Log10 $\\sigma^2_e$ estimate"),x="Package", y = TeX("Log10 $\\sigma^2_e$ estimate"))+
  scale_y_log10() +
  geom_hline(yintercept=0, colour="darkgray", linewidth=0.5)+
  #facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")

ggsave(filename = "output/Fig.test_sigma2.e_log10.png", width = 10, height = 7)


options(digits = 3, scipen = 5)
# bias 
s2.e_plot_bias <- res1 |> 
  #filter(model != "INLA") |>
  ggplot(aes(x=factor(model), y=s2_resid_bias, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title=TeX("Bias $\\sigma^2_e$ estimate"),x="Package", y = TeX("Bias $\\sigma^2_e$ estimate"))+
  geom_hline(yintercept=0, colour="darkgray", linewidth=0.5)+
  #facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")

ggsave(filename = "output/Figure_sigma2_resid_bias.png", width = 10, height = 7)


# RMSE
s2.e_plot_rmse <- res1 |> 
  #filter(model != "INLA") |>
  ggplot(aes(x=factor(model), y=rmse_s2_resid, fill=model, color=model)) + 
  geom_violin()+
  scale_color_manual(values=c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373")) +
  scale_fill_manual(values=alpha(c("#56B4E9", "#7AB47C", "#FBBF24", "#B47AA5", "#E57373"), 0.4)) + 
  labs(title="Plot of RMSE in residual variance estimates",x="Package", y = "RMSE")+
  #facet_wrap(~species_lab, ncol=3, labeller=label_parsed)+
  #facet_wrap(~nrep, ncol=4)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="none")

ggsave(filename = "output/Figure_sigma2_resid_rmse.png", width = 10, height = 7)




# Combine plots of variance estimates
fig_variance <- s2.sp_plot_bias + 
  s2.p_plot_bias + 
  s2.e_plot_bias + 
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(tag_levels='A',
                  theme = theme(plot.background = element_rect(fill = "white", colour = NA)))


ggsave(filename = "output/Figure_variance_estimates_bias.png", width = 10, height = 12)






######### Derive Monte Carlo Standard Errors (MCSE) of simulations ############


## number of simulation iterations of filtering results
nsim1 <- nrow(res1)/5
nsim2 <- nrow(res2)/3




# MCSE for b1 estimate
res1$mcse_b1 <- res1$rmse_b1 / sqrt(64000)



