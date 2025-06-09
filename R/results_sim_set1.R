################### Load results #############################

library(ggplot2);library(devtools);library(cowplot);library(ggdark);theme_set(dark_theme_bw())
library(ggdist);library(tidyverse);library(dplyr)

# load results set A (repeated measures)
load("~/PhD/3_PhyloTMB/sims/sim_results_set1a.RDATA")
res.a <- dat

# load results set B (one measure per species)
load("~/PhD/3_PhyloTMB/sims/sim_results_set1b.RDATA")
res.b <- dat

options(digits = 4, scipen = 5)

################### Plot results #############################


### 1. Run time -------------------------------------

# Plot run times vs packages (y-axis on log10 scale)
runtime_plot <- res.a |> 
  #filter(model %in% c("glmmTMB", "INLA", "brms", "MCMCglmm")) |> 
  ggplot(aes(x=factor(model), y=run_time, fill=model, color=model)) + 
  geom_violin()+
  scale_y_log10()+
  scale_color_manual(values=c("#B47AA5", "#E69F00", "#56B4E9", "#7AB47C", "#FBBF24")) +
  scale_fill_manual(values=alpha(c("#B47AA5", "#E69F00", "#56B4E9", "#7AB47C", "#FBBF24"), 0.4)) + 
  labs(title="Plot of run time",x="Package", y = "Run time (seconds)")+
  facet_wrap(~species_size)+
  geom_boxplot(width=0.1)+
  theme_bw() +
  theme(legend.position="bottom")

ggsave(filename = "output/Figure_runtime.png", width = 10, height = 6)




runtime_plot <- results |> 
  filter(model %in% c("glmmTMB", "INLA", "phyr")) |> 
  ggplot(aes(x=factor(model), y=run_time,
             fill=factor(model), color=factor(model))) + 
  geom_violin()+
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#FBBF24")) +
  scale_fill_manual(values=alpha(c("#E69F00", "#56B4E9", "#FBBF24"), 0.4)) + 
  labs(title="Plot of run time (non MCMC models)",x="Package", y = "Run time (seconds)")+
  geom_boxplot(width=0.1)+
  theme_classic()

ggsave(filename = "output/runtime_nonMCMC.png")







### 2. Fixed effect estimate -----------------------------------

# Bias




# MSE






### 3. Non-phylo variance estimate -----------------------------------







### 4. Phylo variance estimate -----------------------------------







### 5. Residual variance estimate -----------------------------------


