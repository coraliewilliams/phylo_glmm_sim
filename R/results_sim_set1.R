################### Load results #############################

library(ggplot2);library(devtools);library(cowplot);library(ggdark);theme_set(dark_theme_bw())
library(ggdist);library(tidyverse);library(dplyr)

# load results (all set 1)
load("~/PhD/3_PhyloTMB/sims/sim_results_set1a.RDATA")

results <- dat


################### Plot results #############################

# Plot run times vs packages (y-axis on log10 scale)
runtime_plot <- results |> 
  #filter(model %in% c("glmmTMB", "INLA", "brms", "MCMCglmm")) |> 
  ggplot(aes(x=model, y=run_time, fill=model, color=model)) + 
  geom_violin()+
  scale_y_log10()+
  scale_color_manual(values=c("#B47AA5", "#E69F00", "#56B4E9", "#7AB47C", "#FBBF24")) +
  scale_fill_manual(values=alpha(c("#B47AA5", "#E69F00", "#56B4E9", "#7AB47C", "#FBBF24"), 0.4)) + 
  labs(title="Plot of run time",x="Package", y = "Run time (seconds)")+
  geom_boxplot(width=0.1)+
  theme_classic()


ggsave(filename = "output/runtime_set1a_log10_runtime.png")



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


