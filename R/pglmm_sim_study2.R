rm(list=ls())

# Load packages -----------------------------------------

library("pacman")
p_load(parallel, ape, MCMCglmm, brms, phyr, plyr, dplyr, R.utils,
       broom.mixed, ggstance, MASS, stringr, remotes)

# build glmmTMB package from Maeve's prop2 branch (March 2023)
#remotes::install_github("mmcgillycuddy/glmmTMB", ref="prop2_covstruc", subdir="glmmTMB")
library(glmmTMB)
library(INLA)


# Functions -------------------------------------

# function to get run time 
rtime <- function(...) unname(system.time(capture.output(...))["elapsed"])


# Set up ----------------------------------------------------

### number of iterations
iters <- 50

### k.species - number of species
k.species <- c(250, 500, 1000)

### scenarios for random effects:
# set sigma2 for species and phylo effect with large difference to be able to pick up when they don't estimate properly
sigma2.n <- c(0.05)   # size of non-phylogenetic variance component
sigma2.p <- c(0.3)   # size of phylogenetic variance component
sigma2.e <- 0.2      # estimate level variance component


# set up model conditions data frame
all.conds <- expand.grid(k.species = k.species,
                         sigma2.n = sigma2.n, 
                         sigma2.p = sigma2.p,
                         sigma2.e = sigma2.e)



# Simulation  ----------------------------------------------------

# get the number of scenarios/conditions to test
n_cond <- length(all.conds$k.species)

# Set up sim_dat to store all simulated dataframes as a list
sim_dat <- vector("list", length = iters * n_cond)

# Set number of rows for each result block
n_models <- 3
n_rows <- n_models * iters

# Set up dataframe to store results from all iterations
res <- data.frame(
  model = character(n_rows),
  species_size = integer(n_rows),
  sample_size = integer(n_rows),
  iteration = integer(n_rows),
  b0 = numeric(n_rows),
  b1 = numeric(n_rows),
  var.u = numeric(n_rows),
  var.p = numeric(n_rows),
  var.e = numeric(n_rows),
  run_time = numeric(n_rows),
  mu = numeric(n_rows),
  mu_ci_low = numeric(n_rows),
  mu_ci_high = numeric(n_rows),
  s2_sp = numeric(n_rows),
  s2_phylo = numeric(n_rows),
  s2_resid = numeric(n_rows),
  mu_bias = numeric(n_rows),
  mu_mse = numeric(n_rows),
  mu_cov = logical(n_rows),
  mu_ci_width = numeric(n_rows),
  s2_sp_bias = numeric(n_rows),
  s2_phylo_bias = numeric(n_rows),
  s2_resid_bias = numeric(n_rows),
  s2_sp_mse = numeric(n_rows),
  s2_phylo_mse = numeric(n_rows),
  s2_resid_mse = numeric(n_rows),
  stringsAsFactors = FALSE
)

# Pre-allocate dataframes
results <- vector("list")
simdat <- vector("list")

# initialise 
li <- 1 # to save sim dataframe
i <- 1 # iterations
j <- 1 # conditions

# Loop through each conditions 
for (j in 1:n_cond) {
  
  cond  <- all.conds[j,]
  ri <- 1
  
  # Loop iters times
  for (i in 1:iters) {
    
    cat("scenario = ", j, ", iter = ", i, "\n")
    set.seed(123 + i)
    
    #### Simulate data ------------------------------------------------------
    
    ### k.obs - number of observations per species (left skewed distribution)
    k.obs <- round(rbeta(n = cond$k.species, shape1 = 1.5, shape2 = 3) * 99) + 1
    
    ### create species id variable
    sp.id <- rep(seq_len(cond$k.species), times=k.obs)
    
    ### total number of observations
    n <- length(sp.id)
    
    ### simulate simple dataframe with covariate (x) variable
    x <- runif(n, 10, 20) 
    dat <- data.frame(obs = 1:n, x = x, species = sp.id)
    
    ### simulate tree and obtain phylo matrix
    tree <- rtree(cond$k.species, tip.label = seq_len(cond$k.species))
    tree <- compute.brlen(tree, power=1) 
    phylo.mat <- vcv(tree, corr = TRUE) ##covariance matrix for meta-analysis, correlation
    phylo.mat <- phylo.mat[order(as.numeric(rownames(phylo.mat))), order(as.numeric(rownames(phylo.mat)))]
    
    
    ### Simulate response variable (phen) based on cofactor and phylogenetic matrix
    u.n <- rnorm(cond$k.species, 0, sqrt(cond$sigma2.n))[sp.id]
    u.p <- mvrnorm(1, mu=rep(0, cond$k.species), Sigma=cond$sigma2.p*phylo.mat)[sp.id]
    b0 <- 1
    b1 <- 1.5
    ei <- rnorm(n, 0, sqrt(sigma2.e))
    
    ### get estimates of y
    yi <- b0 + b1*x + u.n + u.p + ei
    
    ### append all to dataframe
    dat <- cbind(dat, u.n, u.p, ei, b0, b1, yi)
    
    # create phylo variable 
    dat$phylo <- dat$species
    # format species variable for models
    dat$species <- factor(dat$species)
    
    
    #### Run models ------------------------------------------------------
    
    
    ### glmmTMB model ###
    
    # add variable g
    dat$g <- 1
    # run model with propto
    time.glmmTMB <- rtime(model_glmmTMB <- glmmTMB(
      yi ~ x + (1|species) + propto(0 + species|g,phylo.mat),
      data = dat,
      REML=TRUE))
    
    fit.glmmTMB <- summary(model_glmmTMB)
    #VCV.TMB <- VarCorr(model_glmmTMB)
    #VCV.TMB$cond$g
    
    
    
    ### MCMCglmm model ###
    
    # get precision phylo matrix and order rows
    phylo.prec.mat <- inverseA(tree, nodes = "TIPS", scale = TRUE)$Ainv
    phylo.prec.mat <- phylo.prec.mat[order(as.numeric(rownames(phylo.prec.mat))), order(as.numeric(rownames(phylo.prec.mat)))]

    # set priors with two random effects
    prior <- list(G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000), 
                         G2=list(V=1,nu=1,alpha.mu=0,alpha.V=1000)),
                  R=list(V=1,nu=0.02))
    
    # run model
    time.mcmc <- rtime(model_mcmc <- MCMCglmm(
      yi~x, random=~species+phylo,
      family="gaussian", 
      ginverse=list(phylo=phylo.prec.mat),
      prior=prior, data=dat, nitt=13000,
      burnin=3000, thin=10))
    
    fit.mcmc <- summary(model_mcmc)
    
    
    
    ### INLA model ###
    
    # only allows for one covariate name per f()-term
    
    # set up penalizing complexity priors (FOR NOW DON'T ADD)
    #pcprior = list(prec = list(prior="pc.prec", param = c(20, 0.1)))
    
    # run model (use generic0)
    time.inla <- rtime(model_inla <- inla(
      yi ~ x + f(species, model = "iid") + 
        f(phylo, model = "generic0",
          Cmatrix = phylo.prec.mat),
      family = "gaussian",
      data = dat))
    
    fit.inla <- summary(model_inla)
    
    
    
    
    #### Get performance measures -----------------------------------------------

    
    #---- Fixed effect results: estimate + SE 
    
    # glmmTMB
    coefs_tmb <- as.data.frame(confint(model_glmmTMB, parm="beta_")) ### this is slow -- is it possible to speed up?
    # MCMCglmm 
    coefs_mcmc <- as.data.frame(tidy(model_mcmc, effects="fixed", conf.int=TRUE))
    # INLA
    coefs_inla <- as.data.frame(fit.inla$fixed)
    
    
    
    #---- Random effect results: estimate + SE 
    
    
    # get glmmTMB random effect estimates 
    re_tmb <- as.data.frame(confint(model_glmmTMB, parm="theta_")) ### this is super slow -- is it possible to speed up?
    species_tmb <- re_tmb[1, ]
    phylo_tmb <- re_tmb[2, ]
    # combine into dataframe
    sigma2_tmb <- data.frame(
      model = "glmmTMB",
      group = c("phylo", "species", "Residual"),
      term = "var",
      estimate = c(phylo_tmb$Estimate^2, species_tmb$Estimate^2, sigma(model_glmmTMB)^2),
      std.error = NA, # Replace with the calculated SE if available
      conf.low = c(phylo_tmb$`2.5 %`, species_tmb$`2.5 %`, NA), # Replace with the residual var CI if available
      conf.high = c(phylo_tmb$`97.5 %`, species_tmb$`97.5 %`, NA) # Replace with the residual var CI if available
    )
    
    
    # get MCMCglmm random effect estimates
    sigma2_mcmc <- tidy(model_mcmc, effects="ran_pars", conf.int=TRUE)
    sigma2_mcmc <- sigma2_mcmc %>%
      mutate(model="MCMCglmm") %>% 
      dplyr::select(model, group, term, estimate, std.error, conf.low, conf.high)
    
    
    # get INLA random effect estimates
    re_inla <- 1/fit.inla$hyperpar
    sigma2_inla <- data.frame(
      model = "INLA",
      group = c( "Residual", "species", "phylo"),
      term = "var",
      estimate = re_inla$mean,
      std.error = fit.inla$hyperpar$sd,
      conf.low = re_inla$`0.025quant`,
      conf.high = re_inla$`0.975quant`
    )
    
    
    # merge fixed results together
    s2 <- data.frame(rbind(sigma2_tmb,
                              sigma2_mcmc, 
                              sigma2_inla),
                     stringsAsFactors = FALSE)
    
    # get subsets for each group
    s2_phylo <- s2 %>% filter(group=="phylo")
    s2_sp <- s2 %>% filter(group=="species")
    s2_res <- s2 %>% filter(group=="Residual")
    
    
    
    #### Save run results ---------------------------------------------
    
    # Combine the results for this iteration
    iter_results <- data.frame(
      model = c("glmmTMB", "MCMCglmm", "INLA"),
      species_size = cond$k.species,
      sample_size = n,
      iteration = i,
      b0 = rep(dat$b0[1], 3),
      b1 = rep(dat$b1[1], 3),
      var.u = rep(cond$sigma2.n, 3),
      var.p = rep(cond$sigma2.p, 3),
      var.e = rep(cond$sigma2.e, 3),
      run_time = c(time.glmmTMB, time.mcmc, time.inla),
      mu = c(coefs_tmb$Estimate[2], 
             coefs_mcmc$estimate[2], 
             coefs_inla$mean[2]),
      mu_ci_low = c(coefs_tmb$`2.5 %`[2],
                    coefs_mcmc$conf.low[2], 
                    coefs_inla$`0.025quant`[2]),
      mu_ci_high = c(coefs_tmb$`97.5 %`[2],
                     coefs_mcmc$conf.high[2], 
                     coefs_inla$`0.975quant`[2]),
      s2_sp = unname(s2_sp$estimate),
      s2_phylo = unname(s2_phylo$estimate),
      s2_resid = unname(s2_res$estimate), 
      stringsAsFactors = FALSE
    )
    
    
    # Derive coverage, bias and CI width for the estimates
    iter_results$mu_bias <- iter_results$mu - dat$b1[1]
    iter_results$mu_mse <- (iter_results$mu - dat$b1[1])^2
    iter_results$mu_cov <- iter_results$mu_ci_low < dat$b1[1] & iter_results$mu_ci_low > dat$b1[1]
    iter_results$mu_ci_width <- iter_results$mu_ci_high - iter_results$mu_ci_low
    iter_results$s2_sp_bias <- iter_results$s2_sp - iter_results$var.u
    iter_results$s2_phylo_bias <- iter_results$s2_phylo - iter_results$var.p
    iter_results$s2_resid_bias <- iter_results$s2_resid - iter_results$var.e
    iter_results$s2_sp_mse <- (iter_results$s2_sp - iter_results$var.u)^2
    iter_results$s2_phylo_mse <- (iter_results$s2_phylo - iter_results$var.p)^2
    iter_results$s2_resid_mse <- (iter_results$s2_resid - iter_results$var.e)^2
    
    
    # Save results as RDATA
    save(iter_results, file=sprintf("/srv/scratch/z5394590/phyloglmm/R/pilot/output/study2/res_iter%d_cond%d.RData", i, j))

    # Add the iteration results of models to the overall results
    res[ri:(ri + n_models - 1), ] <- iter_results
    ri <- ri + n_models
    li <- li + 1
    
  }
  
  # Combine the results for this condition
  results[[j]] <- res
}

# Save results as RDATA
save(results, file="/srv/scratch/z5394590/phyloglmm/R/pilot/output/results_sim_study2.RData")


