rm(list=ls())

# Load packages -----------------------------------------

library("pacman")
p_load(parallel, ape, MCMCglmm, brms, phyr, plyr, dplyr, R.utils, glmmTMB, INLA,
       broom.mixed, ggstance, MASS, stringr, remotes, bayestestR, performance)


######################
#### Function ------------------------------------------------------
######################

#function to get model run time and handle any errors or warnings
run_model_safely <- function(expr) {
  warn <- err <- NULL
  rtime <- NA
  value <- withCallingHandlers(
    tryCatch({
      rtime <- unname(system.time(val <- eval(expr))["elapsed"])
      val
    }, error = function(e) {
      err <<- e
      NULL
    }), warning = function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value = value, warning = warn, error = err, rtime = rtime)
}



######################
#### Sim set-up ------------------------------------------------------
######################

### load job array
tab <- read.csv("job_array_sript1a.csv")
#tab <- tab1a

### get job number from pbs script
job <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
#job <- 1 for testing

### get parameters for current job
name <- tab$name[tab$job_number == job] 
scen <- tab$scenario[tab$job_number == job] 
seed <- tab$sim[tab$job_number == job]
dir <- tab$save_location[tab$job_number == job]
k.species <- tab$k.species[tab$job_number == job]  
n.reps <- tab$n.reps[tab$job_number == job]  
sigma2.s <- tab$sigma2.s[tab$job_number == job]   
sigma2.p <- tab$sigma2.p[tab$job_number == job] 
sigma2.e <- tab$sigma2.e[tab$job_number == job] 


# Set up dataframe to store results from all iterations
res <- data.frame(
  model = character(),
  species_size = integer(),
  sample_size = integer(),
  iteration = integer(),
  b0 = numeric(),
  b1 = numeric(),
  var.u = numeric(),
  var.p = numeric(),
  var.e = numeric(),
  run_time = numeric(),
  mu = numeric(),
  mu_ci_low = numeric(),
  mu_ci_high = numeric(),
  s2_sp = numeric(),
  s2_phylo = numeric(),
  s2_resid = numeric(),
  mu_bias = numeric(),
  mu_mse = numeric(),
  mu_cov = logical(),
  mu_ci_width = numeric(),
  s2_sp_bias = numeric(),
  s2_phylo_bias = numeric(),
  s2_resid_bias = numeric(),
  s2_sp_mse = numeric(),
  s2_phylo_mse = numeric(),
  s2_resid_mse = numeric(),
  stringsAsFactors = FALSE
)





######################
#### Simulate data ------------------------------------------------------
######################

### set seed
set.seed(seed)

### set up k.obs (number of observations per species)
if (is.null(n.reps)) {
  k.obs <- round(rbeta(n = k.species, shape1 = 1.25, shape2 = 3) * 99) + 1
} else {
  k.obs <- n.reps
}

### total number of observations
n <- k.obs*k.species

### create species id variable
sp.id <- rep(seq_len(k.species), times=k.obs)

### simulate simple dataframe with covariate (x) variable
x <- runif(n, 10, 20) 
dat <- data.frame(obs = 1:n, x = x, species = sp.id)

### simulate tree and obtain phylo matrix
tree <- rtree(k.species, tip.label = seq_len(k.species))
tree <- compute.brlen(tree, power=1) 
phylo.mat <- vcv(tree, corr = TRUE) ## we want correlation matrix (bounded by -1 and 1)
phylo.mat <- phylo.mat[order(as.numeric(rownames(phylo.mat))), order(as.numeric(rownames(phylo.mat)))]


### Simulate response variable (phen) based on cofactor and phylogenetic matrix
u.s <- rnorm(k.species, 0, sqrt(sigma2.s))[sp.id]
u.p <- mvrnorm(1, mu=rep(0, k.species), Sigma=sigma2.p*phylo.mat)[sp.id]
b0 <- 1
b1 <- 1.5
ei <- rnorm(n, 0, sqrt(sigma2.e))

### get estimates of y
yi <- b0 + b1*x + u.s + u.p + ei

### append all to dataframe
dat <- cbind(dat, u.s, u.p, ei, b0, b1, yi)
dat$species <- factor(dat$species) # format species variable for models
dat$phylo <- dat$species # phylo ID variable (same as species)
dat$sp <- dat$species # create sp variable (for phyr)
dat$g <- 1 # add variable g constant (for glmmTMB)

# save simulated data as R file
save(list = "dat", file = paste0("data/simdat_", job, ".RDATA"))





##################
#### Run models ----------------------------------------------------------------------------------
#################


### phyr model ###

# pglmm model - (1 | sp__) will construct two random terms, one with phylogenetic covariance 
# matrix and another with non-phylogenetic (identity) matrix. 
# --> Phylogenetic correlations can be dropped by removing the __ underscores.
res_phyr <- run_model_safely({
  pglmm(yi ~ x + (1|sp__),
        cov_ranef = list(sp = tree),
        data = dat,
        REML = TRUE)
})

model_phyr <- res_phyr$value
model_phyr_warning <- res_phyr$warning
model_phyr_error <- res_phyr$error
time.phyr <- res_phyr$rtime
conv_phyr <- !is.null(model_phyr)




### glmmTMB model ###
# allows for same random effect names for propto
res_tmb <- run_model_safely({
  glmmTMB(yi ~ x + (1|species) + propto(0 + species|g, phylo.mat),
          data = dat,
          REML = TRUE)
})

model_glmmTMB <- res_tmb$value
model_glmmTMB_error <- res_tmb$error
model_glmmTMB_warning <- res_tmb$warning
time.glmmTMB <- res_tmb$rtime
conv_tmb <- if (!is.null(model_glmmTMB)) model_glmmTMB$sdr$pdHess else NA



### brms model ###
# brms does not allow for duplicated group-level effects
res_brms <- run_model_safely({
  brm(yi ~ x + (1|species) + (1|gr(phylo, cov = phylo.mat)), #phylo.mat is the correlation matrix
      data = dat,
      family = gaussian(),
      chains = 4, iter = 2000, # defaults
      cores = 4, # equal to numebr of chains
      data2 = list(phylo.mat = phylo.mat))
})


model_brm <- res_brms$value
model_brm_error <- res_brms$error
model_brm_warning <- res_brms$warning
time.brms <- res_brms$rtime
# make a warning if at least one Rhat value is above 1.01 (Vehtari et al, 2021)
conv_brms <- if (!is.null(model_brm)) max(rhat(model_brm)) < 1.01 else NA
# store smallest ESS of model
ess_brms <- if (!is.null(model_brm)) {
  min(c(
    effective_sample(model_brm, effects = "fixed")$ESS,
    effective_sample(model_brm, effects = "random")$ESS
  ), na.rm = TRUE)
} else NA




### MCMCglmm model ###

# get precision phylo matrix and order rows
phylo.prec.mat <- inverseA(tree, nodes = "TIPS", scale = TRUE)$Ainv
phylo.prec.mat <- phylo.prec.mat[order(as.numeric(rownames(phylo.prec.mat))),
                                 order(as.numeric(rownames(phylo.prec.mat)))]

# set recommended priors with two random effects
prior <- list(G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000), 
                     G2=list(V=1,nu=1,alpha.mu=0,alpha.V=1000)),
              R=list(V=1,nu=0.02))

# run model
res_mcmc <- run_model_safely({
  MCMCglmm(yi ~ x,
           random = ~species + phylo,
           family = "gaussian",
           ginverse = list(phylo = phylo.prec.mat),
           prior = prior,
           data = dat,
           nitt = 13000, burnin = 3000, thin = 10)
})


model_mcmc <- res_mcmc$value
model_mcmc_error <- res_mcmc$error
model_mcmc_warning <- res_mcmc$warning
time.mcmc <- res_mcmc$rtime

#check convergence using Heidelberger diagnostic
conv_mcmc <- if (!is.null(model_mcmc)) {
  fullchain <- cbind(as.mcmc(model_mcmc$Sol), as.mcmc(model_mcmc$VCV))
  all(heidel.diag(fullchain)[,1] == 1)
} else NA
# store smallest ESS of model
ess_mcmc <- if (!is.null(model_mcmc)) {
  min(c(
    effective_sample(model_mcmc, effects = "fixed")$ESS,
    effective_sample(model_mcmc, effects = "random")$ESS
  ), na.rm = TRUE)
} else NA



### INLA model ###
# only allows for one covariate name per f()-term

# set up recommended penalizing complexity priors 
pcprior = list(prec = list(prior="pc.prec", param = c(20, 0.1)))

# run model (use generic0)
res_inla <- run_model_safely({
  inla(yi ~ x +
         f(species, model = "iid") +
         f(phylo, model = "generic0",
           Cmatrix = phylo.prec.mat,
           hyper = pcprior),
       family = "gaussian", data = dat)
})
model_inla <- res_inla$value
model_inla_error <- res_inla$error
model_inla_warning <- res_inla$warning
time.inla <- res_inla$rtime
conv_inla <- !is.null(model_inla)  #theres no convergence metric but check that it ran




# save all models objects
save(list = c("model_phyr", "model_glmmTMB", "model_brm", "model_mcmc", "model_inla"),
     file = paste0("models/models_", job, ".RDATA"))





###############
#### Store model estimates -------------------------------------------
###############

#------------- Fixed effect results: estimate + SE 

# phyr
coefs_phyr <- if (!is.null(model_phyr)) {
  cf <- as.data.frame(fixef(model_phyr))
  cf$conf.low[2] <- cf$Value[2] - cf$Std.Error[2] * 1.96
  cf$conf.high[2] <- cf$Value[2] + cf$Std.Error[2] * 1.96
  cf
} else data.frame(Value = NA, Std.Error = NA, conf.low = NA, conf.high = NA)


# glmmTMB
coefs_tmb <- if (!is.null(model_glmmTMB)) {
  as.data.frame(confint(model_glmmTMB, parm = "beta_"))
} else data.frame(Estimate = NA, `2.5 %` = NA, `97.5 %` = NA)

# brms
coefs_brm <- if (!is.null(model_brm)) {
  as.data.frame(tidy(model_brm, effects = "fixed", conf.int = TRUE))
} else data.frame(estimate = NA, conf.low = NA, conf.high = NA)

# MCMCglmm 
coefs_mcmc <- if (!is.null(model_mcmc)) {
  as.data.frame(tidy(model_mcmc, effects = "fixed", conf.int = TRUE))
} else data.frame(estimate = NA, conf.low = NA, conf.high = NA)

# INLA
coefs_inla <- if (!is.null(model_inla)) {
  as.data.frame(model_inla$summary.fixed)
} else data.frame(mean = NA, `0.025quant` = NA, `0.975quant` = NA)




#---- Random effect variance results: estimate + SE 


# get phyr variance estimates (model output is on SD scale)
sigma2_phyr <- if (!is.null(model_phyr)) {
  data.frame(
    model = "phyr",
    group = c("phylo", "species", "Residual"),
    term = "var",
    estimate = c(model_phyr$ss[2]^2, #phylogenetic 
                 model_phyr$ss[1]^2, #non-phylogenetic
                 model_phyr$ss[3]^2),#residual 
    std.error = NA, conf.low = NA, conf.high = NA  ##currently not available to my knowledge, could bootstrap tho..
  )
} else data.frame(
  model = "phyr",
  group = c("phylo", "species", "Residual"),
  term = "var",
  estimate = NA, std.error = NA, conf.low = NA, conf.high = NA
)



# get glmmTMB variance estimates (by default on standard deviation scale)
sigma2_tmb <- if (!is.null(model_glmmTMB)) {
  re_tmb <- as.data.frame(confint(model_glmmTMB, parm = "theta_"))
  sd_est <- re_tmb$Estimate
  sd_se <- (re_tmb$`97.5 %` - re_tmb$`2.5 %`) / (2 * 1.96)
  var_est <- sd_est^2
  var_se <- 2 * sd_est * sd_se
  var_ci_low <- var_est - 1.96 * var_se
  var_ci_high <- var_est + 1.96 * var_se
  
  resid_sd <- sigma(model_glmmTMB)
  resid_var <- resid_sd^2
  resid_ci <- confint(model_glmmTMB, parm = "sigma")
  resid_sd_se <- (resid_ci[2] - resid_ci[1]) / (2 * 1.96)
  resid_var_se <- 2 * resid_sd * resid_sd_se
  resid_var_ci <- c(resid_var - 1.96 * resid_var_se, resid_var + 1.96 * resid_var_se)
  
  data.frame(
    model = "glmmTMB",
    group = c("species", "phylo", "Residual"),
    term = "var",
    estimate = c(var_est, resid_var),
    std.error = c(var_se, resid_var_se),
    conf.low = c(var_ci_low, resid_var_ci[1]),
    conf.high = c(var_ci_high, resid_var_ci[2])
  )
} else data.frame(
  model = "glmmTMB",
  group = c("species", "phylo", "Residual"),
  term = "var",
  estimate = NA, std.error = NA, conf.low = NA, conf.high = NA
)





# get brms random effect variance estimates (standard deviation scale)
sigma2_brms <- if (!is.null(model_brm)) {
  tidy(model_brm, effects = "ran_pars") |>
    mutate(model = "brms",
           term = "var",
           estimate = estimate^2) |>  # convert SD to variance
    select(model, group, term, estimate, std.error, conf.low, conf.high)
} else data.frame(
  model = "brms",
  group = c("species", "phylo", "Residual"),
  term = "var",
  estimate = NA, std.error = NA, conf.low = NA, conf.high = NA
)



# get MCMCglmm random effect estimates (variance scale)
sigma2_mcmc <- if (!is.null(model_mcmc)) {
  tidy(model_mcmc, effects = "ran_pars", conf.int = TRUE) |>
    mutate(model = "MCMCglmm",
           group = str_replace(group, "animal", "phylo")) |> 
    select(model, group, term, estimate, std.error, conf.low, conf.high)
} else data.frame(
  model = "MCMCglmm",
  group = c("species", "phylo", "Residual"),
  term = "var",
  estimate = NA, std.error = NA, conf.low = NA, conf.high = NA
)


# get INLA random effect estimates (precision scale i.e. inverse variance) using 1000 posterior samples
## 1/fit.inla$hyperpar - before i used this but the mean of the inveerse is not the same as the inverse of the means (biased when posterio is skewed?)
sigma2_inla <- if (!is.null(model_inla)) {
  prec_samples <- inla.hyperpar.sample(1000, model_inla, improve.marginals = TRUE)
  var_samples <- as.data.frame(1 / prec_samples)  # convert precision to variance
  data.frame(
    model = "INLA",
    group = c("Residual", "species", "phylo"),
    term = "var",
    estimate = apply(var_samples, 2, mean),
    std.error = apply(var_samples, 2, sd),
    conf.low = apply(var_samples, 2, quantile, probs = 0.025),
    conf.high = apply(var_samples, 2, quantile, probs = 0.975)
  )
} else data.frame(
  model = "INLA",
  group = c("Residual", "species", "phylo"),
  term = "var",
  estimate = NA, std.error = NA, conf.low = NA, conf.high = NA
)




# merge fixed results together
s2 <- as.data.frame(rbind(sigma2_phyr,
                          sigma2_tmb,
                          sigma2_brms,
                          sigma2_mcmc, 
                          sigma2_inla))

# get subsets for each group
s2_phylo <- s2 %>% filter(group=="phylo")
s2_sp <- s2 %>% filter(group=="species")
s2_res <- s2 %>% filter(group=="Residual")



#####################
## Save run results ---------------------------------------------
#####################

# Combine the results for this iteration
iter_results <- data.frame(
  model = c("phyr", "glmmTMB", "brms", "MCMCglmm", "INLA"),
  species_size = k.species,
  sample_size = n,
  scenario = scen,
  seed = seed,
  b0 = rep(b0, 5),
  b1 = rep(b1, 5),
  sigma2.s = rep(sigma2.s, 5),
  sigma2.p = rep(sigma2.p, 5),
  sigma2.e = rep(sigma2.e, 5),
  run_time = c(time.phyr, time.glmmTMB, time.brms, time.mcmc, time.inla),
  iter_results$conv <- c(conv_phyr, conv_tmb, conv_brms, conv_mcmc, conv_inla),
  iter_results$ess <- c(ess_phyr, ess_tmb, ess_brms, ess_mcmc, ess_inla),
  mu = c(coefs_phyr$Value[2],
         coefs_tmb$Estimate[2], 
         coefs_brm$estimate[2], 
         coefs_mcmc$estimate[2], 
         coefs_inla$mean[2]),
  mu_ci_low = c(coefs_phyr$conf.low[2],
                coefs_tmb$`2.5 %`[2],
                coefs_brm$conf.low[2], 
                coefs_mcmc$conf.low[2], 
                coefs_inla$`0.025quant`[2]),
  mu_ci_high = c(coefs_phyr$conf.high[2],
                 coefs_tmb$`97.5 %`[2],
                 coefs_brm$conf.high[2], 
                 coefs_mcmc$conf.high[2], 
                 coefs_inla$`0.975quant`[2]),
  s2_sp = s2_sp$estimate,
  s2_phylo = s2_phylo$estimate,
  s2_resid = s2_res$estimate, 
  s2_sp_se = s2_sp$std.error,
  s2_phylo_se = s2_phylo$std.error,
  s2_resid_se = s2_res$std.error,
  s2_sp_ci_low = s2_sp$conf.low,
  s2_phylo_ci_low = s2_phylo$conf.low,
  s2_resid_ci_low = s2_res$conf.low,
  s2_sp_ci_high = s2_sp$conf.high,
  s2_phylo_ci_high = s2_phylo$conf.high,
  s2_resid_ci_high = s2_res$conf.high,
  stringsAsFactors = FALSE
)


# Derive coverage, bias and CI width for the estimates
iter_results$mu_bias <- iter_results$mu - dat$b1[1]
iter_results$mu_mse <- (iter_results$mu - dat$b1[1])^2
iter_results$mu_cov <- iter_results$mu_ci_low < dat$b1[1] & iter_results$mu_ci_low > dat$b1[1]
iter_results$mu_ci_width <- iter_results$mu_ci_high - iter_results$mu_ci_low
iter_results$s2_sp_bias <- iter_results$s2_sp - iter_results$sigma2.s
iter_results$s2_phylo_bias <- iter_results$s2_phylo - iter_results$sigma2.p
iter_results$s2_resid_bias <- iter_results$s2_resid - iter_results$sigma2.e
iter_results$s2_sp_mse <- (iter_results$s2_sp - iter_results$sigma2.s)^2
iter_results$s2_phylo_mse <- (iter_results$s2_phylo - iter_results$sigma2.p)^2
iter_results$s2_resid_mse <- (iter_results$s2_resid - iter_results$sigma2.e)^2


# Save results as RDATA
#save(iter_results, file=sprintf("/srv/scratch/z5394590/phyloglmm/R/pilot/output/study1/res_iter%d_cond%d.RData", i, j))

# Save results as RDATA
save(results, file="/srv/scratch/z5394590/phyloglmm/R/pilot/output/results_sim_study1.RData")


