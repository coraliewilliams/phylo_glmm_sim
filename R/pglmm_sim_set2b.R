rm(list=ls())

# Load packages -----------------------------------------

library("pacman")
p_load(parallel, ape, MCMCglmm, phyr, plyr, dplyr, R.utils, glmmTMB, INLA,
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
tab <- read.csv("job_array_set2b.csv")
#tab <- tab2b

### get job number from pbs script
job <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
#job <- 1

### get parameters for current job
name <- tab$name[tab$job_number == job] 
scen <- tab$scenario[tab$job_number == job] 
seed <- tab$sim[tab$job_number == job]
dir <- tab$save_location[tab$job_number == job]
k.species <- tab$k.species[tab$job_number == job]  
n.reps <- tab$n.reps[tab$job_number == job]  
sigma2.p <- tab$sigma2.p[tab$job_number == job] 
sigma2.e <- tab$sigma2.e[tab$job_number == job] 


# Set up dataframe to store results of current iteration
result <- data.frame(
  model = character(),
  model_error = character(),
  model_warning = character(),
  convergence = logical(),
  ESS = numeric(),
  species_size = integer(),
  sample_size = integer(),
  scenario = character(),
  seed = integer(),
  b0 = numeric(),
  b1 = numeric(),
  sigma2.p = numeric(),
  sigma2.e = numeric(),
  run_time = numeric(),
  mu = numeric(),
  mu_ci_low = numeric(),
  mu_ci_high = numeric(),
  s2_phylo = numeric(),
  s2_resid = numeric(),
  s2_phylo_se = numeric(),
  s2_resid_se = numeric(),
  s2_phylo_ci_low = numeric(),
  s2_resid_ci_low = numeric(),
  s2_phylo_ci_high = numeric(),
  s2_resid_ci_high = numeric(),
  mu_bias = numeric(),
  mu_mse = numeric(),
  mu_cov = logical(),
  mu_ci_width = numeric(),
  s2_phylo_bias = numeric(),
  s2_resid_bias = numeric(),
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
k.obs <- n.reps

### create species id variable
sp.id <- rep(seq_len(k.species), times=k.obs)

### total number of observations
n <- length(sp.id)

### simulate simple dataframe with covariate (x) variable
x <- runif(n, 10, 20) 
x <- scale(x) # scale x to have mean 0 and sd 1)
dat <- data.frame(obs = 1:n, x = x, species = sp.id)

### simulate tree and obtain phylo matrix
tree <- rtree(k.species, tip.label = seq_len(k.species))
tree <- compute.brlen(tree, power=1) 
phylo.mat <- vcv(tree, corr = TRUE) ## we want correlation matrix (bounded by -1 and 1)
phylo.mat <- phylo.mat[order(as.numeric(rownames(phylo.mat))), order(as.numeric(rownames(phylo.mat)))]

### get precision phylo matrix and order rows (for MCMCglmm and INLA)
phylo.prec.mat <- inverseA(tree, nodes = "TIPS", scale = TRUE)$Ainv
phylo.prec.mat <- phylo.prec.mat[order(as.numeric(rownames(phylo.prec.mat))),
                                 order(as.numeric(rownames(phylo.prec.mat)))]


### Simulate response variable (phen) based on cofactor and phylogenetic matrix
u.p <- mvrnorm(1, mu=rep(0, k.species), Sigma=sigma2.p*phylo.mat)[sp.id]
b0 <- 1
b1 <- 1.5
ei <- rnorm(n, 0, sqrt(sigma2.e))

### get estimates of y
yi <- b0 + b1*x + u.p + ei

### append all to dataframe
dat <- cbind(dat, u.p, ei, b0, b1, yi)
dat$phylo <- dat$species # phylo ID variable (same as species) - needs to be numeric to work with INLA
dat$species <- factor(dat$species) # format species variable for models
dat$sp <- dat$species # create sp variable (for phyr)
dat$g <- 1 # add variable g constant (for glmmTMB)



# save simulated data as R file
save(list = c("dat", "phylo.mat", "phylo.prec.mat"), file = paste0("data/b/simdat_", job, ".RDATA"))





##################
#### Run models ----------------------------------------------------------------------------------
#################



### glmmTMB model ###
# allows for same random effect names for propto
res_tmb <- run_model_safely({
  glmmTMB(yi ~ x + propto(0 + species|g, phylo.mat),
          data = dat,
          REML = TRUE)
})

model_glmmTMB <- res_tmb$value
model_glmmTMB_error <- if (is.null(res_tmb$error)) NA else res_tmb$error$message
model_glmmTMB_warning <- if (is.null(res_tmb$warning)) NA else res_tmb$warning$message
time.glmmTMB <- res_tmb$rtime
conv_tmb <- if (!is.null(model_glmmTMB)) model_glmmTMB$sdr$pdHess else FALSE





### MCMCglmm model ###

# set recommended priors with two random effects
prior <- list(G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000)),
              R=list(V=1,nu=0.02))

# run model
res_mcmc <- run_model_safely({
  MCMCglmm(yi ~ x,
           random = ~phylo,
           family = "gaussian",
           ginverse = list(phylo = phylo.prec.mat),
           prior = prior,
           data = dat,
           nitt = 303000, # increase default by x30
           burnin = 3000, # default
           thin = 10) # default
})


model_mcmc <- res_mcmc$value
model_mcmc_error <- if (is.null(res_mcmc$error)) NA else res_mcmc$error$message
model_mcmc_warning <- if (is.null(res_mcmc$warning)) NA else res_mcmc$warning$message
time.mcmc <- res_mcmc$rtime

#check convergence using Heidelberger diagnostic
conv_mcmc <- if (!is.null(model_mcmc)) {
  fullchain <- cbind(as.mcmc(model_mcmc$Sol), as.mcmc(model_mcmc$VCV))
  all(heidel.diag(fullchain)[,1] == 1)
} else FALSE
# store smallest ESS of model
ess_mcmc <- if (!is.null(model_mcmc)) {
  min(c(
    effective_sample(model_mcmc, effects = "fixed")$ESS,
    effective_sample(model_mcmc, effects = "random")$ESS
  ), na.rm = TRUE)
} else NA






### INLA model ###
# only allows for one covariate name per f()-term

# run model (use generic0)
res_inla <- run_model_safely({
  inla(yi ~ x + 
         f(phylo, ## this needs to be a numeric to work
           model = "generic0",
           Cmatrix = phylo.prec.mat),
       family = "gaussian",
       data = dat)
})

model_inla <- res_inla$value
model_inla_error <- if (is.null(res_inla$error)) NA else res_inla$error$message
model_inla_warning <- if (is.null(res_inla$warning)) NA else res_inla$warning$message
time.inla <- res_inla$rtime
conv_inla <- model_inla$ok



###############
#### Store model estimates -------------------------------------------
###############

#------------- Fixed effect results: estimate + SE 


# glmmTMB
coefs_tmb <- if (!is.null(model_glmmTMB)) {
  as.data.frame(confint(model_glmmTMB, parm = "beta_"))
} else data.frame(Estimate = NA, `2.5 %` = NA, `97.5 %` = NA)

# MCMCglmm 
coefs_mcmc <- if (!is.null(model_mcmc)) {
  as.data.frame(tidy(model_mcmc, effects = "fixed", conf.int = TRUE))
} else data.frame(estimate = NA, conf.low = NA, conf.high = NA)

# INLA
coefs_inla <- if (!is.null(model_inla)) {
  as.data.frame(model_inla$summary.fixed)
} else data.frame(mean = NA, `0.025quant` = NA, `0.975quant` = NA, check.names = F)




#---- Random effect variance results: estimate + SE 

# get glmmTMB variance estimates (standard deviation scale)
sigma2_tmb <- if (!is.null(model_glmmTMB)) {
  re_tmb <- as.data.frame(confint(model_glmmTMB, parm = "theta_"))
  sd_est <- re_tmb$Estimate
  sd_se <- (re_tmb$`97.5 %` - re_tmb$`2.5 %`) / (2 * 1.96)
  var_est <- sd_est^2
  var_se <- 2 * sd_est * sd_se ##using delta method to get variance SE (assuming normality and symmetry)
  var_ci_low <- var_est - 1.96 * var_se
  var_ci_high <- var_est + 1.96 * var_se
  
  resid_sd <- sigma(model_glmmTMB)
  resid_var <- resid_sd^2
  resid_ci <- confint(model_glmmTMB, parm = "sigma")
  resid_sd_se <- (resid_ci[2] - resid_ci[1]) / (2 * 1.96)
  resid_var_se <- 2 * resid_sd * resid_sd_se ##using delta method to get variance SE (assuming normality and symmetry)
  resid_var_ci <- c(resid_var - 1.96 * resid_var_se, resid_var + 1.96 * resid_var_se)
  
  data.frame(
    model = "glmmTMB",
    group = c("phylo", "Residual"),
    term = "var",
    estimate = c(var_est, resid_var),
    std.error = c(var_se, resid_var_se),
    conf.low = c(var_ci_low, resid_var_ci[1]),
    conf.high = c(var_ci_high, resid_var_ci[2])
  )
} else data.frame(
  model = "glmmTMB",
  group = c("phylo", "Residual"),
  term = "var",
  estimate = NA, std.error = NA, conf.low = NA, conf.high = NA
)



# get MCMCglmm random effect estimates (variance scale)
sigma2_mcmc <- if (!is.null(model_mcmc)) {
  tidy(model_mcmc, effects = "ran_pars", conf.int = TRUE) |>
    mutate(model = "MCMCglmm",
           group = str_replace(group, "animal", "phylo")) |> 
    dplyr::select(model, group, term, estimate, std.error, conf.low, conf.high)
} else data.frame(
  model = "MCMCglmm",
  group = c("phylo", "Residual"),
  term = "var",
  estimate = NA, std.error = NA, conf.low = NA, conf.high = NA
)




# get INLA random effect estimates (precision scale i.e. inverse variance) using 1000 posterior samples
## 1/fit.inla$hyperpar - before i used this but the mean of the inveerse is not the same as the inverse of the means (biased when posterio is skewed?)
sigma2_inla <- if (!is.null(model_inla)) {
  prec_samples <- inla.hyperpar.sample(1000, model_inla, improve.marginals = TRUE) ##is it all the samples?
  var_samples <- as.data.frame(1 / prec_samples)  # convert precision to variance
  data.frame(
    model = "INLA",
    group = c("Residual", "phylo"),
    term = "var",
    estimate = apply(var_samples, 2, mean),
    std.error = apply(var_samples, 2, sd),
    conf.low = apply(var_samples, 2, quantile, probs = 0.025),
    conf.high = apply(var_samples, 2, quantile, probs = 0.975)
  )
} else data.frame(
  model = "INLA",
  group = c("Residual", "phylo"),
  term = "var",
  estimate = NA, std.error = NA, conf.low = NA, conf.high = NA
)




# merge fixed results together
s2 <- as.data.frame(bind_rows(sigma2_tmb, 
                              sigma2_mcmc, 
                              sigma2_inla))

# get subsets for each group
s2_phylo <- s2 %>% filter(group=="phylo")
s2_res <- s2 %>% filter(group=="Residual")



#####################
## Save run results ---------------------------------------------
#####################

# Combine the results for this iteration
result <- data.frame(
  model = c("glmmTMB", "MCMCglmm", "INLA"),
  model_error = c(model_glmmTMB_error, model_mcmc_error, model_inla_error),
  model_warning = c(model_glmmTMB_warning, model_mcmc_warning, model_inla_warning),
  convergence = c(conv_tmb, conv_mcmc, conv_inla),
  ESS = c(NA, ess_mcmc, NA),
  species_size = rep(k.species, 3),
  sample_size = rep(n, 3),
  scenario = rep(scen, 3),
  seed = rep(seed, 3),
  b0 = rep(b0, 3),
  b1 = rep(b1, 3),
  sigma2.p = rep(sigma2.p, 3),
  sigma2.e = rep(sigma2.e, 3),
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
  s2_phylo = s2_phylo$estimate,
  s2_resid = s2_res$estimate, 
  s2_phylo_se = s2_phylo$std.error,
  s2_resid_se = s2_res$std.error,
  s2_phylo_ci_low = s2_phylo$conf.low,
  s2_resid_ci_low = s2_res$conf.low,
  s2_phylo_ci_high = s2_phylo$conf.high,
  s2_resid_ci_high = s2_res$conf.high,
  stringsAsFactors = FALSE
)


# Derive coverage, bias and CI width for the estimates
result$mu_bias <- result$mu - dat$b1[1]
result$mu_mse <- (result$mu - dat$b1[1])^2
result$mu_cov <- result$mu_ci_low < dat$b1[1] & result$mu_ci_high > dat$b1[1]
result$mu_ci_width <- result$mu_ci_high - result$mu_ci_low
result$s2_phylo_bias <- result$s2_phylo - result$sigma2.p
result$s2_resid_bias <- result$s2_resid - result$sigma2.e
result$s2_phylo_mse <- (result$s2_phylo - result$sigma2.p)^2
result$s2_resid_mse <- (result$s2_resid - result$sigma2.e)^2


# Save results as RDATA
save(list = "result", file = paste0("results/b/res_", job, ".RDATA"))

