###############################################################################
# SIMULATION PARAMETERS AND JOB ARRAY SET UP
###############################################################################


###############
## SCRIPT 1a ## --------------------------------------------------------------------------
###############

############### Unbalanced + balanced design 
name <- "script1a"

# model parameters
k.species <- c(25, 50, 100)  # number of species
n.reps <- c(5, 10, 30, NULL) # number of repeated measures per specie
sigma2.s <- c(0.05, 0.25)    # variance of non-phylogenetic variance component (species level)
sigma2.p <- c(0.05, 0.25)    # variance of phylogenetic variance component (species level)
sigma2.e <- 0.2              # variance of sampling error (estimate level level)


# number of replications per scenario
repl <- 4000
sim <- rep(1:repl)

# make table of all scenarios
tab1a <- expand.grid(sim = sim, name = name, n.reps = n.reps, k.species = k.species,
                     sigma2.s = sigma2.s, sigma2.p = sigma2.p, sigma2.e = sigma2.e)

# add columns for to store results and job number
tab1a$save_location <- rep("/srv/scratch/z5394590/", each=nrow(tab1a))
tab1a$job_number <- c(1:nrow(tab1a))
conds <- length(k.species) * length(sigma2.s) * length(sigma2.p) * length(sigma2.e)
tab1a$scenario <- rep(1:conds, each = repl)

# save as csv file
write.csv(tab1a, "output/job_array_script1a.csv", row.names = FALSE)



############### No repetitions 
name <- "script1b"

# model parameters
k.species <- c(25, 50, 100)  # number of species
n.reps <- 1 # one measures per specie
sigma2.s <- c(0.05, 0.25)    # variance of non-phylogenetic variance component (species level)
sigma2.p <- c(0.05, 0.25)    # variance of phylogenetic variance component (species level)


# number of replications per scenario
repl <- 4000
sim <- rep(1:repl)

# make table of all scenarios
tab1b <- expand.grid(sim = sim, name = name, n.reps = n.reps, k.species = k.species,
                     sigma2.s = sigma2.s, sigma2.p = sigma2.p)

# add columns for to store results and job number
tab1b$save_location <- rep("/srv/scratch/z5394590/", each=nrow(tab1b))
tab1b$job_number <- c(1:nrow(tab1a))
conds <- length(k.species) * length(sigma2.s) * length(sigma2.p) 
tab1b$scenario <- rep(1:conds, each = repl)

# save as csv file
write.csv(tab1b, "output/job_array_script1b.csv", row.names = FALSE)
