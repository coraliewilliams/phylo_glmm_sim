###############################################################################
# SIMULATION PARAMETERS AND JOB ARRAY SET UP
###############################################################################


#############
## SET 1a  ## --------------------------------------------------------------------------
#############

############### Unbalanced + balanced design 
name <- "set1a"

# model parameters
k.species <- c(25, 50, 100)  # number of species
n.reps <- c(5, 10, 30, 999) # number of repeated measures per specie (999 refers to unbalanced design)
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
tab1a$save_location <- rep("/srv/scratch/z5394590/pglmm/set1a", each=nrow(tab1a))
tab1a$job_number <- c(1:nrow(tab1a))
conds <- length(k.species) * length(sigma2.s) * length(sigma2.p) * length(sigma2.e)
tab1a$scenario <- rep(1:conds, each = repl)

# save as csv file
write.csv(tab1a, "output/job_array_set1a.csv", row.names = FALSE)



############### No repeated measures 
name <- "set1b"

# model parameters
n.reps <- 1 # one measures per specie

# make table of all scenarios
tab1b <- expand.grid(sim = sim, name = name, n.reps = n.reps, k.species = k.species,
                     sigma2.e = sigma2.e, sigma2.p = sigma2.p)

# add columns for to store results and job number
tab1b$save_location <- rep("/srv/scratch/z5394590/", each=nrow(tab1b))
tab1b$job_number <- c(1:nrow(tab1b))
conds <- length(k.species) * length(sigma2.e) * length(sigma2.p) 
tab1b$scenario <- rep(1:conds, each = repl)

# save as csv file
write.csv(tab1b, "output/job_array_set1b.csv", row.names = FALSE)



#############
## SET 2a  ## --------------------------------------------------------------------------
#############

############### Unbalanced + balanced design 
name <- "set2a"

# model parameters
k.species <- c(200, 400, 800)  # number of species
n.reps <- c(5, 10, 30, 999) # number of repeated measures per specie (999 refers to unbalanced design)
sigma2.s <- c(0.05, 0.25)    # variance of non-phylogenetic variance component (species level)
sigma2.p <- c(0.05, 0.25)    # variance of phylogenetic variance component (species level)
sigma2.e <- 0.2              # variance of sampling error (estimate level level)


# number of replications per scenario
repl <- 500
sim <- rep(1:repl)

# make table of all scenarios
tab2a <- expand.grid(sim = sim, name = name, n.reps = n.reps, k.species = k.species,
                     sigma2.s = sigma2.s, sigma2.p = sigma2.p, sigma2.e = sigma2.e)

# add columns for to store results and job number
tab2a$save_location <- rep("/srv/scratch/z5394590/pglmm/set2a", each=nrow(tab2a))
tab2a$job_number <- c(1:nrow(tab2a))
conds <- length(k.species) * length(sigma2.s) * length(sigma2.p) * length(sigma2.e)
tab2a$scenario <- rep(1:conds, each = repl)

# save as csv file
write.csv(tab2a, "output/job_array_set2a.csv", row.names = FALSE)



############### No repeated measures
name <- "set2b"

# model parameters
n.reps <- 1 # one measures per specie

# make table of all scenarios
tab2b <- expand.grid(sim = sim, name = name, n.reps = n.reps, k.species = k.species,
                     sigma2.e = sigma2.e, sigma2.p = sigma2.p)

# add columns for to store results and job number
tab2b$save_location <- rep("/srv/scratch/z5394590/", each=nrow(tab2b))
tab2b$job_number <- c(1:nrow(tab2b))
conds <- length(k.species) * length(sigma2.e) * length(sigma2.p) 
tab2b$scenario <- rep(1:conds, each = repl)

# save as csv file
write.csv(tab2b, "output/job_array_set2b.csv", row.names = FALSE)
