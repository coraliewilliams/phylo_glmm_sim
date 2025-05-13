#--- Testing models on default bayesian default tuning settings




#### Test 1
# brms default iter; MCMCglmm default nitt
# load test 1 results
load("~/Projects/phylo_glmm_sim/output/tests/sim_results_test1.RDATA")

# distribution of ESS looks very skewed to the right
hist(dat$ESS)

# get number of iterations
iters <- length(dat$model)/5

# get percentage of ESS < 400 for brms
table(dat$ESS[which(dat$model=="brms")]>400)/iters*100
# get percentage of ESS < 400 for MCMCglmm
table(dat$ESS[which(dat$model=="MCMCglmm")]>400)/iters*100

# check run time across all models 
boxplot(dat$run_time ~ dat$model)



#### Test 2
# brms default iter x2; MCMCglmm default nitt x5
# load test 2 results
load("~/Projects/phylo_glmm_sim/output/tests/sim_results_test2.RDATA")

# distribution of ESS looks very skewed to the right
hist(dat$ESS)

# get number of iterations
iters <- length(dat$model)/5

# get percentage of ESS < 400 for brms
table(dat$ESS[which(dat$model=="brms")]>400)/iters*100
# get percentage of ESS < 400 for MCMCglmm
table(dat$ESS[which(dat$model=="MCMCglmm")]>400)/iters*100

# check run time across all models 
boxplot(dat$run_time ~ dat$model)





#### Test 3
# brms default iter x2; MCMCglmm default nitt x30
# load test 3 results
load("~/Projects/phylo_glmm_sim/output/tests/sim_results_test3.RDATA")

# distribution of ESS looks very skewed to the right
hist(dat$ESS)

# get number of iterations
iters <- length(dat$model)/5

# get percentage of ESS < 400 for brms
table(dat$ESS[which(dat$model=="brms")]>400)/iters*100
hist(dat$ESS[which(dat$model=="brms")])
# get percentage of ESS < 400 for MCMCglmm
table(dat$ESS[which(dat$model=="MCMCglmm")]>400)/iters*100
hist(dat$ESS[which(dat$model=="MCMCglmm")])


# check run time across all models 
boxplot(dat$run_time ~ dat$model)




#### Test 4
# brms default iter x10; MCMCglmm default nitt x30
load("~/Projects/phylo_glmm_sim/output/tests/sim_results_test4.RDATA")


#### on-goong batch1 testing
load("~/Projects/phylo_glmm_sim/output/tests/sim_results_batch1.RDATA")



# distribution of ESS looks very skewed to the right
hist(dat$ESS)

# get number of iterations
iters <- length(dat$model)/5

# get percentage of ESS < 400 for brms
table(dat$ESS[which(dat$model=="brms")]>400)/iters*100
table(dat$ESS[which(dat$model=="brms"& dat$species_size==100)]>400)
table(dat$convergence[which(dat$model=="brms")])/iters*100
hist(dat$ESS[which(dat$model=="brms")])

# get percentage of ESS < 400 for MCMCglmm
table(dat$ESS[which(dat$model=="MCMCglmm")]>400)/iters*100
table(dat$convergence[which(dat$model=="MCMCglmm")])/iters*100
hist(dat$ESS[which(dat$model=="MCMCglmm")])


# check run time across all models 
boxplot(dat$run_time ~ dat$model)





###################
# Small tests (previous to the above test 3) to get optimum MCMCglmm settings i.e. to reach at least 80% of models with minimum ESS of 400

## 3a: nitt=103000 (x10); thin=100 (x10)
load("~/Projects/phylo_glmm_sim/output/tests/sim_results_test_3a.RDATA")
iters <- length(dat$model)/5 # based on n=2821 sims
table(dat$ESS[which(dat$model=="brms")]>400)/iters*100
table(dat$ESS[which(dat$model=="MCMCglmm")]>400)/iters*100


## 3b: nitt=503000 (x50); thin=200 (x20)
load("~/Projects/phylo_glmm_sim/output/tests/sim_results_test_3b.RDATA")
iters <- length(dat$model)/5 # based on n=82 sims
table(dat$ESS[which(dat$model=="brms")]>400)/iters*100
table(dat$ESS[which(dat$model=="MCMCglmm")]>400)/iters*100


## 3b: nitt=203000 (x20); thin=150 (x15)
load("~/Projects/phylo_glmm_sim/output/tests/sim_results_test_3c.RDATA")
iters <- length(dat$model)/5 # based on n=145 sims
table(dat$ESS[which(dat$model=="brms")]>400)/iters*100
table(dat$ESS[which(dat$model=="MCMCglmm")]>400)/iters*100
boxplot(dat$run_time ~ dat$model)


## 3d: nitt=403000 (x40); thin=10 (default)
load("~/Projects/phylo_glmm_sim/output/tests/sim_results_test_3c.RDATA")
iters <- length(dat$model)/5 # based on n=145 sims
table(dat$ESS[which(dat$model=="brms")]>400)/iters*100
table(dat$ESS[which(dat$model=="MCMCglmm")]>400)/iters*100


## 4a: brms default iter x5; MCMCglmm default nitt x30
load("~/Projects/phylo_glmm_sim/output/tests/sim_results_test_4a.RDATA")
iters <- length(dat$model)/5 # based on n=9777 sims
table(dat$ESS[which(dat$model=="brms")]>400)/iters*100
hist(dat$ESS[which(dat$model=="brms")]) ## brms hsitrogram still seems to have high number of low ESS
table(dat$ESS[which(dat$model=="MCMCglmm")]>400)/iters*100



####################