# k.species - number of species
k.species <- c(10^8)

# simulate number of observations per species
k.obs <- round(rbeta(n = k.species, shape1 = 1.5, shape2 = 3) * 99) + 1
k.obs
hist(k.obs)
min(k.obs)
max(k.obs)
mean(k.obs)
as.numeric(names(sort(-table(k.obs)))[1])
median(k.obs)


##### covariate x
x <- runif(10^6, 10, 20) 
hist(x)
median(x)
