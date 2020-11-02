# this script justifies our selection of n = 40 studies to sample

# we first simulate the probability of a result being correct
# from a hypergeomtric distribution

# next, we simulate a population probability 
# based on n results per paper
# from a binomial distribution

# we plot 
# 1) the prior and posterior distribution
# for a p simulated from a hypergeometric at n = 40
# 2) the sd of the posterior distribution at various values of n from 1 to 100, 
# for ps simulated from a hypergeometric
# 3) the prior and posterior distribution 
# for a p simulated from a binomial at n = 40
# 4) the sd of the posterior distribution at various values of n from 1 to 100
# for a p simulated from a binomial

library(rethinking)
library(rstan)

# make a folder to store plots
unlink("sample-output", recursive = TRUE)
dir.create("sample-output")

# 1

      ### hypergeometric sampling distribution ###
      # for 1 "code" per paper

      # make a more intuitive function for dhyper/rhyper than R's:

dhyper2 <- function(x, N, p, k, log = FALSE){
  m <- N * p
  n <- N * (1 - p)
  dhyper(x, m = m, n = n, k = k, log = log)
}

rhyper2 <- function(nn, N, p, k) {
  m <- N * p
  n <- N * (1 - p)
  rhyper(nn, m = m, n = n, k = k)
}

# verify it works:

m <- 20
n <- 10
p <- m / (m + n)
N <- m + n

k <- 5

dhyper(0:5, m = m, n = n, k = k)
sum(dhyper(0:5, m = m, n = n, k = k))
dhyper2(0:5, N = N, p = p, k = k)
sum(dhyper2(0:5, N = N, p = p, k = k))

identical(dhyper2(0:5, N = N, p = p, k = k), 
          dhyper(0:5, m = m, n = n, k = k))



      # determine appropriate sample size
      # for the subset from ~150 total papers with materials

# a single example

# total pop size
N <- 150 

# assign a sample size
# demonstrate with chosen n 40
sample_size <- 40

# values of p to iterate over
p_grid <- seq(0, 1, length.out = (N + 1))

# simulate our "true" p from the beta
sim_p <- rbeta2(1, prob = 0.8, theta = 10)

# simulate a total observation from the hypergeom
observation <- rhyper2(1,
                       N = N, 
                       p = sim_p,
                       k = sample_size)

# assign the prior
prior <- dbeta2(p_grid, prob = 0.8, theta = 10)
prior <- prior/sum(prior)

likelihoods <- dhyper2(x = observation, 
                       N = N, 
                       p = p_grid,
                       k = sample_size)

# compute the posterior
posterior <- likelihoods * prior / sum(likelihoods * prior)

png("./sample-output/hypergeom-n40.png", 
    res = 250,
    height = 1400, 
    width = 1400)

plot(p_grid, 
     prior,
     type = "l",
     ylim = c(0, 0.1), 
     lty = 2, 
     xlab = "posterior probability (solid) shrinks substantially from prior (dotted)", 
     ylab = "")

points(p_grid, 
       posterior, 
       type = "l")

dev.off()

# we see updating from prior with size 40 
# no matter whether our true p is within our prior or not

# estimate the posterior mean
sum(p_grid * posterior)

# estimate the posterior standard deviation

pri <- sample(p_grid, 10000, replace = TRUE, prob = prior)
post <- sample(p_grid, 10000, replace = TRUE, prob = posterior)

mean(pri < 0.5)
mean(post < 0.5)
sd(post)
HPDI(post)



      # wrap into a function
      # plot SD at different sample sizes

sim_sample_hyper <- function(N = 150,
                             sample_size, 
                             p = 0.7, 
                             theta = 10) 
{
  
  p_grid <- seq(0, 1, length.out = (N + 1))
  
  sim_p <- rbeta2(1, prob = p, theta = theta)
  observation <- rhyper2(1, N = N, p = sim_p, k = sample_size)
  
  prior <- dbeta2(p_grid, prob = p, theta = theta)
  prior <- prior/sum(prior)
  
  likelihoods <- dhyper2(x = observation,
                         N = N,
                         p = p_grid, 
                         k = sample_size)
  
  posterior <- likelihoods * prior / sum(likelihoods * prior)
  
  sum(p_grid * posterior)
  
  pri <- sample(p_grid, 10000, replace = TRUE, prob = prior)
  post <- sample(p_grid, 10000, replace = TRUE, prob = posterior)
  
  p_est <- mean(post)
  p_sd <- sd(post)
  p_HPDI <- HPDI(post)
  
  output <- list(
    p = p,
    theta = theta,
    sim_p = sim_p,
    sample_size = sample_size,
    observation = observation,
    p_est = p_est,
    p_sd = p_sd,
    p_HPDI = p_HPDI
  )
  
  return(output)
  
}


# seq of sample sizes to simulate over
x <- rep(seq(10, 100, by = 1), 10)

results <- mapply(sim_sample_hyper, 
                  sample_size = x, 
                  SIMPLIFY = "array")

str(results)

sds <- NA
for(i in 1:length(x)) sds[i] <- results[, i]$p_sd

png("./sample-output/hypergeom-sim.png", 
    res = 250,
    height = 1400, 
    width = 1400)

plot(x, 
     sds,
     ylim = c(0, 0.1),
     xlim = c(0, 100),
     ylab = "SD on posterior probability", 
     xlab = "sample size (out of total n ~ 150)",
     tck = 0.02, 
     cex = 0.5)
abline(v = 40, lty = 2, col = "red")
axis(3, labels = FALSE, tck = 0.02)
axis(4, labels = FALSE, tck = 0.02)

dev.off()

# at n = 40, we get an sd of about ~ 0.04 - 0.06


# 2

      ### sampling from the binomial ###
      # for multiple results per paper

simulate_data <- function(n_papers) 
  
  {
  
    a <- rnorm(1, 0.7, 0.2)
    # prior from hypergeometric model, in log odds space
  
    # simulate a
    logit_p <- a 
    
    # simulate number of results per paper
    # poisson + 1 so none are 0
    n_results <- as.integer(rpois(n_papers, 2) + 1)
  
    n_reproduced <- rbinom(n_papers, 
                           n_results, 
                           logistic(logit_p))
  
    pars <- list(
      n_papers = n_papers,
      a = a
    )
  
    dat <- list(
      n_reproduced = n_reproduced,
      n_results = n_results
    )
  
    out <- list(
      pars = pars,
      dat = dat
    )
  
  return(out)
  
}


# fit stats model
m <- alist(
  n_reproduced ~ dbinom(n_results, p),
  logit(p) <- a,
  a ~ dnorm(0, 0.5)
  )


# check at sample size = 40
sim <- simulate_data(40)

true_value <- inv_logit(sim$pars$a)

m1 <- ulam(m,
           data = sim$dat, 
           control = list(adapt_delta = 0.99), 
           iter = 5000)

post1 <- extract.samples(m1)

mean(inv_logit(rnorm(5000, 0, 0.5))) # prior mean
HPDI(inv_logit(rnorm(5000, 0, 0.5))) # prior sd

mean(inv_logit(post1$a))
HPDI(inv_logit(post1$a))

# n = 40 shows updating from prior

png("./sample-output/binom-n40.png", 
    res = 250, 
    height = 1400, 
    width = 1400)

dens(logistic(rnorm(1000, 0, 0.5)),
     ylim = c(0, 15), 
     lty = 2, 
     ylab = "",
     xlab = "posterior probability (solid) shrinks from prior (dotted)")
dens(logistic(post1$a), add = TRUE)
text(x = 0.2, 
     y = 12, 
     labels = paste0("simulated p = ",
                     round(true_value, 2), 
                     "\nat n = 40") 
     )

dev.off()

# simulate over a bunch of n's
# for a given n_papers, the dataset needs to be created many times
# each time fit the model
# and report the results in aggregate

plist <- seq(from = 10, to = 100, by = 1)
plist <- rep(plist, each = 10) 

dat_list <- lapply(plist, simulate_data)

mod_list <- lapply(dat_list, 
                      function(data) { 
                      ulam(m, 
                           data = data$dat, 
                           control = list(adapt_delta = 0.99), 
                           iter = 10000,
                           chains = 4, 
                           cores = 4)
                      }
                   )

str(mod_list)

# plot sds

compute_width <- function(mod) {
  
  post <- extract.samples(mod)
  sd <- sd(inv_logit(post$a))
  return(sd)
  
}

sdlist <- lapply(mod_list, compute_width)

png("./sample-output/binom-sim.png", 
    res = 250, 
    height = 1400, 
    width = 1400)

plot(plist, 
     unlist(sdlist), 
     xlab = "sample size (out of total n ~ 150)", 
     ylab = "SD on posterior probability of p",
     tck = 0.02, 
     cex = 0.5)
abline(v = 40, lty = 2, col = "red")
axis(3, labels = FALSE, tck = 0.02)
axis(4, labels = FALSE, tck = 0.02)

dev.off()
