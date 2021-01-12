library(rethinking)

unlink("output", recursive = TRUE)
dir.create("output")

# our analysis estimates the average probability of success across all papers

# for a set of n_papers:
# simulate a mean probability across all papers and results (parameter a) --> this is what we estimate
# simulate a number of results per paper
# simulate a number of successful reproductions per paper, with probability (logistic(a))
# fit the stats model and demonstrate at n = 40, we can recover a 

simulate_data <- function(n_papers, a) {
  
  # transform a to prob
  p <- inv_logit(a) 
  
  # simulate number of results per paper
  # poisson + 1 so none are 0
  n_results <- as.integer(rpois(n_papers, 2) + 1)
  
  n_reproduced <- rbinom(n_papers, 
                         n_results, 
                         p)
  
  pars <- list(
    n_papers = n_papers,
    p = p)
  
  data <- list(
    n_papers = n_papers,
    n_reproduced = n_reproduced,
    n_results = n_results
  )
  
  out <- list(pars = pars,
              data = data)
  return(out)
  
}

model <- "
data {
  int <lower = 1> n_papers;
  int <lower = 1> n_results [n_papers];
  int <lower = 0, upper = 10> n_reproduced [n_papers];
}

parameters{

  real a; 
  real <lower = 0> paper_sigma;
  vector [n_papers] z_a_paper;

}

transformed parameters{
  vector [n_papers] a_paper;
  a_paper = z_a_paper * paper_sigma;
}

model{
  vector [n_papers] p;

  a ~ normal(0, 2.5);

  paper_sigma ~ exponential(1);
  
  z_a_paper ~ normal(0, 1);

    for (i in 1:n_papers) {   
      p[i] = inv_logit(a + a_paper[i]);
    }
    n_reproduced ~ binomial(n_results, p);
}
"

# fix a at 0, ~ prob of 0.5
# can assign a distribution too ofc

a <- 0

# set the seed so the figure dimensions are the same when reproducing
set.seed(19886)

      # show updating from prior at n = 40

png("output/fig_n40_response.png", 
    height = 1400, 
    width = 1400,
    res = 250)

# plot prior

par(mar = c(5, 5, 2, 2))

plot(NULL, 
     xlab = "probability", 
     ylab = "density", 
     ylim = c(0, 15), 
     xlim = c(0, 1), 
     cex.lab = 1.5)

dens(inv_logit(rnorm(1e3, 0, 2.5)), 
     lty = 2, 
     lwd = 2, 
     add = TRUE)

# plot true a
abline(v = inv_logit(a), 
       lwd = 2, 
       col = "red")

# plot posterior at n = 40, n = 5, n = 100

n_list <- c(5, 40, 100)

for (i in 1:3) {
  sim <- simulate_data(n_papers = n_list[i], a = 0)
  data <- sim$data
  
  fit <- stan(model_code = model, 
              data = data, 
              iter = 1e5,
              chains = 4, 
              cores = 4, 
              seed = 19886)
  
  post <- extract.samples(fit)
  
  dens(inv_logit(post$a), 
       add = TRUE,
       lwd = i+1, 
       xlab = "probability")
} 

axis(3, labels = FALSE, tck = 0.02)
axis(4, labels = FALSE, tck = 0.02)

dev.off()


      # set of "n"s to loop over

n_list <- 10:200

dat_list <- lapply(n_list, function(x) simulate_data(a = 0, n_papers = x))
mod_list <- lapply(dat_list, function(data) { 
                     stan(model_code = model, 
                          data = data$dat, 
                          control = list(adapt_delta = 0.99), 
                          iter = 10000,
                          chains = 4, 
                          cores = 4)
                   }
)


compute_width <- function(mod) {
  
  post <- extract.samples(mod)
  sd <- sd(inv_logit(post$a))
  return(sd)
  
}

sdlist <- lapply(mod_list, compute_width)

png("./output/n_40_HPDI.png", 
    res = 250, 
    height = 1400, 
    width = 1400)

par(mar = c(5, 5, 2, 2))

plot(n_list, 
     sapply(sdlist, diff), 
     xlab = "sample size", 
     ylab = "width of HPDI",
     tck = 0.02, 
     cex.lab = 1.5)

abline(v = 40, lty = 2, lwd = 2, col = "red")
axis(3, labels = FALSE, tck = 0.02)
axis(4, labels = FALSE, tck = 0.02)

dev.off()


