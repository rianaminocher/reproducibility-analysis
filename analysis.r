      # r script to fit models and produce figures for reproducibility study

unlink("output", recursive = TRUE)
dir.create("output")

# load packages
library(rethinking)
library(rstan)

# add functions
texttab <- function(tab,
                    params,
                    head)
{
  tab <- round(tab, 2)
  tab <- cbind(params, tab, stringsAsFactors = FALSE)
  tab <- rbind(head, tab)
  tab <- apply(tab, 1, function(x) {
    tmp <- paste(x, collapse = " & ")
    paste(tmp, "\\\\", sep = " ")
  })
}

# load data
d1 <- read.csv("input/pubs_anon.csv")
d2 <- read.csv("input/results_anon.csv")

      # numbers for figure 1

sum(d1$downloaded) == 62
sum(d1$emailed & !d1$downloaded) == 473
sum(d1$reply_received & d1$emailed & !d1$downloaded) == 315
sum(d1$data_sent) == 105

gate_sums <- list(d1$downloaded,
                  d1$emailed & !d1$downloaded,
                  d1$reply_received & d1$emailed & !d1$downloaded,
                  d1$data_sent,
                  d1$data_sent | d1$downloaded)

gate_sums <- lapply(gate_sums, function(x) {round(sum(x) / nrow(d1), 2) * 100})


      ### full sample analysis

# A = prob of data at age 0
# b = percentage decline in availability with time

# p = A*e^(-b * age in years)
# log(p) = log(A) + (-b * age in years)

# A = baseline expectation at time = 0
# expect relatively high; 0.3 ~ 1.0

# b is a rate; 0.1 is 10% decline per unit time (year)
# expect data does decay over time to some extent
# like machines with moving parts 

# both A and b must be positive
# A must be bounded between 0 and 1

dens(rbeta(1e3, 4, 2)) # most prob 0.3-1 and peaks at 0.7
dens(rlnorm(1e3, -2, 0.5)) # most prob between 0 and 0.4

      # plot prior

png(filename = "output/pps_data_decay_by_year.png",
    res = 250,
    height = 1400,
    width = 1400)

prior <- data.frame(A = rbeta(1e3, 4, 2),
                    b = rlnorm(1e3, -2, 0.5))

agelist <- 0:50 # sim over age 1 to 50 years

log_p <- sapply(agelist,
                function(t) log_p <- log(prior$A) + (-prior$b * t))

mu_log_p <- apply(log_p, 2, mean)

plot(NULL,
     xlim = c(0, 50),
     ylim = c(0, 1),
     ylab = "probability (material available)",
     xlab = "years")

# plot 200 simulated decay curves

for (i in 1:200) {
  lines(agelist,
        exp(log_p[i, ]),
        col = col.alpha("indianred2", 0.3))
}

lines(agelist, exp(mu_log_p), col = "indianred4", lwd = 2)

dev.off()

# prior assumes decline with prob at 0 being ~0.2 to 1, no strange relationships

      # fit model m1

# stan code
m1_code <-
"data {
  int <lower = 1> N;
  int y[N];
  real x[N];
}

parameters {
  real <lower = 0, upper = 1> A;
  real <lower = 0, upper = 1> b;
}

model {
  vector[N] p;
  A ~ beta(4, 2);
  b ~ lognormal(-2, 0.5);
  for ( i in 1:N ) {
    p[i] = exp(log(A) - b * x[i]);
  }
  y ~ binomial(1, p);
}"

# data list
data <- list(N = nrow(d1),
             y = as.integer(d1$data_available),
             x = as.integer(2018 - d1$year))

# fit model
m1 <- stan(model_code = m1_code,
           data = data,
           chains = 4,
           cores = 4,
           iter = 10000,
           control = list(adapt_delta = 0.99))


# save model results
head <- c("parameter",
          "mean",
          "sd",
          "5.5\\%",
          "94.5\\%",
          "n eff",
          "Rhat")
tabm1 <- precis(m1, 2)
tabm1 <- texttab(tabm1,
                 params = c("$A$", "$b$"),
                 head = head)

tabm1[1] <- paste(tabm1[1], "\\hline")
writeLines(tabm1, "./output/tablem1.txt")

# save model fit object
saveRDS(m1, file = "output/m1.rds")



      # fit model m2 by study type

# p = a * exp(-b[study] * age)
# log(p) <- log(A) - b[study] * age
# vary the relationship b [% decline with time]

d1$study <- paste(d1$species,
                 d1$data_type,
                 sep = " ")

length(unique(d1$study)) == 4 # 4 study types

m2_code <-
"data {
  int <lower = 1> N;
  int <lower = 1> N_study;
  int <lower = 0, upper = 1> y[N];
  int <lower = 1, upper = N_study> study[N];
  int x[N];
}

parameters {
  vector <lower = 0, upper = 1> [N_study] b;
  real <lower = 0, upper = 1> A;
  real <lower = 0, upper = 1> mu;
  real <lower = 0> sigma;
}

model {
  vector[N] p;
  A ~ beta(4, 2);
  mu ~ lognormal(-2, 0.5);
  sigma ~ exponential(1);
  for (j in 1:N_study) {
  b[j] ~ normal(mu, sigma);
  }
  for (i in 1:N) {
    p[i] = exp(log(A) - b[study[i]] * x[i]);
  }
  y ~ binomial(1, p);
}"

data <- list(N = nrow(d1),
             N_study = length(unique(d1$study)),
             y = as.integer(d1$data_available),
             x = as.integer(2018 - d1$year),
             study = as.integer(as.factor(d1$study)))

m2 <- stan(model_code = m2_code,
           data = data,
           chains = 4,
           cores = 4,
           iter = 10000,
           control = list(adapt_delta = 0.99))

tabm2 <- precis(m2, 2)
tabm2 <- texttab(tabm2,
                 params = c("$b$[1]",
                            "$b$[2]",
                            "$b$[3]",
                            "$b$[4]",
                            "$A$",
                            "$b_\\mu$",
                            "$b_\\sigma$"),
                 head = head)
tabm2[1] <- paste(tabm2[1], "\\hline")
writeLines(tabm2, "./output/tablem2.txt")

saveRDS(m2, file = "output/m2.rds")


      # plot data decay

png(filename = "output/data_decay.png",
    res = 250,
    height = 1200,
    width = 2400)

par(mfrow = c(1, 2),
    mar = c(4, 4, 1, 1))

# panel 1

post <- extract.samples(m1)

agelist <- 0:26

cols <- c("hotpink", "lightseagreen", "thistle", "coral")

# note: 1 out of 2 studies in year 27 isn't 0

log_p <- sapply(agelist,
                function(age) {log_p <- log(post$A) - post$b * age})

mu_log_p <- apply(log_p, 2, mean)
HPDI_log_p <- apply(log_p, 2, HPDI)

plot(agelist,
     exp(mu_log_p),
     ylim = c(0, 1),
     xlim = c(0, 26),
     type = "l",
     ylab = "probability (material available)",
     xlab = "years since publication",
     col = cols[1],
     lwd = 1.5,
     cex.axis = 1,
     cex.lab = 1.3)

shade(exp(HPDI_log_p),
      agelist,
      col = col.alpha(cols[1], 0.3))

lines(agelist,
      exp(HPDI_log_p[1, ]),
      lty = 2,
      col = cols[1])

lines(agelist,
      exp(HPDI_log_p[2, ]),
      lty = 2,
      col = cols[1])

years <- sort(unique(d1$year))

freq <- sapply(years, function(x) {
  sum(d1[d1$year == x, ]$data_available) /
    nrow(d1[d1$year == x, ])})
# compute the number of 1s / total sample that year

freqs <- data.frame(years = years,
                    freq = freq)

freqs$years <- 2018 - freqs$years

freqs <- freqs[order(freqs$years), ]

points(freqs$years[1:length(agelist)],
       freqs$freq[1:length(agelist)],
       cex = 0.8)

# add half-life

half_life <- log(0.5) / (-post$b)

sd(half_life)
text(23, 0.9, paste("t(½) =", round(mean(half_life), 2)))

mtext(side = 3, adj = 0, "a)", cex = 1.3)

# panel 2
# 1 - hum obs, 2 - hum exp, 3 - non exp, 4 - non obs

post <- extract.samples(m2)

half_life <- apply(post$b, 2, function(b) -log(0.5) / b)

half_life_mu <- apply(half_life, 2, mean)
half_life_HPDI <- apply(half_life, 2, function(x) HPDI(x, prob = 0.945))

plot(density(half_life[, 1]),
     ylim = c(0, 1),
     xlim = c(0, 28),
     col = col.alpha(cols[1], 0.5),
     main = "",
     ylab = "",
     xlab = "years since publication",
     cex.lab = 1.3)

polygon(density(half_life[, 1]),
     ylim = c(0, 0.8),
     xlim = c(1, 25),
     col = col.alpha(cols[1], 0.3),
     border = col.alpha(cols[1], 1))

polygon(density(half_life[, 2]),
        col = col.alpha(cols[2], 0.3),
        border = col.alpha(cols[2], 1))

polygon(density(half_life[, 3]),
        col = col.alpha(cols[3], 0.3),
        border = col.alpha(cols[3], 1))

polygon(density(half_life[, 4]),
        col = col.alpha(cols[4], 0.3),
        border = col.alpha(cols[4], 1))

xleft <- c(20, 22, 20, 22)
ybottom <- c(0.80, 0.80, 0.75, 0.75)
xright <- c(22, 24, 22, 24)
ytop <- c(0.85, 0.85, 0.80, 0.80)

rect(xleft, ybottom, xright, ytop,
     col = sapply(cols[c(1, 3, 2, 4)],
     function(cols) col.alpha(cols, 0.5)))

text(20.5, 0.88, labels = "hum")
text(23.5, 0.877, labels = "non")
text(18.5, 0.82, labels = "exp")
text(18.5, 0.77, labels = "obs")
text(13.5, 0.30, paste("t(½) =", round(mean(half_life[, 1]), 2)))
text(9, 0.85, paste("t(½) =", round(mean(half_life[, 3]), 2)))

mtext(side = 3, adj = 0, "b)", cex = 1.3)

dev.off()



      # plot study type freq vs year

# (a) frq. of human vs nonhuman over time (1 line)
# (b) frq. of experimental vs observational (1 line)
# (c) frq. of each of the 4 types (4 lines)

png(filename = "output/typebyyear_supp.png",
    res = 250,
    height = 1000,
    width = 3000)

par(mfrow = c(1, 3), mar = c(4.5, 4.5, 1, 0.8))

cols <- c("hotpink", "lightseagreen", "thistle", "coral")

typenames <- c("species", "data_type")
legendnames <- c("hum", "exp")

years <- sort(unique(d1$year))

freq <- table(d1[, typenames[[1]]], d1$year)
freq <- freq[1, ] / (freq[1, ] + freq[2, ])
plot(years,
     freq,
     col = cols[[1]],
     type = "b",
     lty = 1,
     lwd = 2,
     ylab = "frequency",
     xlab = "year",
     cex.lab = 1.5)

legend("topright",
       legendnames[[1]],
       fill = cols[[1]],
       cex = 1.2)

freq <- table(d1[, typenames[[2]]], d1$year)
freq <- freq[1, ] / (freq[1, ] + freq[2, ])
plot(years,
     freq,
     col = cols[[2]],
     type = "b",
     lty = 1,
     lwd = 2,
     ylab = "",
     xlab = "year",
     cex.lab = 1.5)

legend("topright",
       legendnames[[2]],
       fill = cols[[2]],
       cex = 1.2)

studies <- print(unique(d1$study))

freq <- table(d1$year, d1$study)
freqout <- matrix(nrow = nrow(freq), ncol = ncol(freq))
for (i in 1:4){
  freqout[, i] <- freq[, i] / apply(freq, 1, sum)
}

plot(years,
     freqout[, 1],
     ylab = "",
     xlab = "year",
     col = cols[[1]],
     type = "b",
     lty = 1,
     lwd = 2,
     cex.lab = 1.5)

for (i in 2:4) {
  points(years,
       freqout[, i],
       col = cols[[i]],
       type = "b",
       lty = 1,
       lwd = 2)
}

legend("topright",
       legend = c("hum obs", "hum exp", "non exp", "non obs"),
       fill = cols,
       cex = 1.2)

dev.off()



      ### subsample analysis

# when did we have code?
length(which(tapply(d2$scripted, d2$key, sum) >= 1))

# logistic regression, n results per paper
# cluster by paper id

# constrain prior probability to be equal from 0 and 1



      # fit model m3

# 4 outcomes: reproduce, reanalysis, analysis, data

d2$data_understandable <- d2$reproduced | d2$failure_category %in% c("incomplete data", "analysis unclear", "reanalysis result differs")

d2$data_complete <- d2$reproduced | d2$failure_category %in% c("analysis unclear", "reanalysis result differs")

# create a composite outcome
d2$data_understandable_complete <- d2$data_understandable & d2$data_complete

d2$analysis_understandable <- d2$reproduced | d2$failure_category %in% c("reanalysis result differs")

outcomes <- c("data_understandable_complete",
              "analysis_understandable",
              "reproduced")

# use a uniform prior here instead, to show updating

m3_code <-
"data {
  int <lower = 1> N;
  int <lower = 1> N_papers;
  int <lower = 0, upper = 1> y[N];
  int <lower = 1> paper[N];
}

parameters {
  vector[N_papers] a;
  real <lower = 0, upper = 1> mu;
  real <lower = 0> sigma;
}

model {
  mu ~ normal(0, 1.5);
  sigma ~ exponential(1);
  for (j in 1:N_papers){
    a[j] ~ normal(mu, sigma);
  }
  for (i in 1:N) {
    y[i] ~ bernoulli_logit(a[paper[i]]);
  }
}"

data_list <- list(length = length(outcomes))

for (i in 1:length(outcomes)) {
  data_list[[i]] <- list(N = nrow(d2),
                         N_papers = length(unique(d2$key)),
                         y = as.integer(d2[, outcomes[[i]]]),
                         paper = as.integer(as.factor(d2$key)))
}

# fit main model

m3 <- stan(model_code = m3_code,
          data = data_list[[3]],
          chains = 4,
          cores = 4,
          iter = 10000,
          control = list(adapt_delta = 0.99))

saveRDS(m3, file = "output/m3.rds")


      # plot posterior/prior "reproduced"

cols <- c("paleturquoise3", "hotpink2")

png(filename = "output/results_repro.png",
    res = 250,
    height = 1400,
    width = 1400)

par(mfrow = c(1, 1))
prior_mu <- rnorm(1e4, 0, 1.5) # use uniform prior to show updating

plot(density(inv_logit(prior_mu)),
     ylim = c(0, 7),
     xlim = c(0, 1),
     lty = 2,
     lwd = 2,
     col = col.alpha(cols[1], 1),
     main = "",
     ylab = "",
     xlab = "probability of successful reproduction",
     cex.lab = 1.3)

polygon(density(inv_logit(prior_mu)),
        ylim = c(0, 0.8),
        xlim = c(1, 30),
        col = col.alpha(cols[1], 0.3),
        border = NA)

post <- extract.samples(m3)

lines(density(inv_logit(post$mu)),
      lwd = 2,
      col = col.alpha(cols[2], 1))

polygon(density(inv_logit(post$mu)),
        col = col.alpha(cols[2], 0.2),
        border = NA)

abline(v = inv_logit(mean(post$mu)),
       lty = 2,
       lwd = 2)

dev.off()



# fit model to all outcomes

m3_allfits <- lapply(data_list, function(data) {
  stan(model_code = m3_code,
       data = data,
       chains = 4,
       cores = 4,
       iter = 10000,
       control = list(adapt_delta = 0.99))
})

samples <- lapply(m3_allfits, extract.samples)

mu <-
  lapply(samples, function(x) {
    c(round(inv_logit(mean(x$mu)), 3),
      round(inv_logit(HPDI(x$mu, prob = 0.945)), 3))
  })

# save model results for all

for (i in 1:length(m3_allfits)) {
  tab <- precis(m3_allfits[[i]], 2)
  tab <- texttab(tab,
                 params = paste(rownames(precis(m3, 2))),
                 head = c("parameter",
                          "mean",
                          "sd",
                          "5.5\\%",
                          "94.5\\%",
                          "n eff",
                          "Rhat"))
  tab[1] <- paste(tab[1], "\\hline")
  tab <- gsub("mu", "$\\mu$", tab)
  tab <- gsub("sigma", "$\\sigma$", tab)
  writeLines(tab, paste0("./output/tablemsub", i, ".txt"))
}



      # model m4 influence of publication age

# with a gaussian process on the beta vector for year

# we could do 2018:1991, to estimate years
# that are in the full sample

years <- 2018:1991
# make distance matrix
yeard <- as.matrix(dist(years, method = "manhattan"))
year_id <- match(d2$year, years)

m4_code <-
"functions {
  matrix cov_GPL2(matrix x, real eta, real rho, real delta) {
  int N = dims(x)[1];
  matrix[N, N] K;
  for (i in 1:(N-1)) {
     K[i, i] = eta + delta;
     for (j in (i + 1):N) {
          K[i, j] = eta * exp(-rho * square(x[i,j]) );
          K[j, i] = K[i, j];
        }
      }
  K[N, N] = eta + delta;
  return K;
  }
}

data {
  int N_papers;
  int N_years;
  int N;
  int y[N];
  int year[N];
  int paper[N];
  matrix[N_years, N_years] yeard;
}

parameters {
  vector[N_papers] a;
  vector[N_years] beta;
  real mu;
  real<lower=0> sigma;
  real<lower=0> eta;
  real<lower=0> rho;
}

model {
  vector[N] p;
  matrix[N_years, N_years] SIGMA;
  rho ~ exponential(0.5);
  eta ~ exponential(2);
  SIGMA = cov_GPL2(yeard, eta, rho, 0.01);
  beta ~ multi_normal(rep_vector(0, N_years), SIGMA);
  mu ~ normal(0, 1.5);
  sigma ~ exponential(1);
  for (j in 1:N_papers){
    a[j] ~ normal(mu, sigma);
  }
  for (i in 1:N) {
    y[i] ~ bernoulli_logit(a[paper[i]] + beta[year[i]]);
  }
}"

data <- list(N = nrow(d2),
             N_papers = length(unique(d2$key)),
             N_years = length(unique(years)),
             y = as.integer(d2$reproduce),
             paper = as.integer(as.factor(d2$key)),
             year = year_id,
             yeard = yeard)

m4 <- stan(model_code = m4_code,
           data = data,
           iter = 10000,
           chains = 4,
           cores = 4,
           control = list(adapt_delta = 0.99))



      # plot results for supp

png(filename = "output/analytic_by_year.png",
    res = 250,
    height = 1400,
    width = 1400)

post <- extract.samples(m4)
beta_est <- apply(post$beta, 2, mean)
beta_ci <- apply(post$beta, 2, HPDI)

plot(years,
     logistic(beta_est),
     ylim = c(0, 1),
     type = "l",
     ylab = "probability of analytical reproducibility",
     xlab = "age of publication")

shade(logistic(beta_ci), years)

lines(years,
      logistic(beta_ci[1, ]),
      type = "l",
      lty = 2)

points(years,
       logistic(beta_ci[2, ]),
       type = "l",
       lty = 2)

dev.off()

# we don't see an effect of year
# and we don't have much confidence in the estimate
# its because of this we can combine P(D) and P(R|D) w/o conditioning on age


      ### combined probability calculation

# p1 = p(D)
# p2 = p(R|D)

# P(D & R) = p(D) * p(R|D) + p(!D) * p(R|!D)
# p(!D) * p(R|!D) is assumed to be zero!

post <- extract.samples(m3)

p2 <- c(inv_logit(mean(post$mu)),
        inv_logit(HPDI(post$mu, prob = 0.945)))

# for p1 we will fit a simple no predictor logistic model
# use a flat prior, like for p2 normal(0, 1.5)

m5_code <-
"data{
  int N;
  int y[N];
}

parameters{
  real mu;
}

model{
  vector[N] p;
  mu ~ normal(0, 1.5);
  for (i in 1:N) {
  p[i] = mu;
  }
  y ~ bernoulli_logit(p);
}"

data <- list(N = nrow(d1),
             y = as.integer(d1$data_available))

m5 <- stan(model_code = m5_code,
           data = data,
           chains = 4,
           cores = 4,
           iter = 10000,
           control = list(adapt_delta = 0.99))

post <- extract.samples(m5)

p1 <- c(inv_logit(mean(post$mu)),
        inv_logit(HPDI(post$mu, prob = 0.945)))

combined_prob <- p1 * p2
