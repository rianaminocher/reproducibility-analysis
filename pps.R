library(rethinking)
library(mvtnorm)

png(filename = "output/supp_pps.png",
    res = 250,
    height = 2000,
    width = 2000)

layout(matrix(c(1, 2, 1, 3), nrow = 2, ncol = 2))

par(mar = c(5, 15, 4, 15))

A <- rbeta(1e3, 4, 2)
b <- rlnorm(1e3, -2, 1)

age <- 0:63

p <- sapply(age, function(age) A * exp(-b * age))

plot(age,
     apply(p, 2, mean), 
     type = "l", 
     ylim = c(0, 1),
     col = "deeppink2",
     lwd = 2, 
     xlab = "years since publication", 
     ylab = "p")

for (i in 1:100) {
  lines(age, p[i, ], col = col.alpha("deeppink2", 0.1))
}

shade(apply(p, 2, HPDI), age, col = col.alpha("deeppink2", 0.1))

par(xpd = TRUE)
text(0, 1.1, "a)")
par(xpd = FALSE)

par(mar = c(5, 4, 4, 2))

a_paper <- rnorm(1e3, 0, 2.5)
type_sigma <- rexp(1e3, 1)
b_type <- rnorm(1e3, 0, type_sigma)

rho <- rbeta(1e3, 10, 0.5)
alpha <- rexp(1e3, 10)
delta <- rexp(1e3, 1)

age <- 1:63
n_age <- length(age)
age_id <- 1:length(age)
age_dist <- as.matrix(dist(age), method = "manhattan")

K <- matrix(NA, nrow = 63, ncol = 63)
for (i in 1:63) {
  for (j in 1:63) {
        K[i, j] <- mean(rho * exp(-alpha * age_dist[i, j]^2))
  }
  K[i, i] <- 1
}

k_age <- rmvnorm(1e3, rep(0, n_age), mean(delta) * K)

q <- sapply(age, function(age) inv_logit(a_paper + b_type + k_age[, age]))

plot(age,
     apply(q, 2, mean), 
     type = "l", 
     ylim = c(0, 1),
     col = "deeppink2",
     lwd = 2, 
     xlab = "years since publication", 
     ylab = "p")

shade(apply(q, 2, HPDI), age, col = col.alpha("pink", 0.5))

for (i in 1:100) {
  lines(age, q[i, ], col = col.alpha("deeppink2", 0.1))
}

par(xpd = TRUE)
text(0, 1.1, "b)")
par(xpd = FALSE)

cor <- sapply(age_dist[1, ], function(age) rho * exp(-alpha * age^2))

plot(age_dist[1, ], 
     apply(cor, 2, mean),
     ylim = c(0, 1), 
     type = "l", 
     col = "deeppink2", 
     ylab = "correlation between parameters", 
     xlab = "distance between parameters (years)")

for (i in 1:100) {
  lines(age_dist[1, ], cor[i, ], col = col.alpha("deeppink2", 0.1))
}

shade(apply(cor, 2, HPDI), age_dist[1, ], col = col.alpha("deeppink", 0.2))

par(xpd = TRUE)
text(0, 1.1, "c)")
par(xpd = FALSE)

dev.off()
