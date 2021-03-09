library(rethinking)
library(rstan)
library(xtable)

d <- read.csv("input/anon_database.csv")



      # prep data

n <- nrow(d)

n_sub <- sum(!is.na(d$n_results)) # subsample size
index_sub <- d[!is.na(d$n_results), ]$key # store subsample ids

age_range <- (max(d$year): min(d$year))
age_dist <- as.matrix(dist(age_range, method = "manhattan"))
age_id <- match(d$year, age_range)

n_results <- d[index_sub, ]$n_results
x <- d[index_sub, c("n_results", "data_complete", "analysis_clear", "reproduced")]
x <- t(x)
x <- as.matrix(x)



      # fit model

data <- list(n = n, 
             d = as.integer(d$data_available),
             age = as.integer(2018 - d$year),
             type = as.integer(as.factor(d$type)),
             n_sub = n_sub,
             index_sub = index_sub,
             n_stages = 3, 
             x = x,
             age_id = age_id, 
             n_age = max(age_id),
             age_dist = age_dist)

fit <- stan(file = "model.stan", 
            data = data, 
            chains = 4, 
            cores = 4, 
            iter = 1e4)

# note on model syntax:
# p corresponds to p1 in text
# q1, q2, q3 correspond to p2, p3, p4 in text!

# expect some divergent transitions, but this is not a problem, as the chains converge well.

post <- extract.samples(fit)



      # make plots

png("output/pred_main.png",
    res = 250, 
    height = 3000, 
    width = 4000)

cols <- c("#3B0F70FF", "#782281FF", "#B63679FF", "#ED5A5FFF", "#FE9F6DFF")

layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = FALSE), c(10, 9), c(5, 4), TRUE)

# Panel A: data decay

par(mar = c(5, 5, 4, 0))
par(xpd = FALSE)

age <- 0:63
# compute p across all ages
p <- sapply(age, function(age) apply(post$A, 1, mean) * exp(-(apply(post$b, 1, mean)) * age))

plot(age[1:26], 
     apply(p, 2, mean)[1:26],
     type = "l", 
     lwd = 2,
     xlim = c(0, 25), 
     ylim = c(0, 1), 
     ylab = expression("probability of data recovery p"[1]), 
     xlab = "years since publication", 
     col = cols[5], 
     cex.lab = 2.5, 
     cex.axis = 2)

shade(apply(p, 2, HPDI)[, 1:26], 
      age[1:26], 
      col = col.alpha(cols[5], 0.3))

lines(age[1:26],
      apply(p, 2, HPDI)[1, ][1:26],
      lty = 2, 
      lwd = 2,
      col = cols[5])

lines(age[1:26],
      apply(p, 2, HPDI)[2, ][1:26],
      lty = 2, 
      lwd = 2,
      col = cols[5])

# add observed data

obs <- data.frame(p = tapply(d$data_available, d$year, mean))

# scale points by total observations for each year
sizes <- tapply(d$data_available, d$year, length)
points(2018 - as.numeric(rownames(obs)),
       obs$p, 
       cex = sizes/8)

par(xpd = TRUE)
text("A", x = 0, y = 1.1, cex = 3.5)
par(xpd = FALSE)

# Panel B: 

# marginalising over the observed age distribution

par(mar = c(5, 2, 4, 20))

r <- apply(post$p, 1, mean) * # p across dates observed in our sample
  inv_logit(post$a[, 1]) * # baseline probability
  inv_logit(post$a[, 2]) *
  inv_logit(post$a[, 3])

ps <- list(r, 
           inv_logit(post$a[, 3]), 
           inv_logit(post$a[, 2]), 
           inv_logit(post$a[, 1]),
           post$p)

plot(x = NULL,
     y = NULL, 
     ylim = c(0, 60),
     xlim = c(0, 1), 
     xlab = "p",
     ylab = "", 
     yaxt = "n", 
     cex.lab = 2.5, 
     cex.axis = 2)

offset <- c(0, 15, 25, 35, 45)

for (i in 1:5) {
  dens <- density(ps[[i]])
  dens$y <- (dens$y / max(dens$y)) * 8
  dens$y <- dens$y + offset[i]
  polygon(dens, 
          border = col.alpha(cols[i], 0.8),
          col = col.alpha(cols[i], 0.8),
          lwd = 2)
}

abline(v = lapply(ps, mean), 
       col = cols, 
       lty = 2, 
       lwd = 2)

par(xpd = TRUE)

legend <- c(expression("p"[1]),
            expression("p"[2]),
            expression("p"[3]),
            expression("p"[4]),
            "p(r)")

legend(x = 1.1,
       y = 50,
       legend = legend, 
       cex = 2,
       fill = sapply(cols[c(5, 4, 3, 2, 1)], function(col) col.alpha(col, 0.8)))

text("B", x = 0, y = 70, cex = 3.5)
par(xpd = FALSE)

# Panel C:
par(mar = c(5, 2, 4, 10))

plot(NULL,
     ylim = c(0, 5),
     xlim = c(0, 22),
     main = "",
     ylab = "",
     xlab = "years since publication",
     yaxt = "n",
     cex.lab = 2.5, 
     cex.axis = 2)

for (i in 1:4) {
  points(mean(-log(0.5)/ post$b[, i]), i, 
         col = cols[5], 
         lwd = 2, 
         cex = 2)
}

for (i in 1:4) {
  arrows(y0 = i, 
         y1 = i, 
         x0 = HPDI(-log(0.5)/ post$b[, i])[1],
         x1 = HPDI(-log(0.5)/ post$b[, i])[2],
         angle = 90, 
         length = 0, 
         col = cols[5], 
         lwd = 2) 
}

text(x = c(16, 16, 13, 13),
     y = c(1:4), 
     c("exp", "obs"), 
     cex = 1.7)

text(x = c(20, 18), 
     y = c(1.5, 3.5), 
     c("human", "nonhuman"), 
     cex = 1.7)

text(x = 18, y = 1.5, "]", cex = 4.5, family = 'Helvetica Neue UltraLight')
text(x = 15, y = 3.5, "]", cex = 4.5, family = 'Helvetica Neue UltraLight')

par(xpd = TRUE)
text("C", x = 0, y = 6, cex = 3.5)

dev.off()

# print model results just p's

ps[[5]] <- apply(ps[[5]], 1, mean) # average across samples 
s<- precis(ps[5:1], prob = 0.89)[1:4]
rownames(s) <- c("$p_{1}$", "$p_{2}$", "$p_{3}$", "$p_{4}$", "$p(r)$")

print(xtable(s), 
      file = "output/tab1.txt", 
      only.contents = TRUE, 
      sanitize.rownames.function = function(x) {x})
        
# print model results full

s <- precis(fit, pars = c("A", "b", "a", "paper_sigma", "rho", "alpha", "delta"), 2, prob = 0.89)
rownames(s) <- c("$\\alpha_{1}$", "$\\alpha_{2}$", "$\\alpha_{3}$", "$\\alpha_{4}$",
                 "$\\lambda_{1}$", "$\\lambda_{2}$", "$\\lambda_{3}$", "$\\lambda_{4}$",
                 "$\\phi_{2}$", "$\\phi_{3}$", "$\\phi_{4}$",
                 "$\\psi_{2}$", "$\\psi_{3}$", "$\\psi_{4}$",
                 "$\\rho_{2}$", "$\\rho_{3}$", "$\\rho_{4}$",
                 "$\\upsilon_{2}$", "$\\upsilon_{3}$", "$\\upsilon_{4}$", 
                 "$\\delta_{2}$", "$\\delta_{3}$", "$\\delta_{4}$")
s <- as.data.frame(s)
colnames(s) <- c("mean", "sd", "lower", "upper", "n_eff", "Rhat")

print(xtable(s), 
      file = "output/tabfull.txt", 
      only.contents = TRUE,
      sanitize.rownames.function = function(x) {x})

# decomposing r by age

age <- 0:63
p <- sapply(age, function(age) apply(post$A, 1, mean) * exp(-(apply(post$b, 1, mean)) * age))

q <- array(dim = c(nrow(post$A), 3, 64))
age <- 1:64
for(i in 1:3) {
  q[ , i, ] <- sapply(age, function(age) inv_logit(post$a[, i] 
                                                   + post$k_age[, i, age]))
}

r <- p * q[, 1, ] * q[, 2, ] * q[, 3, ]

apply(r, 2, mean)
apply(r, 2, HPDI)

png("output/supp_decay_age.png", 
    height = 1800, 
    width = 2000, 
    res = 250)

par(mar = c(5, 5, 2, 2))

plot(0:63, 
     apply(r, 2, mean), 
     type = "l", 
     ylim = c(0, 1), 
     xlim = c(0, 63), 
     ylab = "overall probability of reproduction p(r)", 
     xlab = "years since publication", 
     col = "#000004FF", 
     lwd = 2, 
     cex.lab = 2, 
     cex.axis = 1.5)
shade(apply(r, 2, HPDI), 
      0:63, 
      col = col.alpha("#000004FF", alpha = 0.3))
lines(0:63,
      apply(r, 2, HPDI)[1, ],
      lty = 2, 
      col = col.alpha("#000004FF", alpha = 0.7), 
      lwd = 2)
lines(0:63,
      apply(r, 2, HPDI)[2, ],
      lty = 2, 
      col = col.alpha("#000004FF", alpha = 0.7),
      lwd = 2)

dev.off()

# k_age caterpillar plot

png("output/supp_age.png", 
    height = 1500, 
    width = 2000, 
    res = 250)

par(mfrow = c(1, 3))

for(j in 1:3) {
  plot(apply(post$k_age[, j, ], 2, mean),
       1:64, 
       xlim = c(-7, 7), 
       xlab = "", 
       ylab = "age")
  
  for (i in 1:64) {
    arrows(y1 = i, 
           y0 = i,
           x1 = apply(post$k_age[, j, ], 2, PI)[, i][1], 
           x0 = apply(post$k_age[, j, ], 2, PI)[, i][2], 
           angle = 90, 
           length = 0)
  } 
}

dev.off()

# supplementary figure including sample sizes

png("output/supp_size.png", 
    height = 1500, 
    width = 3000, 
    res = 250)

par(mfrow = c(1, 2), 
    mar = c(5, 5, 4, 2))

obs <- data.frame(p = tapply(d$key, d$year, length))

plot(2018 - as.numeric(rownames(obs)), 
     obs$p,
     ylab = "number of papers", 
     xlab = "years since publication", 
     type = "l", 
     lty = 2, 
     lwd = 1.5, 
     cex.lab = 2, 
     cex.axis = 1.5)

par(xpd = TRUE)
text(0, 40, "A", cex = 2)
par(xpd = FALSE)

type_tab <- as.matrix(t(table(d$type, d$year)))
cols <- c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")

plot(2018 - as.numeric(names(type_tab[, 1])),
     type_tab[, 1],
     xlim = c(0, 63),
     ylim = c(0, 25),
     ylab = "number of papers",
     xlab = "years since publication", 
     type = "l", 
     col = cols[1],
     lwd = 1.5,
     cex.lab = 2, 
     cex.axis = 1.5)

for (i in 2:4) {
  lines(2018 - as.numeric(names(type_tab[, i])),
        type_tab[, i], 
        col = cols[i], 
        lwd = 2)
}

legend <- c("hum exp", "hum obs", "non exp", "non obs")
legend("topright",
       fill = cols,
       legend = legend, 
       cex = 1.3)

par(xpd = TRUE)
text(0, 28, "B", cex = 2)
par(xpd = FALSE)

dev.off()

# save traceplots

png("output/trace_A.png", 
    res = 250, 
    height = 1400, 
    width = 1400)

traceplot(fit, pars = "A")

dev.off()

png("output/trace_A.png", 
    res = 250, 
    height = 1400, 
    width = 1400)

traceplot(fit, pars = "A")

dev.off()

png("output/trace_b.png", 
    res = 250, 
    height = 1400, 
    width = 1400)

traceplot(fit, pars = "b")

dev.off()

png("output/trace_p234.png", 
    res = 250, 
    height = 800, 
    width = 1600)

traceplot(fit, pars = "a")

dev.off()

