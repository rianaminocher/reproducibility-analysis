# perform robustness checks

# 1: estimate p1 with a logistic, rather than an exponential decay model
# 2: estimate overall reproducibility with a 2-stage model

library(rethinking)
library(xtable)

d <- read.csv("../input/anon_database.csv")



# 1: p1 with a logistic

# the model code is identical to the main model code: model.stan
# except p1 is not estimated as A*exp(-b*t)
# simply a logistic regression
# logit(p1) = a + bt
# a and b both vary with type of study included

# prep data

n <- nrow(d)

n_sub <- sum(!is.na(d$n_results)) # subsample size
index_sub <- d[!is.na(d$n_results), ]$key # store subsample ids

age_range <- (max(d$year): min(d$year))
age_dist <- as.matrix(dist(age_range, method = "manhattan"))
age_id <- match(d$year, age_range)

n_results <- d[index_sub, ]$n_results
x <- d[index_sub, c("n_results", 
                    "data_complete", 
                    "analysis_clear",
                    "reproduced")]
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

fit1 <- stan(file = "robustness-checks-m1.stan", 
             data = data, 
             chains = 4, 
             cores = 4, 
             iter = 1e4)

precis(fit1, 3, pars = c("A", "b", "a_paper"))

post <- extract.samples(fit1)

# plot predicted prob by year

p <- sapply(0:63, function(x) inv_logit(apply(post$A, 1, mean) + apply(post$b, 1, mean) * x))

png("output/logistic_decay.png", 
    res = 250, 
    height = 1000, 
    width = 1200)

plot(0:63, 
     apply(p, 2, mean),
     type = "l", 
     lwd = 2,
     xlim = c(0, 63), 
     ylim = c(0, 1), 
     ylab = expression("probability of data recovery p"[1]), 
     xlab = "years since publication", 
     col = "blue", 
     cex.lab = 1, 
     cex.axis = 1)

shade(apply(p, 2, HPDI), 
      0:63, 
      col = col.alpha("blue", 0.3))

lines(0:63,
      apply(p, 2, HPDI)[1, ],
      lty = 2, 
      lwd = 2,
      col = "blue")

lines(0:63,
      apply(p, 2, HPDI)[2, ],
      lty = 2, 
      lwd = 2,
      col = "blue")

dev.off()


# plot Figure 2B
# for each type, plot the density of the beta value

p <- array(dim = c(4, 2e4, 64))

for (i in 1:4) {
  p[i, , ] <- sapply(0:63, function(x) inv_logit(post$A[, i] + post$b[, i] * x))
}

png("output/logistic_type.png", 
    res = 250, 
    height = 1000, 
    width = 1200)

cols <- c("#2171B5", "#EFF3FF",  "#BDD7E7", "#6BAED6")

plot(0:63, 
     apply(p[1, , ], 2, mean),
     type = "l", 
     lwd = 2,
     xlim = c(0, 63), 
     ylim = c(0, 1), 
     ylab = expression("probability of data recovery p"[1]), 
     xlab = "years since publication", 
     col = cols[1], 
     cex.lab = 1, 
     cex.axis = 1)

for (i in 2:4){
  lines(0:63, 
        apply(p[i, , ], 2, mean), 
        col = cols[i])
}

for (i in 1:4) {
  shade(apply(p[i, , ], 2, HPDI), 
        0:63, 
        col = col.alpha(cols[i], 0.3))
}

for (i in 1:4) {
  lines(0:63,
        apply(p[i, , ], 2, HPDI)[1, ],
        lty = 2, 
        lwd = 2,
        col = cols[i])
}

for (i in 1:4) {
  lines(0:63,
        apply(p[i, , ], 2, HPDI)[2, ],
        lty = 2, 
        lwd = 2,
        col = cols[i])
}

legend("topright", 
       legend = c("hum exp", "hum obs", "non exp", "non obs"), 
       fill = cols)

dev.off()

b_tab <- precis(fit1, 2, pars = "b")
print(xtable(b_tab), 
      file = "output/logistic_bs.txt", 
      only.contents = TRUE, 
      sanitize.rownames.function = function(x) {x})

# calculate the estimated p's

r <- apply(inv_logit(post$p), 1, mean) * # p across dates observed in our sample
  inv_logit(post$a[, 1]) * # baseline probability
  inv_logit(post$a[, 2]) *
  inv_logit(post$a[, 3])

ps <- list(r, 
           inv_logit(post$a[, 3]), 
           inv_logit(post$a[, 2]), 
           inv_logit(post$a[, 1]),
           apply(inv_logit(post$p), 1, mean))

s <- precis(ps[5:1], prob = 0.89)[1:4]
rownames(s) <- c("$p_{1}$", "$p_{2}$", "$p_{3}$", "$p_{4}$", "$p(r)$")

print(xtable(s), 
      file = "output/logistic_ps.txt", 
      only.contents = TRUE, 
      sanitize.rownames.function = function(x) {x})

# save results

s <- precis(fit1, pars = c("A", "b", "a", "paper_sigma"), 2, prob = 0.89)
rownames(s) <- c("$\\alpha_{1}$", "$\\alpha_{2}$", "$\\alpha_{3}$", "$\\alpha_{4}$",
                 "$\\lambda_{1}$", "$\\lambda_{2}$", "$\\lambda_{3}$", "$\\lambda_{4}$",
                 "$\\phi_{2}$", "$\\phi_{3}$", "$\\phi_{4}$",
                 "$\\psi_{2}$", "$\\psi_{3}$", "$\\psi_{4}$")
s <- as.data.frame(s)
colnames(s) <- c("mean", "sd", "lower", "upper", "n_eff", "Rhat")

print(xtable(s), 
      file = "output/logistic_full.txt", 
      only.contents = TRUE,
      sanitize.rownames.function = function(x) {x})



# 2: estimate probability with a 2-stage model

# here we model p1, the prob of recovery, with exp decay
# and p2, the prob n results reproduce, with a logistic function
# we ignore intermediate stages of failure
# we can multiply p1*p2 to produce p(r) and demonstrate nothing changes

x <- d$reproduced[index_sub]

data <- list(n = n, 
             d = as.integer(d$data_available),
             age = as.integer(2018 - d$year),
             type = as.integer(as.factor(d$type)),
             n_sub = n_sub,
             index_sub = index_sub,
             x = x,
             n_results = d$n_results[index_sub],
             age_id = age_id, 
             n_age = max(age_id),
             age_dist = age_dist)

fit2 <- stan(file = "robustness-checks-m2.stan", 
             data = data, 
             chains = 4, 
             cores = 4, 
             iter = 1e4)

precis(fit2, 2)

post <- extract.samples(fit2)

r <- apply(post$p, 1, mean) * inv_logit(post$a)

ps <- list(r = r, 
           p2 = inv_logit(post$a),
           p1 = apply(post$p, 1, mean))

s <- precis(ps[3:1], prob = 0.89)[1:4]
rownames(s) <- c("$p_{1}$", "$p_{2}$", "$p(r)$")

print(xtable(s), 
      file = "output/2stage_ps.txt", 
      only.contents = TRUE, 
      sanitize.rownames.function = function(x) {x})

# plot ps

png("output/2_stage_dens.png", 
    res = 250, 
    height = 1000, 
    width = 1200)

par(mar = c(5, 4, 4, 6))

plot(x = NULL,
     y = NULL, 
     ylim = c(0, 30),
     xlim = c(0, 1), 
     xlab = "p",
     ylab = "", 
     yaxt = "n")

offset <- c(0, 10, 20)
cols <- c("orange", "pink", "lightblue")

for (i in 1:3) {
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
            "p(r)")

legend(x = 1.1,
       y = 30,
       legend = legend, 
       fill = sapply(cols[c(3, 2, 1)], function(col) col.alpha(col, 0.8)))

par(xpd = FALSE)

dev.off()



### code analysis

# get the code availability "scripted" codes from the original results
# model reproducibility (p2) single stage logistic regression
# logit(p) = a + a_id + a_code

x <- d$reproduced[index_sub]

data <- list(n_sub = n_sub,
             x = x,
             n_results = d$n_results[index_sub],
             code = d$code[index_sub])

fit_code <- stan(file = "robustness-checks-code-available.stan", 
                 data = data)

post <- extract.samples(fit_code)

traceplot(fit_code, pars = "b_code")

png("output/code_analysis.png", 
    res = 250, 
    height = 1200, 
    width = 2000)

par(mfrow = c(1, 2))

plot(NULL, 
     xlim = c(1, 40), 
     ylim = c(0, 1), 
     xlab = "study id", 
     ylab = "prob. results reproduce")

for (i in 1:40) {
  p <- inv_logit(post$a + post$a_paper[, i] + post$b_code * data$code[i])
  points(i, 
         mean(p), 
         col = ifelse(data$code[i] == 1, "blue", "red"))
  arrows(x0 = i, x1 = i, 
         y0 = HPDI(p)[1], 
         y1 = HPDI(p)[2], 
         length = 0, 
         col = ifelse(data$code[i] == 1, "blue", "red"))
}

points(data$x / data$n_results)

par(xpd = TRUE)
mtext(at = c(1.2), "A")
par(xpd = FALSE)

dens(inv_logit(post$a), col = "red")
dens(inv_logit(post$a + post$b_code), col = "blue", add = TRUE)

par(xpd = TRUE)
mtext(at = c(0.1, 5), "B")
par(xpd = FALSE)

dev.off()

s <- precis(fit_code, prob = 0.89)
rownames(s) <- c("a", "b\\_code", "paper\\_sigma")
print(xtable(s), 
      file = "output/code_tab.txt", 
      only.contents = TRUE, 
      sanitize.rownames.function = function(x) {x})

data <- d[d$key %in% index_sub, ]
data <- data[data$code_available == 1, ]
sum(data$n_results)
sum(data$reproduced)

data <- d[d$key %in% index_sub, ]
data <- data[data$code_available == 0, ]
sum(data$n_results)
sum(data$reproduced)



# TO DO:
# 3: estimate p1, p2, p3, p4 with individual glms, fit with R's stats glm function

# to strip away any mystery in the stan models
# i estimate p1, p2, p3, p4 with the "glm" function from R's stats package
# showing that the results are largely the same

# p1:
# data_available ~ a + b*age

data <- data.frame(available = d$data_available, 
                   age = 2018 - d$year)

m1 <- glm(available ~ age, data = data, family = binomial)

png("output/glm_decay.png",
    res = 250,
    height = 1000, 
    width = 1200)

p1 <- sapply(0:63, function(x) inv_logit(m1$coefficients[1] + m1$coefficients[2]*x))
plot(0:63, 
     p1, 
     type = "l", 
     ylim = c(0, 1), 
     xlab = "years since publication")

dev.off()

# p2, p3, p4

# maybe i need to disaggregate
# figure how to specify + variation by papers

# data_complete ~ a + b*age + random effect for paper
# analysis_complete ~ a + b*age + random effect for paper
# results_reproduced ~ a + b*age + random effect for paper

