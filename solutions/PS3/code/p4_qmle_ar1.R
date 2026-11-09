# Problem 4: Hessian, OPG and sandwich standard errors for the Gaussian QMLE of an AR(1)
# Run from solutions/PS3/:  Rscript code/p4_qmle_ar1.R

set.seed(2026)

# b) Estimates and the three standard errors for Theta = (c, phi, sigma^2)
qmle_ar1 <- function(y) {
  n <- length(y) - 1                      # T - 1 observations, conditioning on y_1
  X <- cbind(1, y[-length(y)]); yy <- y[-1]
  b <- solve(crossprod(X), crossprod(X, yy))
  e <- drop(yy - X %*% b); s2 <- mean(e^2)
  # Scores, one row per t
  S <- cbind(X * e / s2, -1 / (2 * s2) + e^2 / (2 * s2^2))
  # Average Hessian
  H <- matrix(0, 3, 3)
  H[1:2, 1:2] <- -crossprod(X) / (n * s2)
  H[1:2, 3] <- H[3, 1:2] <- -colSums(X * e) / (n * s2^2)
  H[3, 3] <- 1 / (2 * s2^2) - mean(e^2) / s2^3
  O <- crossprod(S) / n
  Hi <- solve(H)
  se <- cbind(hessian = sqrt(diag(solve(-H)) / n),
              opg = sqrt(diag(solve(O)) / n),
              sandwich = sqrt(diag(Hi %*% O %*% Hi) / n))
  list(theta = c(b, s2), se = se)
}

sim_ar1 <- function(T, phi, shocks, burn = 200) {
  e <- shocks(T + burn)
  as.numeric(stats::filter(e, phi, method = "recursive"))[(burn + 1):(burn + T)]
}
gauss <- function(n) rnorm(n)
t7 <- function(n) sqrt(5 / 7) * rt(n, df = 7)      # unit variance, kurtosis 5

truth <- c(0, 0.5, 1); M <- 2000
res <- NULL
for (dist in c("Normal", "t7")) {
  for (T in c(200, 1000)) {
    cover <- array(NA, c(M, 3, 3))
    for (m in 1:M) {
      y <- sim_ar1(T, 0.5, if (dist == "Normal") gauss else t7)
      f <- qmle_ar1(y)
      cover[m, , ] <- abs(f$theta - truth) <= 1.96 * f$se
    }
    cv <- apply(cover, c(2, 3), mean)
    res <- rbind(res, data.frame(dist = dist, T = T, param = c("c", "phi", "sigma2"),
                                 hessian = cv[, 1], opg = cv[, 2], sandwich = cv[, 3]))
  }
}
print(res, digits = 3, row.names = FALSE)

# Limiting coverage for sigma^2 with kurtosis kappa: the true variance is (kappa - 1) sigma^4
kappa <- 5
cat(sprintf("\nPredicted limiting coverage for sigma^2 (kappa = 5): Hessian %.3f, OPG %.3f\n",
            2 * pnorm(1.96 * sqrt(2 / (kappa - 1))) - 1, 2 * pnorm(1.96 * sqrt(4 / (kappa - 1)^2)) - 1))

rows <- with(res, sprintf("%s & %d & $%s$ & %.1f & %.1f & %.1f \\\\",
                          c(Normal = "Normal", t7 = "$t_7$")[dist], T,
                          c(c = "c", phi = "\\phi", sigma2 = "\\sigma^2")[param],
                          100 * hessian, 100 * opg, 100 * sandwich))
writeLines(c(rows[1:6], "\\midrule", rows[7:12]), "tables/p4_coverage.tex")
