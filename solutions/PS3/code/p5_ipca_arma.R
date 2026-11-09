# Problem 5: ARMA(p,q) models for monthly IPCA, model selection and forecasting
# Run from solutions/PS3/:  Rscript code/p5_ipca_arma.R

set.seed(2026)
df <- read.csv("../../problem_sets/PS2/ipca_igpm.csv")
df$date <- as.Date(paste0(df$Date, "-01"))
y <- df$IPCA; T <- length(y)

# a) The series
pdf("figures/p5_ipca.pdf", width = 8, height = 3.6)
par(mar = c(4, 4.2, 1, 1))
plot(df$date, y, type = "l", col = "#0072BC", lwd = 1.2, xlab = "", ylab = "Monthly IPCA (%)")
abline(h = 0, col = "grey60")
dev.off()

# b) Conditional log-likelihood. par = (c, phi_1..phi_p, theta_1..theta_q, sigma^2).
# Every model conditions on the first m0 = 2 observations, with presample shocks set to zero,
# so that all models are fit to the same T - 2 observations and their AIC/BIC are comparable.
m0 <- 2
arma_shocks <- function(par, y, p, q) {
  c0 <- par[1]; phi <- par[1 + seq_len(p)]; theta <- par[1 + p + seq_len(q)]
  e <- numeric(length(y))
  for (t in (m0 + 1):length(y)) {
    e[t] <- y[t] - c0 - sum(phi * y[t - seq_len(p)]) - sum(theta * e[t - seq_len(q)])
  }
  e[(m0 + 1):length(y)]
}
loglik <- function(par, y, p, q) {
  s2 <- par[length(par)]
  if (s2 <= 0) return(-Inf)
  e <- arma_shocks(par[-length(par)], y, p, q)
  -length(e) / 2 * log(2 * pi * s2) - sum(e^2) / (2 * s2)
}

# Stationarity and invertibility: roots of 1 - phi_1 z - ... and 1 + theta_1 z + ... outside the unit circle
roots_ok <- function(phi, theta) {
  ok_ar <- length(phi) == 0 || all(Mod(polyroot(c(1, -phi))) > 1)
  ok_ma <- length(theta) == 0 || all(Mod(polyroot(c(1, theta))) > 1)
  ok_ar && ok_ma
}

# c) Conditional MLE: maximize over (c, phi, theta) with sigma^2 concentrated out,
# from several starting points; then the full Hessian gives the standard errors
fit_arma <- function(y, p, q) {
  n <- length(y) - m0
  conc <- function(b) {
    phi <- b[1 + seq_len(p)]; theta <- b[1 + p + seq_len(q)]
    if (!roots_ok(phi, theta)) return(1e10)
    e <- arma_shocks(b, y, p, q)
    n / 2 * log(mean(e^2))
  }
  # Start: OLS for the AR part on the same sample, zeros or +-0.3 for the MA part
  Z <- embed(y, m0 + 1); X <- cbind(1, Z[, 1 + seq_len(p), drop = FALSE])
  b_ols <- qr.solve(X, Z[, 1])
  starts <- lapply(c(0, 0.3, -0.3), function(s) c(b_ols, rep(s, q)))
  fits <- lapply(starts, function(s) optim(s, conc, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12)))
  best <- fits[[which.min(sapply(fits, `[[`, "value"))]]
  b <- best$par
  s2 <- mean(arma_shocks(b, y, p, q)^2)
  par <- c(b, s2)
  H <- optimHess(par, function(x) loglik(x, y, p, q))
  se <- sqrt(diag(solve(-H)))
  ll <- loglik(par, y, p, q); k <- length(par)
  names(par) <- c("c", if (p) paste0("phi", 1:p), if (q) paste0("theta", 1:q), "sigma2")
  list(p = p, q = q, par = par, se = se, ll = ll, k = k,
       aic = -2 * ll + 2 * k, bic = -2 * ll + k * log(n), n = n,
       spread = diff(range(sapply(fits, `[[`, "value"))))
}

orders <- subset(expand.grid(p = 0:2, q = 0:2), p + q > 0)
fits <- lapply(seq_len(nrow(orders)), function(i) fit_arma(y, orders$p[i], orders$q[i]))

ic <- data.frame(model = sprintf("ARMA(%d,%d)", orders$p, orders$q),
                 k = sapply(fits, `[[`, "k"), loglik = sapply(fits, `[[`, "ll"),
                 AIC = sapply(fits, `[[`, "aic"), BIC = sapply(fits, `[[`, "bic"))
ic$dAIC <- ic$AIC - min(ic$AIC); ic$dBIC <- ic$BIC - min(ic$BIC)
print(ic, digits = 6, row.names = FALSE)
cat("Largest gap between starting points (concentrated -loglik):",
    signif(max(sapply(fits, `[[`, "spread")), 3), "\n")

for (f in fits) {
  cat(sprintf("\nARMA(%d,%d)\n", f$p, f$q))
  print(round(cbind(est = f$par, se = f$se, lo = f$par - 1.96 * f$se, hi = f$par + 1.96 * f$se), 4))
}

# Table for the solutions: one row per model, estimates with 95% intervals below
fmt <- function(f, nm) {
  if (!nm %in% names(f$par)) return(c("", ""))
  i <- which(names(f$par) == nm)
  c(sprintf("$%.3f$", f$par[i]), sprintf("\\scriptsize $[%.2f,\\ %.2f]$", f$par[i] - 1.96 * f$se[i], f$par[i] + 1.96 * f$se[i]))
}
nms <- c("c", "phi1", "phi2", "theta1", "theta2", "sigma2")
rows <- sapply(seq_along(fits), function(j) {
  f <- fits[[j]]
  cells <- sapply(nms, function(nm) fmt(f, nm))
  paste0(ic$model[j], " & ", paste(cells[1, ], collapse = " & "), sprintf(" & %.1f & %.1f \\\\", f$aic, f$bic),
         "\n & ", paste(cells[2, ], collapse = " & "), " & & \\\\[2pt]")
})
writeLines(rows, "tables/p5_models.tex")

# Cross-check with stats::arima, conditional sum of squares on the same sample (n.cond = 2).
# arima reports the mean mu instead of c = mu (1 - phi_1 - ... - phi_p).
cat("\nCross-check with stats::arima (CSS, n.cond = 2): largest coefficient difference per model\n")
for (f in fits) {
  chk <- arima(y, order = c(f$p, 0, f$q), method = "CSS", n.cond = m0)
  phi <- f$par[1 + seq_len(f$p)]
  mine <- c(f$par[-c(1, length(f$par))], f$par[1] / (1 - sum(phi)))
  cat(sprintf("ARMA(%d,%d): %.1e\n", f$p, f$q, max(abs(mine - coef(chk)))))
}

# d) Both criteria pick the same model
stopifnot(which.min(ic$AIC) == which.min(ic$BIC))
ic_best <- fits[[which.min(ic$BIC)]]
cat(sprintf("\nChosen by both AIC and BIC: ARMA(%d,%d)\n", ic_best$p, ic_best$q))
bb <- ic_best$par
cat("Inverse AR roots:", round(Re(1 / polyroot(c(1, -bb[2:3]))), 3),
    " inverse MA root:", round(Re(1 / polyroot(c(1, bb[4]))), 3), "\n")

# A closer look: where does ARMA(2,1) gain over AR(1)? Log-likelihood contributions by period
a1 <- fits[[which(ic$model == "ARMA(1,0)")]]
ll_t <- function(f) { e <- arma_shocks(f$par[-length(f$par)], y, f$p, f$q); s2 <- f$par[length(f$par)]
  -0.5 * log(2 * pi * s2) - e^2 / (2 * s2) }
yrs <- as.integer(format(df$date[(m0 + 1):T], "%Y"))
gain <- tapply(ll_t(ic_best) - ll_t(a1), cut(yrs, c(1994, 1996, 1999, 2009, 2019, 2026),
               labels = c("1995-96", "1997-99", "2000-09", "2010-19", "2020-26")), sum)
cat("Log-likelihood gain of ARMA(2,1) over AR(1) by period:\n"); print(round(gain, 2))

# Robustness: start the sample in January 2000, after the post-Real Plan disinflation
y00 <- y[df$date >= as.Date("2000-01-01")]
f00 <- lapply(seq_len(nrow(orders)), function(i) fit_arma(y00, orders$p[i], orders$q[i]))
ic00 <- data.frame(model = ic$model, loglik = sapply(f00, `[[`, "ll"),
                   dAIC = sapply(f00, `[[`, "aic") - min(sapply(f00, `[[`, "aic")),
                   dBIC = sapply(f00, `[[`, "bic") - min(sapply(f00, `[[`, "bic")))
cat("\nSample from 2000-01 (T =", length(y00), "):\n"); print(ic00, digits = 4, row.names = FALSE)

# Robustness: exact Gaussian likelihood (no presample shocks set to zero), full sample
for (o in list(c(1, 0, 0), c(2, 0, 1))) {
  fm <- arima(y, order = o, method = "ML")
  cat(sprintf("Exact ML ARMA(%d,%d): loglik %.2f, AIC %.2f, BIC %.2f\n", o[1], o[3], fm$loglik, AIC(fm), BIC(fm)))
}

# The model used from here on: the AR(1)
best <- a1
b <- best$par; p <- best$p; q <- best$q
e <- c(rep(0, m0), arma_shocks(b[-length(b)], y, p, q))
mu <- b["c"] / (1 - sum(b[1 + seq_len(p)]))
cat(sprintf("Implied mean %.4f (sample mean %.4f), sigma %.4f\n", mu, mean(y), sqrt(b["sigma2"])))

# e) Fitted values and residuals
fitted <- y - e
idx <- (m0 + 1):T
pdf("figures/p5_fit.pdf", width = 8, height = 6)
par(mfrow = c(2, 1), mar = c(2.5, 4.2, 1.5, 1))
plot(df$date[idx], y[idx], type = "l", col = "grey55", lwd = 1.2, xlab = "", ylab = "Monthly IPCA (%)",
     main = "Actual and fitted values")
lines(df$date[idx], fitted[idx], col = "#0072BC", lwd = 1.2)
legend("topright", bty = "n", lwd = 2, col = c("grey55", "#0072BC"), legend = c("actual", "fitted"))
plot(df$date[idx], e[idx], type = "h", col = "#2D3E50", xlab = "", ylab = "Residual", main = "Residuals")
abline(h = c(-1, 1) * 1.96 * sqrt(b["sigma2"]), lty = 2, col = "#C62828")
dev.off()

res <- e[idx]
cat(sprintf("Residuals: sd %.3f; skewness %.2f; kurtosis %.2f; share outside +-1.96 sigma %.3f\n",
            sd(res), mean((res - mean(res))^3) / sd(res)^3, mean((res - mean(res))^4) / sd(res)^4,
            mean(abs(res) > 1.96 * sqrt(b["sigma2"]))))
yrs <- as.integer(format(df$date[idx], "%Y"))
cat("Residual sd by period:", sprintf("%s: %.3f", c("1995-99", "2000-09", "2010-19", "2020-26"),
    tapply(res, cut(yrs, c(1994, 1999, 2009, 2019, 2026)), sd)), "\n")

# f) Residual ACF
r_acf <- acf(res, lag.max = 24, plot = FALSE)$acf[-1]
pdf("figures/p5_resid_acf.pdf", width = 7, height = 3.6)
par(mar = c(4, 4.2, 1, 1))
plot(1:24, r_acf, type = "h", lwd = 3, col = "#2D3E50", ylim = c(-0.25, 0.3), xlab = "Lag (months)",
     ylab = "Residual autocorrelation", xaxt = "n")
axis(1, at = c(1, 6, 12, 18, 24))
abline(h = 0, col = "grey40"); abline(h = c(-1, 1) * 1.96 / sqrt(length(res)), lty = 2, col = "grey50")
dev.off()
cat("Residual ACF, lags 1-13:", round(r_acf[1:13], 3), "\n")
lb <- Box.test(res, lag = 24, type = "Ljung-Box", fitdf = p + q)
cat(sprintf("Ljung-Box(24): stat %.1f, p-value %.4f\n", lb$statistic, lb$p.value))

# g)-h) Simulated forecast paths: iterate the ARMA forward with N(0, sigma^2) shocks
H <- 24
sim_path <- function() {
  yy <- c(y, numeric(H)); ee <- c(e, rnorm(H, 0, sqrt(b["sigma2"])))
  for (t in (T + 1):(T + H)) {
    yy[t] <- b["c"] + sum(b[1 + seq_len(p)] * yy[t - seq_len(p)]) +
      sum(b[1 + p + seq_len(q)] * ee[t - seq_len(q)]) + ee[t]
  }
  yy[(T + 1):(T + H)]
}
fdates <- seq(df$date[T], by = "month", length.out = H + 1)[-1]
show <- df$date >= as.Date("2020-01-01")

one <- sim_path()
pdf("figures/p5_one_path.pdf", width = 8, height = 3.6)
par(mar = c(4, 4.2, 1, 1))
plot(df$date[show], y[show], type = "l", col = "#0072BC", lwd = 1.5, xlim = range(df$date[show], fdates),
     ylim = range(y[show], one), xlab = "", ylab = "Monthly IPCA (%)")
lines(c(df$date[T], fdates), c(y[T], one), col = "#EB811B", lwd = 1.5)
legend("topright", bty = "n", lwd = 2, col = c("#0072BC", "#EB811B"), legend = c("data", "one simulated path"))
dev.off()

paths <- replicate(1000, sim_path())                       # H x 1000
bands <- apply(paths, 1, quantile, c(0.05, 0.5, 0.95))
fmean <- rowMeans(paths)
# Exact conditional mean of the forecast, for comparison with the simulated mean
point <- { yy <- c(y, numeric(H)); ee <- c(e, numeric(H))
  for (t in (T + 1):(T + H)) yy[t] <- b["c"] + sum(b[1 + seq_len(p)] * yy[t - seq_len(p)]) +
    sum(b[1 + p + seq_len(q)] * ee[t - seq_len(q)]); yy[(T + 1):(T + H)] }
pdf("figures/p5_fan.pdf", width = 8, height = 3.8)
par(mar = c(4, 4.2, 1, 1))
plot(df$date[show], y[show], type = "l", col = "#0072BC", lwd = 1.5, xlim = range(df$date[show], fdates),
     ylim = range(y[show], bands), xlab = "", ylab = "Monthly IPCA (%)")
polygon(c(fdates, rev(fdates)), c(bands[1, ], rev(bands[3, ])), col = adjustcolor("#EB811B", 0.25), border = NA)
lines(c(df$date[T], fdates), c(y[T], fmean), col = "#EB811B", lwd = 2)
abline(h = mu, lty = 3, col = "grey40")
legend("topright", bty = "n", lwd = c(2, 2, 8, 1), lty = c(1, 1, 1, 3),
       col = c("#0072BC", "#EB811B", adjustcolor("#EB811B", 0.25), "grey40"),
       legend = c("data", "mean of simulated paths", "5th-95th percentiles", "unconditional mean"))
dev.off()
cat("\nForecast at h = 1, 3, 12, 24 (mean, 5%, 95%):\n")
print(round(rbind(mean = fmean, exact = point, bands)[, c(1, 3, 12, 24)], 3))
cat(sprintf("Last observation %s: %.2f; last residual %.3f\n", df$Date[T], y[T], e[T]))
phi_b <- b[1 + seq_len(p)]; theta_b <- b[1 + p + seq_len(q)]
psi <- c(1, ARMAtoMA(ar = phi_b, ma = theta_b, lag.max = 2000))
g0 <- b["sigma2"] * sum(psi^2)
cat(sprintf("Unconditional 5-95%% range of y under the model: [%.3f, %.3f]\n",
            mu - 1.645 * sqrt(g0), mu + 1.645 * sqrt(g0)))
cat("psi weights 0-6, 12, 24:", round(psi[c(1:7, 13, 25)], 3), "\n")

# Implied autocorrelations of AR(1) and ARMA(2,1) against the sample ACF
acf_y <- acf(y[-(1:m0)], lag.max = 24, plot = FALSE)$acf[-1]
acf_ar1 <- b["phi1"]^(1:24)
acf_best <- ARMAacf(ar = bb[2:3], ma = bb[4], lag.max = 24)[-1]
pdf("figures/p5_acf_models.pdf", width = 7, height = 3.6)
par(mar = c(4, 4.2, 1, 1))
plot(1:24, acf_y, type = "h", lwd = 3, col = "grey60", ylim = c(-0.1, 0.8), xlab = "Lag (months)",
     ylab = "Autocorrelation", xaxt = "n")
axis(1, at = c(1, 6, 12, 18, 24)); abline(h = 0, col = "grey40")
lines(1:24, acf_ar1, type = "b", pch = 16, col = "#EB811B", lwd = 2)
lines(1:24, acf_best, type = "b", pch = 16, col = "#0072BC", lwd = 2)
legend("topright", bty = "n", lwd = c(3, 2, 2), pch = c(NA, 16, 16), col = c("grey60", "#EB811B", "#0072BC"),
       legend = c("sample ACF of IPCA", "implied by AR(1)", "implied by ARMA(2,1)"))
dev.off()
cat("ACF at lags 1, 3, 6, 9, 12 - sample:", round(acf_y[c(1, 3, 6, 9, 12)], 3),
    " AR(1):", round(acf_ar1[c(1, 3, 6, 9, 12)], 3), " ARMA(2,1):", round(acf_best[c(1, 3, 6, 9, 12)], 3), "\n")
cat("ARMA(2,1) psi weights 1-6:", round(ARMAtoMA(ar = bb[2:3], ma = bb[4], lag.max = 6), 3), "\n")
