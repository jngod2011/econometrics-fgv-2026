# Problem 3 d): the ARMA(1,1) log-likelihood is flat along theta = -phi
# Run from solutions/PS3/:  Rscript code/p3_arma11_likelihood.R

set.seed(2026)
T <- 500
grid <- seq(-0.9, 0.9, by = 0.02)
G <- expand.grid(phi = grid, theta = grid)

# Log-likelihood with sigma^2 = 1 on the whole grid at once: back out the shocks with
# e_t = y_t - phi y_{t-1} - theta e_{t-1}, starting from y_0 = e_0 = 0
loglik_grid <- function(y, phi, theta) {
  e <- numeric(length(phi)); ylag <- 0; ssr <- numeric(length(phi))
  for (t in seq_along(y)) {
    e <- y[t] - phi * ylag - theta * e
    ssr <- ssr + e^2
    ylag <- y[t]
  }
  -length(y) / 2 * log(2 * pi) - ssr / 2
}

heat <- function(ll, main, truth, file) {
  z <- matrix(pmax(ll - max(ll), -60), length(grid), length(grid))   # rows: phi, columns: theta
  best <- G[which.max(ll), ]
  pdf(file, width = 5.6, height = 5)
  par(mar = c(4.2, 4.2, 2, 1))
  image(grid, grid, z, col = hcl.colors(60, "viridis"), zlim = c(-60, 0), useRaster = TRUE,
        xlab = expression(phi), ylab = expression(theta), main = main)
  abline(0, -1, col = "white", lwd = 2, lty = 2)
  points(truth[1], truth[2], pch = 4, cex = 1.6, lwd = 3, col = "#EB811B")
  points(best$phi, best$theta, pch = 1, cex = 1.6, lwd = 2, col = "white")
  legend("bottomleft", bg = "white", box.col = "grey70", lty = c(2, NA, NA), pch = c(NA, 4, 1),
         col = c("grey30", "#EB811B", "grey30"), lwd = 2, cex = 0.85,
         legend = c(expression(theta == -phi), "true", "grid maximum"))
  dev.off()
  best
}

# White noise data: phi = theta = 0
y_wn <- rnorm(T)
ll_wn <- loglik_grid(y_wn, G$phi, G$theta)
best_wn <- heat(ll_wn, "White noise data", c(0, 0), "figures/p3_loglik_wn.pdf")

line_pts <- c(-0.8, -0.4, 0, 0.4, 0.8)
ll_line <- loglik_grid(y_wn, line_pts, -line_pts)
cat("White noise data: log-likelihood on the line theta = -phi\n")
print(data.frame(phi = line_pts, theta = -line_pts, loglik = sprintf("%.10f", ll_line)))
cat(sprintf("Max over the grid: %.4f at (phi, theta) = (%.2f, %.2f)\n", max(ll_wn), best_wn$phi, best_wn$theta))
cat(sprintf("Max over the line (white noise loglik): %.4f\n", ll_line[3]))

# ARMA(1,1) data with phi = 0.5, theta = 0.45 (arima.sim uses a burn-in)
y_arma <- as.numeric(arima.sim(model = list(ar = 0.5, ma = 0.45), n = T))
ll_arma <- loglik_grid(y_arma, G$phi, G$theta)
best_arma <- heat(ll_arma, "ARMA(1,1) data, (0.5, 0.45)", c(0.5, 0.45), "figures/p3_loglik_arma.pdf")
ll_line2 <- loglik_grid(y_arma, line_pts, -line_pts)
cat("\nARMA(1,1) data: log-likelihood on the line theta = -phi\n")
print(data.frame(phi = line_pts, theta = -line_pts, loglik = sprintf("%.4f", ll_line2)))
cat(sprintf("Max over the grid: %.4f at (phi, theta) = (%.2f, %.2f); line value minus max: %.1f\n",
            max(ll_arma), best_arma$phi, best_arma$theta, ll_line2[3] - max(ll_arma)))

# c) standard errors of phi-hat from the information matrix
se_phi <- function(phi, theta, T = 200) sqrt((1 - phi^2) * (1 + phi * theta)^2 / (phi + theta)^2 / T)
cat(sprintf("\nc) se(phi-hat), T = 200: (0.5, 0.45): %.4f; (0.5, -0.45): %.4f; AR(1) with phi = 0.5: %.4f\n",
            se_phi(0.5, 0.45), se_phi(0.5, -0.45), sqrt(0.75 / 200)))
