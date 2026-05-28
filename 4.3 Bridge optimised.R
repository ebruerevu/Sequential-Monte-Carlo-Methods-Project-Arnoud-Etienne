# =========================================================
# Bridge Particle Filter (Vectorised + Faster)
# =========================================================

# -------------------------
# Parameters
# -------------------------

N <- 2^5
n <- 1200

x_0  <- 0
x_30 <- 1.49
x_60 <- -5.91
x_90 <- -1.17

theta <- pi
epsilon <- 0.0259
sigma_square <- 0.3238

dt <- 90 / n
sd_step <- sqrt(dt)

gamma <- 0.25
h <- N * 0.5

# observation locations
obs_idx <- c(401, 801, 1201)

# -------------------------------------------------
# Storage
# -------------------------------------------------

x <- matrix(0, nrow = N, ncol = n + 1)

# ancestors
a <- matrix(0L, nrow = N, ncol = n)

# initialise
x[,1] <- x_0
x[,401] <- x_30
x[,801] <- x_60
x[,1201] <- x_90

# current particles + weights
x_cur <- rep(x_0, N)
w_cur <- rep(1 / N, N)

# -------------------------------------------------
# Helper functions
# -------------------------------------------------

ESS_fun <- function(w) {
  s <- sum(w)
  s^2 / sum(w^2)
}

z_fun <- function(k) {

  tau <- (n + 1 - k) * dt

  sqrt(2 * pi * epsilon * tau) *
    (
      exp(-0.5 * sigma_square * tau)
      + 1 + epsilon
    )
}

# -------------------------------------------------
# Main Particle Filter
# -------------------------------------------------

for (k in 1:(n - 1)) {

  # -----------------------------------------
  # Normalise weights
  # -----------------------------------------

  w_norm <- w_cur / sum(w_cur)

  # -----------------------------------------
  # ESS
  # -----------------------------------------

  ess <- ESS_fun(w_norm)

  # -----------------------------------------
  # Resampling condition
  # -----------------------------------------

  do_resample <- (
    ess < h &&
    (((k + 1)*dt) %% 1 == 0)
  )

  # -----------------------------------------
  # Resampling
  # -----------------------------------------

  if (do_resample) {

    ancestors <- sample.int(
      N,
      size = N,
      replace = TRUE,
      prob = w_norm
    )

    x_prev <- x_cur[ancestors]

    w_prev <- rep(1 / N, N)

    cat("Resample at k =", k, "\n")

  } else {

    ancestors <- 1:N

    x_prev <- x_cur

    w_prev <- w_norm
  }

  # store ancestry
  a[,k] <- ancestors

  # -----------------------------------------
  # Propagation (VECTORISED)
  # -----------------------------------------

  mu_prop <- x_prev +
    sin(x_prev - theta) * dt

  x_new <- rnorm(
    N,
    mean = mu_prop,
    sd = sd_step
  )

  # -----------------------------------------
  # Respect observations
  # -----------------------------------------

  if ((k + 1) %in% obs_idx) {

    if (k + 1 == 401) {
      x_new <- rep(x_30, N)
    }

    if (k + 1 == 801) {
      x_new <- rep(x_60, N)
    }

    if (k + 1 == 1201) {
      x_new <- rep(x_90, N)
    }
  }

  # save trajectories
  x[,k + 1] <- x_new

  # update current particles
  x_cur <- x_new

  # -----------------------------------------
  # Bridge weights
  # -----------------------------------------

  tau <- (n + 1 - k) * dt

  hat_vals <- round(x_new / (2 * pi)) * (2 * pi)

  qvals <- (
    cos(x_90 - hat_vals)
    + 1 + epsilon
  ) *
  exp(
    -(x_90 - hat_vals)^2 /
    (2 * sigma_square * tau)
  ) / z_fun(k)

  # -----------------------------------------
  # Weight update
  # -----------------------------------------

  w_cur <- (qvals ^ gamma) * w_prev
}

# =========================================================
# Estimate normalising constant
# =========================================================

nc <- sum(w_cur)

cat("Estimated normalising constant =", nc, "\n")

# =========================================================
# Plot particle trajectories
# =========================================================

plot(
  x[1,],
  type = "l",
  ylim = c(-15, 15),
  xlab = "t",
  ylab = "x",
  col = "blue"
)

for (i in 2:N) {
  lines(x[i,], col = rgb(0,0,1,0.3))
}

# =========================================================
# Reconstruct ancestral lines
# =========================================================

b <- matrix(0L, nrow = N, ncol = n - 1)

for (i in 1:N) {

  b[i,n-1] <- i

  for (k in 2:(n - 1)) {

    b[i,n-k] <-
      a[
        b[i,n-k+1],
        n-k+1
      ]
  }
}