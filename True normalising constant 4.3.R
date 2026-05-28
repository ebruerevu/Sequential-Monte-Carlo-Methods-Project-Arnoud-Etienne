# ============================================================
# Parameters
# ============================================================
# 2^20 is a too large number for rstudio to handle, so we opted for 2^16

TN <- 2^16
n  <- 1200
dt <- 90 / n

theta <- pi

epsilon <- 0.1
sigma_square <- 1

# observations
obs_times <- c(1, 401, 801, 1201)

obs_values <- c(
  0,
  1.49,
  -5.91,
  -1.17
)

# ============================================================
# Constants
# ============================================================

sd_step <- sqrt(dt)

# ============================================================
# Initial particles
# ============================================================

x_cur <- rep(obs_values[1], TN)

w_cur <- rep(1 / TN, TN)

# ============================================================
# Helper functions
# ============================================================

hat_x <- function(x) {
  round(x / (2*pi)) * 2*pi
}

tz <- function(k) {

  tau <- (n + 1 - k) * dt

  sqrt(2*pi*epsilon*tau) *
    (
      exp(-0.5 * sigma_square * tau)
      + 1 + epsilon
    )
}

tq_vec <- function(x, obs, k) {

  tau <- (n + 1 - k) * dt

  h <- hat_x(x)

  (
    cos(obs - h) + 1 + epsilon
  ) *
    exp(
      -(obs - h)^2 /
      (2 * sigma_square * tau)
    ) / tz(k)
}

# ============================================================
# Estimate normalising constant
# ============================================================

logZ <- 0

for(k in 1:n) {

  # ----------------------------------------------------------
  # Propagate
  # ----------------------------------------------------------

  mu_prop <- x_cur +
    sin(x_cur - theta) * dt

  x_cur <- rnorm(
    TN,
    mean = mu_prop,
    sd   = sd_step
  )

  # ----------------------------------------------------------
  # Observation update
  # ----------------------------------------------------------

  obs_match <- which(obs_times == (k + 1))

  if(length(obs_match) > 0) {

    obs <- obs_values[obs_match]

    w_increment <- tq_vec(
      x_cur,
      obs,
      k
    )

    # normalising constant estimate
    logZ <- logZ + log(mean(w_increment))

    # update weights
    w_cur <- w_increment

    # normalize
    w_cur <- w_cur / sum(w_cur)

    # ESS
    ess <- 1 / sum(w_cur^2)

    # optional resampling
    if(ess < 0.5 * TN) {

      ancestors <- sample.int(
        TN,
        size = TN,
        replace = TRUE,
        prob = w_cur
      )

      x_cur <- x_cur[ancestors]

      w_cur <- rep(1/TN, TN)
    }
  }
}

tnc <- exp(logZ)

print(tnc)