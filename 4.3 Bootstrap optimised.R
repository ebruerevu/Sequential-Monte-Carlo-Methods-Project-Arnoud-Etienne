N <- 2^5
n <- 1200

dt <- 90/n

theta <- pi
epsilon <- 0.0259
sigma_square <- 0.3238

sd_step <- sqrt(dt)

x <- matrix(0, N, n+1)

x[,1] <- 0
x[,401] <- 1.49
x[,801] <- -5.91
x[,1201] <- -1.17

obs_idx <- c(401,801,1201)

z <- function(k){
    tau <- (n + 1 - k) * dt

    sqrt(2*pi*epsilon*tau) *
    (
        exp(-0.5 * sigma_square * tau)
        + 1 + epsilon
    )
}

for(k in 1:(n-1)){

    idx <- which(x[,k+1] == 0)

    mu_prop <- x[idx,k] +
        sin(x[idx,k] - theta) * dt

    x[idx,k+1] <- rnorm(
        length(idx),
        mean = mu_prop,
        sd = sd_step
    )
}

hat_vals <- round(x[,n]/(2*pi))*2*pi

w <- (
    cos(x[,n+1] - hat_vals)
    + 1 + epsilon
) *
exp(
    -(x[,n+1] - hat_vals)^2 /
    (2 * sigma_square * dt)
) / z(n)

nc <- sum(w)

plot(
    x[1,],
    type = "l",
    ylim = c(-15,15),
    xlab = "t",
    ylab = "x"
)

for(i in 2:N){
    lines(x[i,])
}