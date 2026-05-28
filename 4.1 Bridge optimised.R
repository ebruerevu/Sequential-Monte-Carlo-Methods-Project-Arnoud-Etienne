#parameters 4.1
theta_1 <- 0.0187
theta_2 <- 0.2610
theta_3 <- 0.0224
theta_ratio <- theta_1 / theta_2
theta_ratio_2_3 <- theta_3^2/(2 * theta_2)

#Nr. of time steps & step size
N <- 2^5
n <- 100 
dt <- 1/n

steps <- 1:n
exp_forward <- exp(-theta_2 * dt * steps) #one time step, exp() are expensive to continuously compute
sigma_vec <- sqrt(theta_ratio_2_3 *(1 - exp(-2 * theta_2 * dt * steps)))

#Known x-values
x_0 <- 0.07
x_n <- 0.15

#Effective sample size
ESS_fun <- function(w) {
  (sum(w)^2) / sum(w^2)
}

x_cur <- rep(x_0, N)
x <- matrix(0, N, n + 1)
x[,1] <- x_0
x[,n+1] <- x_n
w_cur <- rep(1 / N, N)
a <- matrix(0, N, n + 1)
a[,1] <- 1:N
for (k in 1:(n-1)){
	RESS <- ESS_fun(w_cur)
  	do_resample <- (
          	RESS < h &&
          	(k-1) %% 10 == 0
        )
	if (do_resample){
      	ancestors = sample.int(N, size=N, replace = TRUE, prob=w_cur)
		a[,k] <- ancestors
				
		x_prev <- x_cur[ancestors]
				
      	w_prev <- rep(1 / N, N)
   	}
	else {
		a[,k] <- 1:N
      	x_prev <- x_cur

            w_prev <- w_cur / sum(w_cur)
	}
	# propagate
      mu_prop <- theta_ratio + (x_prev - theta_ratio) * exp_forward[1]
	x_new <- rnorm(N, mean = mu_prop, sd   = sigma_vec[1])
	x[, k + 1] <- x_new
	x_cur <- x_new
	remain <- n + 1 - k

      mu_future <- theta_ratio + (x_new - theta_ratio) * exp_forward[remain]
	w_cur <- dnorm(x_n, mean = mu_future, sd   = sigma_vec[remain]) * w_prev
}
