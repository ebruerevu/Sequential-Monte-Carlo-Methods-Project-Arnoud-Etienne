#parameters

n <- 1200
dt <- 90/n

theta <- pi
epsilon <- 0.0259
sigma_square <- 0.3238
gamma <- 0.25

sd_step <- sqrt(dt)

x_0  <- 0
x_30 <- 1.49
x_60 <- -5.91
x_90 <- -1.17

# observation locations
obs_idx <- c(401, 801)

#Z number of experiments comparing bootstrap and bridge
Z <- 2^12

#Effective sample size function
ESS_fun <- function(w) {
  s <- sum(w)
  s^2 / sum(w^2)
}

#Normalizing constant z
z_fun <- function(k) {

  tau <- ((n + 1) - k) * dt

  sqrt(2 * pi * sigma_square * tau) *
    (
      exp(-0.5 * sigma_square * tau)
      + 1 + epsilon
    )
}

#True normalising constant obtained from the r script "True normalising constant 4.3.R", we use here  
tnc <- 2.497804e-07

#Mean-squared error of log z^{1:Z}
MSElog <- function(c,t){
	c <- pmax(c, .Machine$double.xmin)
	(1/Z * sum((log(c)-log(tnc))^2))^(-1)*t^-1
}	
#Effective sample size of z^{1:Z}
ESSt <- function(c,t){
	(sum(c))^2/(sum(c^2))*t^{-1}
}

#Now the 96 experiments
#Storage
MSElog_time <- matrix(nrow = 2, ncol = 48)
ESS_time <- matrix(nrow = 2, ncol = 48)

N_values = 2^(5:7)
experiment_id <- 1
for (N in N_values){
h <- N*0.5
	for (r in 1:16){
		#Bootstrap
		ncboot <- numeric(Z)
		boot_time <- system.time({
		for (s in 1:Z){
			#x <- matrix(0, nrow = N, ncol = n + 1)
			# initialise
			#x[,1] <- x_0
			#x[,401] <- x_30
			#x[,801] <- x_60
			#x[,1201] <- x_90
			x_cur <- rep(x_0, N)
			for(k in 1:(n-1)) {

    				mu_prop <- x_cur + sin(x_cur - theta) * dt

				x_new <- rnorm(N, mean = mu_prop, sd = sd_step)

    			if ((k + 1) %in% obs_idx) {

    					if (k + 1 == 401) {
      					x_new <- rep(x_30, N)
    					}

    					if (k + 1 == 801) {
      					x_new <- rep(x_60, N)
    					}

  				}
			#Update particle
			x_cur <- x_new
			
			}
			hat_vals <- round(x_new/(2*pi))*2*pi

			w <- (cos(x_90 - hat_vals)+ 1 + epsilon) *exp(-(x_90 - hat_vals)^2 /(2 * sigma_square * dt)) / z_fun(n)
			
			ncboot[s] <- sum(w)
		}
		})
		mean_time_boot <- boot_time[3] / Z
		MSElog_time[1,experiment_id] <- MSElog(ncboot ,mean_time_boot)
		ESS_time[1,experiment_id] <- ESSt(ncboot ,mean_time_boot)
		
		#Bridge
		ncbridge <- numeric(Z)
		bridge_time <- system.time({
		for (s in 1:Z){
			#x <- matrix(0, nrow = N, ncol = n + 1)

			# ancestors
			#a <- matrix(0L, nrow = N, ncol = n)

			# initialise
			#x[,1] <- x_0
			#x[,401] <- x_30
			#x[,801] <- x_60
			#x[,1201] <- x_90

			# current particles + weights
			x_cur <- rep(x_0, N)
			w_cur <- rep(1 / N, N)
			for (k in 1:(n - 1)) {

  				# -----------------------------------------
  				# Normalise weights
  				# -----------------------------------------

  				w_norm <- w_cur / sum(w_cur)

  				# -----------------------------------------
  				# ESS
  				# -----------------------------------------

  				ess <- ESS_fun(w_cur)

  				# -----------------------------------------
  				# Resampling condition
  				# -----------------------------------------

  				do_resample <- (
    					ess < h &&
    					abs(((k + 1)*dt) %% 1) < 1e-10 #To take into account for floating point errors
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

  				} else {

    				ancestors <- 1:N

    				x_prev <- x_cur

    				w_prev <- w_norm
  				}

  				# store ancestry
  				#a[,k] <- ancestors

  				# -----------------------------------------
  				# Propagation (VECTORISED)
  				# -----------------------------------------

  				mu_prop <- x_prev + sin(x_prev - theta) * dt

  				x_new <- rnorm(N, mean = mu_prop, sd = sd_step)

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

  				}

  				# save trajectories
  					#x[,k + 1] <- x_new

  				# update current particles
  					x_cur <- x_new

  				# -----------------------------------------
  				# Bridge weights
  				# -----------------------------------------

  				tau <- (n - k) * dt

  				hat_vals <- round(x_new / (2 * pi)) * (2 * pi)

  				qvals <- (cos(x_90 - hat_vals)+ 1 + epsilon) *
  						exp(-(x_90 - hat_vals)^2 /(2 * sigma_square * tau)) / z_fun(k)
				qvals <- pmax(qvals, .Machine$double.xmin)

  				# -----------------------------------------
  				# Weight update
  				# -----------------------------------------

				logw <- gamma * log(qvals) + log(w_prev)
  				w_cur <- exp(logw)
			}
			ncbridge[s] <- sum(w_cur)
		}
		})
	mean_time_bridge <- bridge_time[3]/Z
	MSElog_time[2,experiment_id] <- MSElog(ncbridge ,mean_time_bridge)
	ESS_time[2,experiment_id] <- ESSt(ncbridge ,mean_time_bridge)
	cat("Completed experiment:",experiment_id," N =", N,"\n")
	experiment_id <- experiment_id + 1
	}
}

#metric plots
par(mfrow=c(1,2))

plot(MSElog_time[1,1:16],MSElog_time[2,1:16],type = "p", xlim = c(0.001, 0.2), ylim = c(0.001, 0.2), log = 'xy',main = "MSE(log z)^(-1)*Mean(t)^(-1)", xlab = "Bootstrap", ylab = "Bridge", col = "blue")
points(MSElog_time[1,17:21],MSElog_time[2,17:21], col = "red")
#points(MSElog_time[1,33:48],MSElog_time[2,33:48], col = "green")
abline(a = 0, b = 1)

plot(log(ESS_time[1,]),log(ESS_time[2,]),type = "p", xlim = c(11,14), ylim = c(11,14), log = 'xy')
abline(a = 0, b = 1)