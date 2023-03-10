# Build gp model for logistic regression
gp_model <- function(dat, auc_var = "auc_cros", est = NULL){
  # dat: dataset contains sample size `n` and corresponding C-statistics `auc_var`
  # auc_var: variable name of contains C-statistics
  # est: estimated posterior summary of external data

  require(R2jags, quietly = TRUE)

  dat <- na.omit(dat)

  dat$auc <- dat[[auc_var]]
  prams <- c("alpha", "beta", "gamma", "sigma_g", "phi", "sigma_y")

  if(is.null(est)){
    dat_lst <- list(x = dat$n,
                    y = dat$auc,
                    d = dist(dat$n) |>  as.matrix(),
                    N = nrow(dat),
                    w = dat$n/max(dat$n))
    mod_file <- "gp-mod.txt"
  }else {
    # With priori
    # Moment matching
    est_beta_param <- function(val) {
      mu <- val[1]
      var <- val[2]
      alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
      beta <- alpha * (1 / mu - 1)
      c(alpha, beta)
    }

    est <- cbind(est, var = est[,"sd"]^2)
    # est <- cbind(est, var = 1) # Set all variance to 1
    est <- est[c("beta", "gamma", "phi", "sigma_g"), c("mean", "var")]

    dat_lst <- list(x = dat$n,
                    y = dat$auc,
                    d = dist(dat$n) |>  as.matrix(),
                    N = nrow(dat),
                    w = dat$n/max(dat$n),
                    # a = c(1/2, 1/2),
                    b = c(est["beta", 1], 1/est["beta", 2]),
                    # b = c(est["beta", 1], 1),
                    c = est_beta_param(est["gamma", ]),
                    s = c(est["sigma_g", ][1], 1/est["sigma_g", ][2]),
                    # s = c(est["sigma_g", ][1], 1),
                    p = c(est["phi", ][1], 1/est["phi", ][2]))
    
    mod_file <- "gp-mod-prior.txt"
  }
  
  jags(
    data = dat_lst,
    parameters.to.save = prams,
    model.file = mod_file,
    n.chains = 1, 
    n.iter = 20000, 
    n.burnin = 10000,
    n.thin = 10,
    quiet = TRUE
  )
}

# Summarise prediction and data
predict_gp <- function(fit, dat, auc_var = "auc_cros"){
  # fit: fitted GP model
  # dat: observed dataset contains sample size `n` and corresponding C-statistics `auc_var`
  # auc_var: variable name of contains C-statistics

  dat$auc <- dat[[auc_var]]
  
  dat <- dat[,.(auc=mean(auc, na.rm = TRUE)), .(n)]
  dfm <- data.frame(x = dat$n,
                    y = dat$auc)

  x_star <- seq(50, 2000, by = 10)
  pred_mcmc(fit, x_star = x_star, dfm)  |> 
    as_tibble() |> 
    rename(n = x) 
}


# Predict
pred_mcmc <- function(mod, x_star, dat){

  alpha <- mod$BUGSoutput$summary["alpha", "mean"]
  beta <- mod$BUGSoutput$summary["beta", "mean"]
  gamma <- mod$BUGSoutput$summary["gamma", "mean"]
  sigma_y <- mod$BUGSoutput$summary["sigma_y", "mean"]
  sigma_g <- mod$BUGSoutput$summary["sigma_g", "mean"]
  phi <- mod$BUGSoutput$summary["phi", "mean"]

  x <- dat$x
  y <- dat$y

  n_obs <- length(x)
  
  mu <- 1 - alpha - beta * x^(-gamma)
  mu_pred <- 1 - alpha - beta * x_star^(-gamma)
  
  Sigma <-  sigma_y*2 * diag(length(x)) + sigma_g^2*exp(-(phi^2)*fields::rdist(x,x)^2)
  Sigma_star <- sigma_g^2*exp(-(phi^2)*fields::rdist(x_star,x)^2)
  Sigma_star_star <- sigma_g^2*exp(-(phi^2)*fields::rdist(x_star,x_star)^2)

  pred_mean <- mu_pred + Sigma_star %*% solve(Sigma) %*% (y - mu)
  pred_var <- Sigma_star_star - Sigma_star %*% solve(Sigma) %*% t(Sigma_star)
  
  data.frame(x = x_star, 
             fit = pred_mean,
             lwr = pred_mean - 1.96 * sqrt(diag(pred_var)), 
             upr = pred_mean + 1.96 * sqrt(diag(pred_var)))
}

