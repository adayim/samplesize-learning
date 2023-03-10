# Fit NSL model
nls_model <- function(dat, auc_var = "auc_cros") {
  # dat: dataset contains sample size `n` and corresponding C-statistics `auc_var`
  # auc_var: variable name of contains C-statistics

  require(boot, quietly = TRUE)
  require(gslnls, quietly = TRUE)

  dat <- na.omit(dat)
  dat$auc <- dat[[auc_var]]

  # Calculate weight
  if(!"w" %in% names(dat))
    dat$w <- dat$n / max(dat$n)

  # Fit model
  gsl_nls(
    fn = auc ~ (1 - inv.logit(alpha)) - (beta * (n^(-inv.logit(gamma)))),
    data = dat,
    start = list(alpha = -1, beta = 3, gamma = 2),
    algorithm = "lm",
    control = gsl_nls_control(maxiter = 10000),
    weights = w
  )
  
}

## Prediction ==============
# Confidence/prediction bands for nonlinear regression 
# (i.e., objects of class nls) are based on a linear approximation
# as described in Bates & Watts (2007). Also known as Delta method.
predict_nls <- function(fit, dat, auc_var = "auc_cros"){
  # fit: fitted NLS model
  # dat: observed dataset contains sample size `n` and corresponding C-statistics `auc_var`
  # auc_var: variable name of contains C-statistics

  dat$auc <- dat[[auc_var]]
  
  new_data <- data.frame(n = seq(50, 2000, by = 10))
  
  predict(fit, newdata = new_data, interval = "prediction", level= 0.95)  %>% 
    as_tibble() %>% 
    mutate(n = new_data$n)
}

