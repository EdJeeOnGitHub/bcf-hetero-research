#' @importFrom stats approxfun lm qchisq quantile sd
#' @importFrom RcppParallel RcppParallelLibs
Rcpp::loadModule(module = "TreeSamples", TRUE)

.ident <- function(...){
# courtesy https://stackoverflow.com/questions/19966515/how-do-i-test-if-three-variables-are-equal-r
  args <- c(...)
  if( length( args ) > 2L ){
    #  recursively call ident()
    out <- c( identical( args[1] , args[2] ) , .ident(args[-1]))
  }else{
    out <- identical( args[1] , args[2] )
  }
  return( all( out ) )
}

.cp_quantile = function(x, num=10000, cat_levels=8){
  nobs = length(x)
  nuniq = length(unique(x))

  if(nuniq==1) {
    ret = x[1]
    warning("A supplied covariate contains a single distinct value.")
  } else if(nuniq < cat_levels) {
    xx = sort(unique(x))
    ret = xx[-length(xx)] + diff(xx)/2
  } else {
    q = approxfun(sort(x),quantile(x,p = 0:(nobs-1)/nobs))
    ind = seq(min(x),max(x),length.out=num)
    ret = q(ind)
  }

  return(ret)
}

.get_chain_tree_files = function(tree_path, chain_id, no_output = FALSE){
  if (is.null(tree_path) | no_output){
    out <- list(
                "con_trees" = toString(character(0)), 
                "mod_trees" = toString(character(0)),
                "var_trees" = toString(character(0)),
                "var_ratio_trees" = toString(character(0))
                )
  } else{
    out <- list("con_trees" = paste0(tree_path,'/',"con_trees.", chain_id, ".txt"), 
                "mod_trees" = paste0(tree_path,'/',"mod_trees.", chain_id, ".txt"),
                "var_trees" = paste0(tree_path,'/',"var_trees.", chain_id, ".txt"),
                "var_ratio_trees" = paste0(tree_path,'/',"var_ratio_trees.", chain_id, ".txt"))
  }
  return(out)
}

.get_do_type = function(n_cores, log_file){
  if(n_cores>1){
    cl <- parallel::makeCluster(n_cores, outfile=log_file)
    # Workers don't inherit an in-script .libPaths(), so without this they can
    # load another installed bcf (e.g. CRAN) that has no var_trees support.
    parallel::clusterCall(cl, function(lib_paths) .libPaths(lib_paths), .libPaths())

    message(sprintf("Running in parallel, saving BCF logs to %s \n", log_file))
    doParallel::registerDoParallel(cl)
    `%doType%`  <- foreach::`%dopar%`
  } else {
    cl <- NULL
    `%doType%`  <- foreach::`%do%`
  }
  
  do_type_config <- list('doType'  = `%doType%`,
                         'n_cores' = n_cores,
                         'cluster' = cl)
  
  return(do_type_config)
}

.cleanup_after_par = function(do_type_config){
  if(do_type_config$n_cores>1){
    parallel::stopCluster(do_type_config$cluster)
  }
}

.parse_vartree <- function(vartree, yscale, x_variance) {
  if (is.null(vartree)) vartree <- list()
  num_trees <- if (is.null(vartree$num_trees)) 40L else as.integer(vartree$num_trees)
  if (!is.finite(num_trees) || num_trees <= 0L) stop("vartree$num_trees must be positive")
  nu <- if (is.null(vartree$nu)) 10.0 else as.numeric(vartree$nu)
  if (!is.finite(nu) || nu <= 2) stop("vartree$nu must be finite and greater than 2")
  lambda <- if (is.null(vartree$lambda)) stats::var(yscale) else as.numeric(vartree$lambda)
  if (!is.finite(lambda) || lambda <= 0) lambda <- 1.0
  inv_m <- 1 / num_trees
  nu_tree <- 2 / (1 - (1 - 2 / nu)^inv_m)
  lambda_tree <- lambda^inv_m
  list(
    num_trees = num_trees,
    nu = nu,
    lambda = lambda,
    nu_tree = nu_tree,
    lambda_tree = lambda_tree,
    numcut = if (is.null(vartree$numcut)) 100L else as.integer(vartree$numcut),
    base = if (is.null(vartree$base)) 0.95 else as.numeric(vartree$base),
    power = if (is.null(vartree$power)) 2.0 else as.numeric(vartree$power),
    sparse = if (is.null(vartree$sparse)) FALSE else isTRUE(vartree$sparse),
    x_variance_ncol = ncol(x_variance)
  )
}

.wtd_var <- function(x, w) {
  w <- as.numeric(w)
  x <- as.numeric(x)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 1) stop("weights must have sum greater than 1")
  mu <- stats::weighted.mean(x, w)
  sum(w * (x - mu)^2) / (sw - 1)
}

.cara_utility <- function(m, sigma2, alpha, exponent_bounds = c(-700, 700)) {
  exponent <- -alpha * m + 0.5 * alpha^2 * sigma2
  exponent <- pmin(pmax(exponent, exponent_bounds[1]), exponent_bounds[2])
  -exp(exponent)
}

#' Fit Bayesian Causal Forests
#'
#' @references Hahn, Murray, and Carvalho (2020). Bayesian regression tree models for causal inference: regularization, confounding, and heterogeneous effects.
#'  https://projecteuclid.org/journals/bayesian-analysis/volume-15/issue-3/Bayesian-Regression-Tree-Models-for-Causal-Inference--Regularization-Confounding/10.1214/19-BA1195.full. 
#'  (Call citation("bcf") from the command line for citation information in Bibtex format.)
#'
#' @details Fits the Bayesian Causal Forest model (Hahn et. al. 2020): For a response
#' variable y, binary treatment z, and covariates x, we return estimates of mu, tau, and sigma in
#' the model
#' \deqn{y_i = \mu(x_i, \pi_i) + \tau(x_i, \pi_i)z_i + \epsilon_i}
#' where \eqn{\pi_i} is an (optional) estimate of the propensity score \eqn{\Pr(Z_i=1 | X_i=x_i)} and
#' \eqn{\epsilon_i \sim N(0,\sigma^2)}
#'
#' Some notes:
#' \itemize{
#'    \item By default, bcf writes each sample (including the trees in the ensemble) for each chain to a text file, 
#'    which is used for prediction by the predict.bcf function. These text files may be large if bcf is run for many samples, 
#'    so we also provide an option to suppress this output by setting no_output = TRUE. If bcf is run with no_output = TRUE, 
#'    it will not be possible to predict from the model after the fact.
#'    \item x_control and x_moderate must be numeric matrices. See e.g. the makeModelMatrix function in the
#'    dbarts package for appropriately constructing a design matrix from a data.frame
#'    \item sd_control and sd_moderate are the prior SD(mu(x)) and SD(tau(x)) at a given value of x (respectively). If
#'    use_muscale = FALSE, then this is the parameter \eqn{\sigma_\mu} from the original BART paper, where the leaf parameters
#'    have prior distribution \eqn{N(0, \sigma_\mu/m)}, where m is the number of trees.
#'    If use_muscale=TRUE then sd_control is the prior median of a half Cauchy prior for SD(mu(x)). If use_tauscale = TRUE,
#'    then sd_moderate is the prior median of a half Normal prior for SD(tau(x)).
#'    \item By default the prior on \eqn{\sigma^2} is calibrated as in Chipman, George and McCulloch (2010).
#' }
#' @param y Response variable
#' @param z Treatment variable
#' @param x_control Design matrix for the prognostic function mu(x)
#' @param x_moderate Design matrix for the covariate-dependent treatment effects tau(x)
#' @param pihat Length n estimates of propensity score
#' @param w An optional vector of weights. When present, BCF fits a model \eqn{y | x ~ N(f(x), \sigma^2 / w)}, where \eqn{f(x)} is the unknown function.
#' @param random_seed A random seed passed to R's set.seed
#' @param n_chains  An optional integer of the number of MCMC chains to run
#' @param n_threads An optional integer of the number of threads to parallelize within chain bcf operations on
#' @param nburn Number of burn-in MCMC iterations
#' @param nsim Number of MCMC iterations to save after burn-in. The chain will run for nsim*nthin iterations after burn-in
#' @param nthin Save every nthin'th MCMC iterate. The total number of MCMC iterations will be nsim*nthin + nburn.
#' @param update_interval Print status every update_interval MCMC iterations
#' @param ntree_control Number of trees in mu(x)
#' @param sd_control SD(mu(x)) marginally at any covariate value (or its prior median if use_muscale=TRUE)
#' @param base_control Base for tree prior on mu(x) trees (see details)
#' @param power_control Power for the tree prior on mu(x) trees
#' @param ntree_moderate Number of trees in tau(x)
#' @param sd_moderate SD(tau(x)) marginally at any covariate value (or its prior median if use_tauscale=TRUE)
#' @param base_moderate Base for tree prior on tau(x) trees (see details)
#' @param power_moderate Power for the tree prior on tau(x) trees (see details)
#' @param no_output logical, whether to suppress writing trees and training log to text files, defaults to FALSE.
#' @param save_tree_directory Specify where trees should be saved. Keep track of this for predict(). Defaults to working directory. Setting to NULL skips writing of trees.
#' @param log_file file where BCF should save its logs when running multiple chains in parallel. This file is not written too when only running one chain. 
#' @param nu Degrees of freedom in the chisq prior on \eqn{sigma^2}
#' @param lambda Scale parameter in the chisq prior on \eqn{sigma^2}
#' @param sigq Calibration quantile for the chisq prior on \eqn{sigma^2}
#' @param sighat Calibration estimate for the chisq prior on \eqn{sigma^2}
#' @param include_pi Takes values "control", "moderate", "both" or "none". Whether to
#' include pihat in mu(x) ("control"), tau(x) ("moderate"), both or none. Values of "control"
#' or "both" are HIGHLY recommended with observational data.
#' @param use_muscale Use a half-Cauchy hyperprior on the scale of mu.
#' @param use_tauscale Use a half-Normal prior on the scale of tau.
#' @param verbose logical, whether to print log of MCMC iterations, defaults to TRUE.
#' @param x_variance Design matrix for the residual variance function.
#' @param vartree Optional list of scalar product-of-trees variance prior parameters.
#' @param variance_model Residual variance model: \code{"homoscedastic"},
#' \code{"shared"}, or \code{"ratio"}.
#' @param standardize Logical; if \code{TRUE}, center and scale the outcome
#' before fitting and transform posterior draws back to the original scale.
#' Set to \code{FALSE} for prior-predictive calibration with fixed scales.
#' @return A fitted bcf object that is a list with elements
#' \item{tau}{\code{nsim} by \code{n} matrix of posterior samples of individual-level treatment effect estimates}
#' \item{mu}{\code{nsim} by \code{n} matrix of posterior samples of prognostic function E(Y|Z=0, x=x) estimates}
#' \item{sigma}{Length \code{nsim} vector of posterior samples of sigma}
#' @examples
#'\dontrun{
#'
#' # data generating process
#' p = 3 #two control variables and one moderator
#' n = 250
#' 
#' set.seed(1)
#'
#' x = matrix(rnorm(n*p), nrow=n)
#'
#' # create targeted selection
#' q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
#'
#' # generate treatment variable
#' pi = pnorm(q)
#' z = rbinom(n,1,pi)
#'
#' # tau is the true (homogeneous) treatment effect
#' tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
#'
#' # generate the response using q, tau and z
#' mu = (q + tau*z)
#'
#' # set the noise level relative to the expected mean function of Y
#' sigma = diff(range(q + tau*pi))/8
#'
#' # draw the response variable with additive error
#' y = mu + sigma*rnorm(n)
#'
#' # If you didn't know pi, you would estimate it here
#' pihat = pnorm(q)
#'
#' bcf_fit = bcf(y, z, x, x, pihat, nburn=2000, nsim=2000)
#'
#' # Get posterior of treatment effects
#' tau_post = bcf_fit$tau
#' tauhat = colMeans(tau_post)
#' plot(tau, tauhat); abline(0,1)
#'
#'}
#'\dontrun{
#'
#' # data generating process
#' p = 3 #two control variables and one moderator
#' n = 250
#' #
#' set.seed(1)
#'
#' x = matrix(rnorm(n*p), nrow=n)
#'
#' # create targeted selection
#' q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
#'
#' # generate treatment variable
#' pi = pnorm(q)
#' z = rbinom(n,1,pi)
#'
#' # tau is the true (homogeneous) treatment effect
#' tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
#'
#' # generate the response using q, tau and z
#' mu = (q + tau*z)
#'
#' # set the noise level relative to the expected mean function of Y
#' sigma = diff(range(q + tau*pi))/8
#'
#' # draw the response variable with additive error
#' y = mu + sigma*rnorm(n)
#'
#' pihat = pnorm(q)
#'
#' # nburn and nsim should be much larger, at least a few thousand each
#' # The low values below are for CRAN.
#' bcf_fit = bcf(y, z, x, x, pihat, nburn=100, nsim=10)
#'
#' # Get posterior of treatment effects
#' tau_post = bcf_fit$tau
#' tauhat = colMeans(tau_post)
#' plot(tau, tauhat); abline(0,1)
#'}
#'
#' @useDynLib bcf
#' @export
bcf <- function(y, z, x_control, x_moderate=x_control, pihat, w = NULL, 
                random_seed = sample.int(.Machine$integer.max, 1),
                n_chains = 4,
                n_threads = max((RcppParallel::defaultNumThreads()-2),1), #max number of threads, minus a arbitrary holdback, over the number of cores
                nburn, nsim, nthin = 1, update_interval = 100,
                ntree_control = 200,
                sd_control = 2*sd(y),
                base_control = 0.95,
                power_control = 2,
                ntree_moderate = 50,
                sd_moderate = sd(y),
                base_moderate = 0.25,
                power_moderate = 3, 
                no_output = FALSE, 
                save_tree_directory = '.',
                log_file=file.path('.',sprintf('bcf_log_%s.txt',format(Sys.time(), "%Y%m%d_%H%M%S"))),
                nu = 3, lambda = NULL, sigq = .9, sighat = NULL,
                include_pi = "control", use_muscale=TRUE, use_tauscale=TRUE, verbose=TRUE,
                x_variance = x_control, vartree = NULL, variance_model = "homoscedastic",
                standardize = TRUE
) {

  
  if(is.null(w)){
    w <- matrix(1, ncol = 1, nrow = length(y))
  }

  pihat = as.matrix(pihat)
  if(!.ident(length(y),
             length(z),
             length(w),
             nrow(x_control),
             nrow(x_moderate),
             nrow(x_variance),
             nrow(pihat))
    ) {
    stop("Data size mismatch. The following should all be equal:
         length(y): ", length(y), "\n",
         "length(z): ", length(z), "\n",
         "length(w): ", length(w), "\n",
         "nrow(x_control): ", nrow(x_control), "\n",
         "nrow(x_moderate): ", nrow(x_moderate), "\n",
         "nrow(x_variance): ", nrow(x_variance), "\n",
         "nrow(pihat): ", nrow(pihat),"\n"
    )
  }

  if(any(is.na(y))) stop("Missing values in y")
  if(any(is.na(z))) stop("Missing values in z")
  if(any(is.na(w))) stop("Missing values in w")
  if(any(is.na(x_control))) stop("Missing values in x_control")
  if(any(is.na(x_moderate))) stop("Missing values in x_moderate")
  if(any(is.na(x_variance))) stop("Missing values in x_variance")
  if(any(is.na(pihat))) stop("Missing values in pihat")
  if(any(!is.finite(y))) stop("Non-numeric values in y")
  if(any(!is.finite(z))) stop("Non-numeric values in z")
  if(any(!is.finite(w))) stop("Non-numeric values in w")
  if(any(!is.finite(x_control))) stop("Non-numeric values in x_control")
  if(any(!is.finite(x_moderate))) stop("Non-numeric values in x_moderate")
  if(any(!is.finite(x_variance))) stop("Non-numeric values in x_variance")
  if(any(!is.finite(pihat))) stop("Non-numeric values in pihat")
  if(!all(sort(unique(z)) == c(0,1))) stop("z must be a vector of 0's and 1's, with at least one of each")
  if(!variance_model %in% c("homoscedastic", "shared", "ratio")) {
    stop("variance_model must be 'homoscedastic', 'shared', or 'ratio'")
  }
  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize)) {
    stop("standardize must be TRUE or FALSE")
  }
  use_ratio <- variance_model == "ratio"
  use_hetero <- variance_model %in% c("shared", "ratio") || !is.null(vartree)
  if(use_hetero && variance_model == "homoscedastic") variance_model <- "shared"

  if(length(unique(y))<5) warning("y appears to be discrete")

  if(nburn<0) stop("nburn must be positive")
  if(nsim<0) stop("nsim must be positive")
  if(nthin<0) stop("nthin must be positive")
  if(nthin>nsim+1) stop("nthin must be < nsim")
  if(nburn<1000) warning("A low (<1000) value for nburn was supplied")
  if(nsim*nburn<1000) warning("A low (<1000) value for total iterations after burn-in was supplied")

  ### TODO range check on parameters

  ###
  x_c = matrix(x_control, ncol=ncol(x_control))
  x_m = matrix(x_moderate, ncol=ncol(x_moderate))
  x_v = matrix(x_variance, ncol=ncol(x_variance))

  if(include_pi=="both" | include_pi=="control") {
    x_c = cbind(x_control, pihat)
  }
  if(include_pi=="both" | include_pi=="moderate") {
    x_m = cbind(x_moderate, pihat)
  }
  cutpoint_list_c = lapply(1:ncol(x_c), function(i) .cp_quantile(x_c[,i]))
  cutpoint_list_m = lapply(1:ncol(x_m), function(i) .cp_quantile(x_m[,i]))

  if (standardize) {
    sdy = sqrt(.wtd_var(y, w))
    muy = stats::weighted.mean(y, w)
    yscale = (y-muy)/sdy
  } else {
    sdy = 1
    muy = 0
    yscale = y
  }
  vartree_params <- .parse_vartree(vartree, yscale, x_v)
  cutpoint_list_v = lapply(1:ncol(x_v), function(i) .cp_quantile(x_v[,i], num = vartree_params$numcut))


  if(is.null(lambda)) {
    if(is.null(sighat)) {
      lmf = lm(yscale~z+as.matrix(x_c), weights = w)
      sighat = summary(lmf)$sigma #sd(y) #summary(lmf)$sigma
    }
    qchi = qchisq(1.0-sigq,nu)
    lambda = (sighat*sighat*qchi)/nu
  }

  dir = tempdir()

  perm = order(z, decreasing=TRUE)

  con_sd = ifelse(abs(2*sdy - sd_control)<1e-6, 2, sd_control/sdy)
  mod_sd = ifelse(abs(sdy - sd_moderate)<1e-6, 1, sd_moderate/sdy)/ifelse(use_tauscale,0.674,1) # if HN make sd_moderate the prior median

  if (n_threads > 1) {
    stop(
      "n_threads > 1 is disabled: on this container RcppParallel/TBB does not ",
      "respect the cgroup CPU quota (RcppParallel::defaultNumThreads() reports ",
      "the host's raw core count, not the quota), and allsuff()/GetSuffBirthWorker/",
      "GetSuffDeathWorker/FitWorker in funs.cpp invoke parallelFor/parallelReduce ",
      "once per tree per MCMC step -- millions of tiny calls per fit. The result is ",
      "not just unhelpful but actively pathological: measured ~5-7x slower wall time ",
      "than n_threads=1, sys time approaching/exceeding user time (thread scheduling ",
      "and malloc-arena contention, not compute), CPU utilization far exceeding the ",
      "requested thread count, and non-deterministic results run-to-run at a fixed ",
      "seed (floating-point summation order depends on TBB's work-stealing schedule). ",
      "Use n_threads = 1 and get parallelism from running more independent processes ",
      "instead (e.g. run_simulation_config_parallel.sh's max_parallel). If you're ",
      "running on different hardware where this may no longer hold, re-benchmark ",
      "before removing this guard.",
      call. = FALSE
    )
  }
  RcppParallel::setThreadOptions(numThreads=n_threads)

  # Hardcoding n_cores = 1, needs more attention to make multi-core works with multi-threading
  n_cores <- 1
  do_type_config <- .get_do_type(n_cores, log_file)
  `%doType%` <- do_type_config$doType
  
  chain_out <- foreach::foreach(iChain=1:n_chains) %doType% {
    
    this_seed = random_seed + iChain - 1
    
    if(verbose) cat("Calling bcfoverparRcppClean From R\n")
    set.seed(this_seed)
    
    tree_files = .get_chain_tree_files(save_tree_directory, iChain, no_output)

    fitbcf = bcfoverparRcppClean(y_ = yscale[perm], z_ = z[perm], w_ = w[perm],
                                 x_con_ = t(x_c[perm,,drop=FALSE]), x_mod_ = t(x_m[perm,,drop=FALSE]), 
                                 x_con_info_list = cutpoint_list_c, 
                                 x_mod_info_list = cutpoint_list_m,
                                 random_des = matrix(1),
                                 random_var = matrix(1),
                                 random_var_ix = matrix(1),
                                 random_var_df = 3,
                                 burn = nburn, nd = nsim, thin = nthin,
                                 ntree_mod = ntree_moderate, ntree_con = ntree_control, 
                                 lambda = lambda, nu = nu,
                                 con_sd = con_sd,
                                 mod_sd = mod_sd, # if HN make sd_moderate the prior median
                                 mod_alpha = base_moderate, 
                                 mod_beta = power_moderate, 
                                 con_alpha = base_control, 
                                 con_beta = power_control,
                                 treef_con_name_ = tree_files$con_trees, 
                                 treef_mod_name_ = tree_files$mod_trees, 
                                 treef_var_name_ = tree_files$var_trees,
                                 treef_var_ratio_name_ = tree_files$var_ratio_trees,
                                 status_interval = update_interval,
                                 use_mscale = use_muscale, use_bscale = use_tauscale, 
                                 b_half_normal = TRUE, verbose_sigma=verbose, 
                                 no_output=no_output,
                                 x_var_ = t(x_v[perm,,drop=FALSE]),
                                 x_var_info_list = cutpoint_list_v,
                                 ntree_var = vartree_params$num_trees,
                                 var_lambda = vartree_params$lambda_tree,
                                 var_nu = vartree_params$nu_tree,
                                 var_alpha = vartree_params$base,
                                 var_beta = vartree_params$power,
                                 use_hetero = use_hetero,
                                 use_ratio = use_ratio)
    
    if(verbose) cat("bcfoverparRcppClean returned to R\n")

    ac = fitbcf$m_post[,order(perm)]

    Tm = fitbcf$b_post[,order(perm)] * (1.0/ (fitbcf$b1 - fitbcf$b0))

    Tc = ac * (1.0/fitbcf$msd) 

    tau_post = sdy*fitbcf$b_post[,order(perm)]

    mu_post  = muy + sdy*(Tc*fitbcf$msd + Tm*fitbcf$b0)
    sigma2_post = sdy^2 * fitbcf$sigma2_post[,order(perm)]
    sigma0_2_post = sdy^2 * fitbcf$sigma0_2_post[,order(perm)]
    sigma1_2_post = sdy^2 * fitbcf$sigma1_2_post[,order(perm)]
    log_var_ratio_post = fitbcf$log_var_ratio_post[,order(perm)]
    
    list(sigma = sdy*fitbcf$sigma,
         yhat = muy + sdy*fitbcf$yhat_post[,order(perm)],
         sigma2 = sigma2_post,
         sigma0_2 = sigma0_2_post,
         sigma1_2 = sigma1_2_post,
         log_var_ratio = log_var_ratio_post,
         sdy = sdy,
         con_sd = con_sd,
         mod_sd = mod_sd,
         muy = muy,
         mu  = mu_post,
         tau = tau_post,
         mu_scale = fitbcf$msd,
         tau_scale = fitbcf$bsd,
         b0 = fitbcf$b0,
         b1 = fitbcf$b1,
         perm = perm,
         include_pi = include_pi,
         variance_model = variance_model,
         vartree = if (use_hetero) vartree_params else NULL,
         random_seed=this_seed, 
         has_file_output=!no_output
    )

  }


  all_sigma = c()
  all_mu_scale = c()
  all_tau_scale = c()

  all_b0 = c()
  all_b1 = c()
  
  all_yhat = c()
  all_mu   = c()
  all_tau  = c()
  all_sigma2 = c()
  all_sigma0_2 = c()
  all_sigma1_2 = c()
  all_log_var_ratio = c()
  
  chain_list=list()

  n_iter = length(chain_out[[1]]$sigma)
  
  for (iChain in 1:n_chains){
    sigma            <- chain_out[[iChain]]$sigma
    mu_scale         <- chain_out[[iChain]]$mu_scale
    tau_scale        <- chain_out[[iChain]]$tau_scale
    
    b0               <- chain_out[[iChain]]$b0
    b1               <- chain_out[[iChain]]$b1

    yhat             <- chain_out[[iChain]]$yhat
    tau              <- chain_out[[iChain]]$tau
    mu               <- chain_out[[iChain]]$mu
    sigma2           <- chain_out[[iChain]]$sigma2
    sigma0_2         <- chain_out[[iChain]]$sigma0_2
    sigma1_2         <- chain_out[[iChain]]$sigma1_2
    log_var_ratio    <- chain_out[[iChain]]$log_var_ratio
    has_file_output  <- chain_out[[iChain]]$has_file_output

    # -----------------------------    
    # Support Old Output
    # -----------------------------
    all_sigma       = c(all_sigma,     sigma)
    all_mu_scale    = c(all_mu_scale,  mu_scale)
    all_tau_scale   = c(all_tau_scale, tau_scale)
    all_b0 = c(all_b0, b0)
    all_b1 = c(all_b1, b1)

    all_yhat = rbind(all_yhat, yhat)
    all_mu   = rbind(all_mu,   mu)
    all_tau  = rbind(all_tau,  tau)
    all_sigma2 = rbind(all_sigma2, sigma2)
    all_sigma0_2 = rbind(all_sigma0_2, sigma0_2)
    all_sigma1_2 = rbind(all_sigma1_2, sigma1_2)
    all_log_var_ratio = rbind(all_log_var_ratio, log_var_ratio)

    # -----------------------------    
    # Make the MCMC Object
    # -----------------------------

    scalar_df <- data.frame("sigma"     = sigma,
                            "tau_bar"   = matrixStats::rowWeightedMeans(tau, w),
                            "mu_bar"    = matrixStats::rowWeightedMeans(mu, w),
                            "yhat_bar"  = matrixStats::rowWeightedMeans(yhat, w),
                            "sigma2_bar" = matrixStats::rowWeightedMeans(sigma2, w),
                            "sigma0_2_bar" = matrixStats::rowWeightedMeans(sigma0_2, w),
                            "sigma1_2_bar" = matrixStats::rowWeightedMeans(sigma1_2, w),
                            "log_var_ratio_bar" = matrixStats::rowWeightedMeans(log_var_ratio, w),
                            "mu_scale"  = mu_scale, 
                            # "tau_scale" = tau_scale,
                            "b0"  = b0, 
                            "b1"  = b1)
    
    # y_df <- as.data.frame(chain$yhat)
    # colnames(y_df) <- paste0('y',1:ncol(y_df))
    # 
    # mu_df <- as.data.frame(chain$mu)
    # colnames(mu_df) <- paste0('mu',1:ncol(mu_df))
    # 
    # tau_df <- as.data.frame(chain$tau)
    # colnames(tau_df) <- paste0('tau',1:ncol(tau_df))
    
    chain_list[[iChain]] <- coda::as.mcmc(scalar_df)
    # -----------------------------    
    # Sanity Check Constants Across Chains
    # -----------------------------
    if(chain_out[[iChain]]$sdy              != chain_out[[1]]$sdy)              stop("sdy not consistent between chains for no reason")
    if(chain_out[[iChain]]$con_sd           != chain_out[[1]]$con_sd)           stop("con_sd not consistent between chains for no reason")
    if(chain_out[[iChain]]$mod_sd           != chain_out[[1]]$mod_sd)           stop("mod_sd not consistent between chains for no reason")
    if(chain_out[[iChain]]$muy              != chain_out[[1]]$muy)              stop("muy not consistent between chains for no reason")
    if(chain_out[[iChain]]$include_pi       != chain_out[[1]]$include_pi)       stop("include_pi not consistent between chains for no reason")
    if(chain_out[[iChain]]$variance_model   != chain_out[[1]]$variance_model)   stop("variance_model not consistent between chains for no reason")
    if(any(chain_out[[iChain]]$perm         != chain_out[[1]]$perm))            stop("perm not consistent between chains for no reason")
    if(chain_out[[iChain]]$has_file_output  != chain_out[[1]]$has_file_output)  stop("has_file_output not consistent between chains for no reason")
  }
  
  fitObj <- list(sigma = all_sigma,
                 yhat = all_yhat,
                 sdy = chain_out[[1]]$sdy,
                 muy = chain_out[[1]]$muy,
                 mu  = all_mu,
                 tau = all_tau,
                 sigma2 = all_sigma2,
                 sigma0_2 = all_sigma0_2,
                 sigma1_2 = all_sigma1_2,
                 log_var_ratio = all_log_var_ratio,
                 mu_scale = all_mu_scale,
                 tau_scale = all_tau_scale,
                 b0 = all_b0,
                 b1 = all_b1,
                 perm = perm,
                 include_pi = chain_out[[1]]$include_pi,
                 standardize = standardize,
                 variance_model = chain_out[[1]]$variance_model,
                 vartree = chain_out[[1]]$vartree,
                 random_seed = chain_out[[1]]$random_seed,
                 coda_chains = coda::as.mcmc.list(chain_list),
                 raw_chains = chain_out, 
                 has_file_output = has_file_output)
  
  attr(fitObj, "class") <- "bcf"
  
  .cleanup_after_par(do_type_config)
  
  return(fitObj)
}

#' Fit Bayesian Causal Forests With Heteroscedastic Residual Variance
#'
#' Fits a BCF model with a shared covariate-dependent residual variance or with
#' separate treatment and control variances linked by a covariate-dependent
#' variance ratio.
#'
#' When tree output is enabled with \code{no_output = FALSE}, serialized
#' variance trees are written alongside the mean-model trees for use by
#' \code{predict(..., type = "sigma2")}.
#'
#' @inheritParams bcf
#' @param ... Additional arguments passed to \code{bcf()}.
#' @param x_variance Design matrix for the residual variance function.
#' @param vartree List of scalar product-of-trees variance prior parameters.
#' @param variance_model Either \code{"shared"} or \code{"ratio"}. The ratio
#' model estimates \code{log(sigma1^2 / sigma0^2)} in addition to baseline variance.
#' @return A fitted \code{bcf_hetero} object with posterior draws \code{mu},
#' \code{tau}, and variance draws. Shared mode returns \code{sigma2}; ratio
#' mode also returns \code{sigma0_2}, \code{sigma1_2}, and \code{log_var_ratio},
#' where \code{log_var_ratio = log(sigma1_2 / sigma0_2)}.
#' @export
bcf_hetero <- function(y, z, x_control, x_moderate=x_control, pihat, ...,
                       x_variance = x_control,
                       vartree = list(num_trees = 40, nu = 10, lambda = NULL,
                                      numcut = 100, sparse = FALSE),
                       variance_model = "shared",
                       standardize = TRUE) {
  if (!variance_model %in% c("shared", "ratio")) {
    stop("variance_model must be 'shared' or 'ratio'")
  }
  fit <- bcf(y = y, z = z, x_control = x_control, x_moderate = x_moderate,
             pihat = pihat, ..., x_variance = x_variance, vartree = vartree,
             variance_model = variance_model, standardize = standardize)
  fit$sigma_mean <- fit$sigma
  fit$sigma <- sqrt(fit$sigma2)
  if (identical(variance_model, "ratio")) {
    fit$sigma0 <- sqrt(fit$sigma0_2)
    fit$sigma1 <- sqrt(fit$sigma1_2)
  }
  attr(fit, "class") <- c("bcf_hetero", "bcf")
  fit
}

#' Variance-Aware CARA Treatment Scores From Heteroscedastic BCF Draws
#'
#' Computes posterior mean treatment scores using
#' \deqn{U(m, \sigma^2; \alpha) =
#' -\exp(-\alpha m + 0.5 \alpha^2 \sigma^2).}
#' Ratio fits use the treatment-specific residual variances. The exponent is
#' clamped to \code{exponent_bounds} before exponentiation.
#'
#' @param fit A fitted object returned by \code{bcf_hetero()}.
#' @param alpha Positive CARA risk-aversion parameter.
#' @param exponent_bounds Numeric length-two bounds applied to the utility
#' exponent before exponentiation to avoid overflow.
#' @return A vector of posterior mean treatment scores, one per observation.
#' @export
bcf_hetero_cara_score <- function(fit, alpha, exponent_bounds = c(-700, 700)) {
  if (is.null(fit$mu) || is.null(fit$tau) || is.null(fit$sigma2)) {
    stop("fit must contain posterior draw matrices mu, tau, and sigma2")
  }
  if (!is.matrix(fit$mu) || !is.matrix(fit$tau) || !is.matrix(fit$sigma2)) {
    stop("fit$mu, fit$tau, and fit$sigma2 must be matrices")
  }
  if (!identical(dim(fit$mu), dim(fit$tau)) || !identical(dim(fit$mu), dim(fit$sigma2))) {
    stop("fit$mu, fit$tau, and fit$sigma2 must have identical dimensions")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha < 0) {
    stop("alpha must be a single finite nonnegative value")
  }
  if (length(exponent_bounds) != 2L || any(!is.finite(exponent_bounds)) ||
      exponent_bounds[1] >= exponent_bounds[2]) {
    stop("exponent_bounds must be a finite increasing length-two vector")
  }

  if (identical(fit$variance_model, "ratio")) {
    if (is.null(fit$sigma0_2) || is.null(fit$sigma1_2)) {
      stop("ratio fits must contain sigma0_2 and sigma1_2")
    }
    if (!identical(dim(fit$mu), dim(fit$sigma0_2)) ||
        !identical(dim(fit$mu), dim(fit$sigma1_2))) {
      stop("fit$mu, fit$sigma0_2, and fit$sigma1_2 must have identical dimensions")
    }
    treated_u <- .cara_utility(fit$mu + fit$tau, fit$sigma1_2, alpha, exponent_bounds)
    control_u <- .cara_utility(fit$mu, fit$sigma0_2, alpha, exponent_bounds)
  } else {
    treated_u <- .cara_utility(fit$mu + fit$tau, fit$sigma2, alpha, exponent_bounds)
    control_u <- .cara_utility(fit$mu, fit$sigma2, alpha, exponent_bounds)
  }
  colMeans(treated_u - control_u)
}

#' Diagnostics for Heteroscedastic BCF Simulation Comparisons
#'
#' Computes posterior residual-variance summaries, optional correlations with
#' true residual variances, variance-aware CARA scores, and rank and top-k
#' comparisons with a reference fit.
#'
#' @param fit A fitted object returned by \code{bcf_hetero()}.
#' @param fit_reference Optional fitted standard BCF object used for overlap and rank comparisons.
#' @param true_sigma2 Optional vector of true residual variances.
#' @param true_sigma0_2 Optional vector of true control residual variances.
#' @param true_sigma1_2 Optional vector of true treatment residual variances.
#' @param alpha Optional CARA risk-aversion parameter. If supplied, variance-aware scores are computed.
#' @param top_k Optional number of observations selected for overlap/Jaccard comparisons.
#' @param exponent_bounds Numeric length-two bounds used by \code{bcf_hetero_cara_score()}.
#' @return A list containing posterior summaries and any requested diagnostics.
#' @export
bcf_hetero_diagnostics <- function(fit, fit_reference = NULL, true_sigma2 = NULL,
                                   true_sigma0_2 = NULL, true_sigma1_2 = NULL,
                                   alpha = NULL, top_k = NULL,
                                   exponent_bounds = c(-700, 700)) {
  if (is.null(fit$sigma2) || !is.matrix(fit$sigma2)) {
    stop("fit must be a heteroscedastic BCF fit with a sigma2 draw matrix")
  }
  n <- ncol(fit$sigma2)
  sigma2_mean <- colMeans(fit$sigma2)

  out <- list(
    sigma2_mean = sigma2_mean,
    mu_mean = if (!is.null(fit$mu)) colMeans(fit$mu) else NULL,
    tau_mean = if (!is.null(fit$tau)) colMeans(fit$tau) else NULL
  )

  if (identical(fit$variance_model, "ratio")) {
    out$sigma0_2_mean <- colMeans(fit$sigma0_2)
    out$sigma1_2_mean <- colMeans(fit$sigma1_2)
    out$log_var_ratio_mean <- colMeans(fit$log_var_ratio)
  }

  if (!is.null(true_sigma2)) {
    if (length(true_sigma2) != n) stop("true_sigma2 must have length ncol(fit$sigma2)")
    finite <- is.finite(true_sigma2) & is.finite(sigma2_mean)
    out$sigma2_correlation <- if (sum(finite) >= 2) {
      stats::cor(true_sigma2[finite], sigma2_mean[finite])
    } else {
      NA_real_
    }
  }

  if (!is.null(true_sigma0_2)) {
    if (!identical(fit$variance_model, "ratio")) {
      stop("true_sigma0_2 requires a ratio variance fit")
    }
    if (length(true_sigma0_2) != n) stop("true_sigma0_2 must have length ncol(fit$sigma2)")
    finite <- is.finite(true_sigma0_2) & is.finite(out$sigma0_2_mean)
    out$sigma0_2_correlation <- if (sum(finite) >= 2) {
      stats::cor(true_sigma0_2[finite], out$sigma0_2_mean[finite])
    } else {
      NA_real_
    }
  }

  if (!is.null(true_sigma1_2)) {
    if (!identical(fit$variance_model, "ratio")) {
      stop("true_sigma1_2 requires a ratio variance fit")
    }
    if (length(true_sigma1_2) != n) stop("true_sigma1_2 must have length ncol(fit$sigma2)")
    finite <- is.finite(true_sigma1_2) & is.finite(out$sigma1_2_mean)
    out$sigma1_2_correlation <- if (sum(finite) >= 2) {
      stats::cor(true_sigma1_2[finite], out$sigma1_2_mean[finite])
    } else {
      NA_real_
    }
  }

  if (!is.null(alpha)) {
    out$cara_score <- bcf_hetero_cara_score(fit, alpha, exponent_bounds)
  }

  if (!is.null(fit_reference)) {
    if (is.null(fit_reference$tau) || !is.matrix(fit_reference$tau)) {
      stop("fit_reference must contain a tau draw matrix")
    }
    if (ncol(fit_reference$tau) != n) {
      stop("fit_reference$tau must have the same number of columns as fit$sigma2")
    }
    hetero_rank_score <- if (!is.null(out$cara_score)) out$cara_score else colMeans(fit$tau)
    reference_rank_score <- colMeans(fit_reference$tau)
    out$rank_correlation <- stats::cor(hetero_rank_score, reference_rank_score,
                                       use = "pairwise.complete.obs")
    out$rank_change <- rank(-hetero_rank_score, ties.method = "average") -
      rank(-reference_rank_score, ties.method = "average")

    if (!is.null(top_k)) {
      top_k <- as.integer(top_k)
      if (!is.finite(top_k) || top_k < 1L || top_k > n) {
        stop("top_k must be an integer between 1 and ncol(fit$sigma2)")
      }
      hetero_top <- order(hetero_rank_score, decreasing = TRUE)[seq_len(top_k)]
      reference_top <- order(reference_rank_score, decreasing = TRUE)[seq_len(top_k)]
      intersection_n <- length(intersect(hetero_top, reference_top))
      union_n <- length(union(hetero_top, reference_top))
      out$top_k <- top_k
      out$top_k_overlap <- intersection_n
      out$top_k_jaccard <- intersection_n / union_n
      out$top_k_hetero <- hetero_top
      out$top_k_reference <- reference_top
    }
  }

  out
}

#' Takes a fitted bcf object produced by bcf() and produces summary stats and MCMC diagnostics.
#' This function is built using the coda package and meant to mimic output from rstan::print.stanfit().
#' It includes, for key parameters, posterior summary stats, effective sample sizes, 
#' and Gelman and Rubin's convergence diagnostics. 
#' By default, those parameters are: sigma (the error standard deviation when the weights
#' are all equal), tau_bar (the estimated sample average treatment effect), mu_bar
#' (the average outcome under control/z=0 across all observations in the sample), and
#' yhat_bat (the average outcome under the realized treatment assignment across all
#' observations in the sample).
#' 
#' We strongly suggest updating the coda package to our 
#' Github version, which uses the Stan effective size computation. 
#' We found the native coda effective size computation to be overly optimistic in some situations
#' and are in discussions with the coda package authors to change it on CRAN.
#' @param object output from a BCF predict run.
#' @param ... additional arguments affecting the summary produced.
#' @param params_2_summarise parameters to summarise.
#' @return No return value, called for side effects
#' @examples
#'\dontrun{
#'
#' # data generating process
#' p = 3 #two control variables and one moderator
#' n = 250
#' 
#' set.seed(1)
#'
#' x = matrix(rnorm(n*p), nrow=n)
#'
#' # create targeted selection
#' q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
#'
#' # generate treatment variable
#' pi = pnorm(q)
#' z = rbinom(n,1,pi)
#'
#' # tau is the true (homogeneous) treatment effect
#' tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
#'
#' # generate the response using q, tau and z
#' mu = (q + tau*z)
#'
#' # set the noise level relative to the expected mean function of Y
#' sigma = diff(range(q + tau*pi))/8
#'
#' # draw the response variable with additive error
#' y = mu + sigma*rnorm(n)
#'
#' # If you didn't know pi, you would estimate it here
#' pihat = pnorm(q)
#'
#' bcf_fit = bcf(y, z, x, x, pihat, nburn=2000, nsim=2000)
#'
#' # Get model fit diagnostics
#' summary(bcf_fit)
#'
#'}
#' @export
summary.bcf <- function(object,
                        ..., 
                        params_2_summarise = c('sigma','tau_bar','mu_bar','yhat_bar')){

  chains_2_summarise <- object$coda_chains[,params_2_summarise]

  message("Summary statistics for each Markov Chain Monte Carlo run")
  print(summary(chains_2_summarise))

  message("\n----\n\n")


  message("Effective sample size for summary parameters")
  
  ef = function(e) {
    if(e$message == "unused argument (crosschain = TRUE)") {
      message("Reverting to coda's default ESS calculation. See ?summary.bcf for details.\n\n")
      print(coda::effectiveSize(chains_2_summarise))
    } else {
      stop(e)
    }
  }
  tryCatch(print(coda::effectiveSize(chains_2_summarise, crosschain = TRUE)),
           error = ef) 
  message("\n----\n\n")
  
  
  if (length(chains_2_summarise) > 1){
    message("Gelman and Rubin's convergence diagnostic for summary parameters")
    print(coda::gelman.diag(chains_2_summarise, autoburnin = FALSE))
    message("\n----\n\n")
  }
  
}


#' Takes a fitted bcf object produced by bcf() along with serialized tree samples and produces predictions for a new set of covariate values
#' 
#' This function takes in an existing BCF model fit and uses it to predict estimates for new data.
#' It is important to note that this function requires that you indicate where the trees from the model fit are saved.
#' You can do so using the save_tree_directory argument in bcf(). Otherwise, they will be saved in the working directory.
#' @param object output from a BCF predict run
#' @param ... additional arguments affecting the predictions produced.
#' @param x_predict_control matrix of covariates for the "prognostic" function mu(x) for predictions (optional)
#' @param x_predict_moderate matrix of covariates for the covariate-dependent treatment effects tau(x) for predictions (optional)
#' @param z_pred Treatment variable for predictions (optional except if x_pre is not empty)
#' @param pi_pred propensity score for prediction
#' @param save_tree_directory directory where the trees have been saved
#' @param x_predict_variance matrix of covariates for the residual variance function in heteroscedastic fits
#' @param type Which posterior prediction components to return
#' @param log_file File to log progress
#' @param n_cores An optional integer of the number of cores to run your MCMC chains on
#' @param verbose Logical; set to FALSE to suppress extra output
#' @return A list with elements: tau (samples of treatment effects), mu
#' (samples of predicted control outcomes), yhat (samples of predicted values),
#' and coda_chains (coda objects for scalar summaries). Heteroscedastic fits
#' also return sigma2 and sigma draws when type is "all" or "sigma2"; ratio
#' fits additionally return sigma0_2, sigma1_2, and log_var_ratio, where
#' log_var_ratio is log(sigma1_2 / sigma0_2). The log standard-deviation ratio
#' is log_var_ratio / 2.
#' @examples
#'\dontrun{
#'
#' # data generating process
#' p = 3 #two control variables and one moderator
#' n = 250
#'
#' x = matrix(rnorm(n*p), nrow=n)
#'
#' # create targeted selection
#' q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
#'
#' # generate treatment variable
#' pi = pnorm(q)
#' z = rbinom(n,1,pi)
#'
#' # tau is the true (homogeneous) treatment effect
#' tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
#'
#' # generate the response using q, tau and z
#' mu = (q + tau*z)
#'
#' # set the noise level relative to the expected mean function of Y
#' sigma = diff(range(q + tau*pi))/8
#'
#' # draw the response variable with additive error
#' y = mu + sigma*rnorm(n)
#'
#' # If you didn't know pi, you would estimate it here
#' pihat = pnorm(q)
#'
#' n_burn = 5000
#' n_sim = 5000
#'
#' bcf_fit = bcf(y               = y,
#'               z               = z,
#'               x_control       = x,
#'               x_moderate      = x,
#'               pihat           = pihat,
#'               nburn           = n_burn,
#'               nsim            = n_sim,
#'               n_chains        = 2,
#'               update_interval = 100,
#'               save_tree_directory = './trees')
#'
#' # Predict using new data
#' 
#' x_pred = matrix(rnorm(n*p), nrow=n)
#' 
#' pred_out = predict(bcf_out=bcf_fit,
#'                    x_predict_control=x_pred,
#'                    x_predict_moderate=x_pred,
#'                    pi_pred=pihat,
#'                    z_pred=z,
#'                    save_tree_directory = './trees')
#'
#'}
#' @export
predict.bcf <- function(object, 
                        x_predict_control,
                        x_predict_moderate,
                        pi_pred,
                        z_pred, 
                        save_tree_directory,
                        x_predict_variance = x_predict_control,
                        type = c("all", "mu", "tau", "sigma2"),
                        log_file=file.path('.',sprintf('bcf_log_%s.txt',format(Sys.time(), "%Y%m%d_%H%M%S"))),
                        n_cores=2, verbose = TRUE,
                        ...) {
    type <- match.arg(type)
                        
    if(any(is.na(x_predict_moderate))) stop("Missing values in x_predict_moderate")
    if(any(is.na(x_predict_control))) stop("Missing values in x_predict_control")
    if(any(is.na(x_predict_variance))) stop("Missing values in x_predict_variance")
    if(any(is.na(z_pred))) stop("Missing values in z_pred")
    if(any(!is.finite(x_predict_moderate))) stop("Non-numeric values in x_pred_moderate")
    if(any(!is.finite(x_predict_control))) stop("Non-numeric values in x_pred_control")
    if(any(!is.finite(x_predict_variance))) stop("Non-numeric values in x_predict_variance")
    if(any(!is.finite(pi_pred))) stop("Non-numeric values in pi_pred")
    if(!all(sort(unique(z_pred)) == c(0,1))) stop("z_pred must be a vector of 0's and 1's, with at least one of each")
    if(!object$has_file_output) stop("No tree samples were serialized during sampling. To enable prediction, re-run bcf with no_output = FALSE \n")

    if((is.null(x_predict_moderate) & !is.null(x_predict_control)) | (!is.null(x_predict_moderate) & is.null(x_predict_control))) {
        stop("If you want to predict, you need to add values to both x_pred_control and x_pred_moderate")
    }

    pi_pred = as.matrix(pi_pred)
    if(!.ident(length(z_pred),
                nrow(x_predict_moderate),
                nrow(x_predict_control),
                nrow(x_predict_variance),
                nrow(pi_pred))
        ) {
        stop("Data size mismatch. The following should all be equal:
            length(z_pred): ", length(z_pred), "\n",
            "nrow(x_pred_moderate): ", nrow(x_predict_moderate), "\n",
            "nrow(x_pred_control): ", nrow(x_predict_control), "\n",
            "nrow(x_predict_variance): ", nrow(x_predict_variance), "\n",
            "nrow(pi_pred): ", nrow(pi_pred), "\n"
        )
    }

    message("Initializing BCF Prediction\n")
    x_pm = matrix(x_predict_moderate, ncol=ncol(x_predict_moderate))
    x_pc = matrix(x_predict_control, ncol=ncol(x_predict_control))
    x_pv = matrix(x_predict_variance, ncol=ncol(x_predict_variance))

    if(object$include_pi=="both" | object$include_pi=="control") {
        x_pc = cbind(x_predict_control, pi_pred)
    }
    if(object$include_pi=="both" | object$include_pi=="moderate") {
        x_pm = cbind(x_predict_moderate, pi_pred)
    }


    message("Starting Prediction \n")

    n_chains = length(object$coda_chains)
    muy = object$muy
    sdy = object$sdy
    
    do_type_config <- .get_do_type(n_cores, log_file=log_file)
    `%doType%` <- do_type_config$doType
    
    chain_out <- foreach::foreach(iChain=1:n_chains) %doType% {
      
      tree_files = .get_chain_tree_files(save_tree_directory, iChain)

      if(verbose) cat("Starting to Predict Chain ", iChain, "\n")
      
      mods = TreeSamples$new()
      mods$load(tree_files$mod_trees)
      Tm = mods$predict(t(x_pm))
      
      cons = TreeSamples$new()
      cons$load(tree_files$con_trees)
      Tc = cons$predict(t(x_pc))

      Sigma2 = NULL
      Sigma0_2 = NULL
      Sigma1_2 = NULL
      LogVarRatio = NULL
      if(identical(object$variance_model, "shared")) {
        vars = TreeSamples$new()
        vars$load(tree_files$var_trees)
        Sigma2 = sdy^2 * vars$predict_prec(t(x_pv))
      }
      if(identical(object$variance_model, "ratio")) {
        vars = TreeSamples$new()
        vars$load(tree_files$var_trees)
        vars_ratio = TreeSamples$new()
        vars_ratio$load(tree_files$var_ratio_trees)
        Sigma0_2 = sdy^2 * vars$predict_prec(t(x_pv))
        SigmaRatio = vars_ratio$predict_prec(t(x_pv))
        Sigma1_2 = Sigma0_2 * SigmaRatio
        LogVarRatio = log(SigmaRatio)
        Sigma2 = Sigma0_2 * ifelse(z_pred == 1, SigmaRatio, 1)
      }
      
      
      list(Tm = Tm,
           Tc = Tc,
           Sigma2 = Sigma2,
           Sigma0_2 = Sigma0_2,
           Sigma1_2 = Sigma1_2,
           LogVarRatio = LogVarRatio)
    }
    
    all_yhat = c()
    all_mu   = c()
    all_tau  = c()
    all_sigma2 = c()
    all_sigma0_2 = c()
    all_sigma1_2 = c()
    all_log_var_ratio = c()
    
    chain_list=list()

    for (iChain in 1:n_chains){
      
      
        # Extract Chain Specific Information
    
        Tm = chain_out[[iChain]]$Tm
        Tc = chain_out[[iChain]]$Tc
        
        this_chain_bcf_out = object$raw_chains[[iChain]]
        
        b1 = this_chain_bcf_out$b1
        b0 = this_chain_bcf_out$b0
        mu_scale = this_chain_bcf_out$mu_scale
        


        # Calculate, tau, y, and mu

        
        mu  = muy + sdy*(Tc*mu_scale + Tm*b0)
        tau = sdy*(b1 - b0)*Tm
        yhat = mu + t(t(tau)*z_pred)
        
        
        # Package Output up
        all_yhat = rbind(all_yhat, yhat)
        all_mu   = rbind(all_mu,   mu)
        all_tau  = rbind(all_tau,  tau)
        if(identical(object$variance_model, "shared")) {
          all_sigma2 = rbind(all_sigma2, chain_out[[iChain]]$Sigma2)
        }
        if(identical(object$variance_model, "ratio")) {
          all_sigma2 = rbind(all_sigma2, chain_out[[iChain]]$Sigma2)
          all_sigma0_2 = rbind(all_sigma0_2, chain_out[[iChain]]$Sigma0_2)
          all_sigma1_2 = rbind(all_sigma1_2, chain_out[[iChain]]$Sigma1_2)
          all_log_var_ratio = rbind(all_log_var_ratio, chain_out[[iChain]]$LogVarRatio)
        }
        
        
        
        scalar_df <- data.frame("tau_bar"   = matrixStats::rowWeightedMeans(tau, w=NULL),
                                "mu_bar"    = matrixStats::rowWeightedMeans(mu, w=NULL),
                                "yhat_bar"  = matrixStats::rowWeightedMeans(yhat, w=NULL))

        chain_list[[iChain]] <- coda::as.mcmc(scalar_df)
    }

   .cleanup_after_par(do_type_config)


    out <- list(tau = all_tau,
                mu = all_mu,
                yhat = all_yhat,
                coda_chains = coda::as.mcmc.list(chain_list))
    if(identical(object$variance_model, "shared")) {
      out$sigma2 <- all_sigma2
      out$sigma <- sqrt(all_sigma2)
    }
    if(identical(object$variance_model, "ratio")) {
      out$sigma2 <- all_sigma2
      out$sigma <- sqrt(all_sigma2)
      out$sigma0_2 <- all_sigma0_2
      out$sigma1_2 <- all_sigma1_2
      out$sigma0 <- sqrt(all_sigma0_2)
      out$sigma1 <- sqrt(all_sigma1_2)
      out$log_var_ratio <- all_log_var_ratio
    }
    if(type == "mu") return(out["mu"])
    if(type == "tau") return(out["tau"])
    if(type == "sigma2") {
      if(!object$variance_model %in% c("shared", "ratio")) {
        stop("type = 'sigma2' requires a heteroscedastic BCF fit")
      }
      if(identical(object$variance_model, "ratio")) {
        return(out[c("sigma2", "sigma", "sigma0_2", "sigma1_2",
                     "sigma0", "sigma1", "log_var_ratio")])
      }
      return(out[c("sigma2", "sigma")])
    }
    out
}
