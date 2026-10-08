#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <algorithm>
#include <cstdio>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"
#include "logging.h"
#include "varfuns.h"
#include "paired_mean.h"
#include "joint_variance.h"
#include "paired_variance.h"

using namespace Rcpp;
// Rstudios check's suggest not ignoring these
// #pragma GCC diagnostic ignored "-Wunused-parameter"
// #pragma GCC diagnostic ignored "-Wcomment"
// #pragma GCC diagnostic ignored "-Wformat"
// #pragma GCC diagnostic ignored "-Wsign-compare"

// y = m(x) + b(x)z + e, e~N(0, sigma^2_y

//x_con is the design matrix for m. It should have n = rows
//x_mod is the design matrix for b. It should have n = rows
//data should come in sorted with all trt first, then control cases

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List bcfoverparRcppClean(NumericVector y_, NumericVector z_, NumericVector w_,
                  NumericVector x_con_, NumericVector x_mod_,
                  List x_con_info_list, List x_mod_info_list,
                  arma::mat random_des, //needs to come in with n rows no matter what(?)
                  arma::mat random_var, arma::mat random_var_ix, //random_var_ix*random_var = diag(Var(random effects))
                  double random_var_df,
                  int burn, int nd, int thin, //Draw nd*thin + burn samples, saving nd draws after burn-in
                  int ntree_mod, int ntree_con,
                  double lambda, double nu, //prior pars for sigma^2_y
                  double con_sd, // Var(m(x)) = con_sd^2 marginally a priori (approx)
                  double mod_sd, // Var(b(x)) = mod_sd^2 marginally a priori (approx)
                  double con_alpha, double con_beta,
                  double mod_alpha, double mod_beta,
                  CharacterVector treef_con_name_, CharacterVector treef_mod_name_,
                  CharacterVector treef_var_name_,
                  CharacterVector treef_var_ratio_name_,
                  int status_interval=100,
                  bool RJ= false, bool use_mscale=true, bool use_bscale=true, bool b_half_normal=true,
                  double trt_init = 1.0, bool verbose_sigma=false, 
                  bool no_output=false,
                  NumericVector x_var_ = NumericVector::create(),
                  List x_var_info_list = List::create(),
                  int ntree_var = 0,
                  double var_lambda = 1.0,
                  double var_nu = 10.0,
                  double var_alpha = 0.95,
                  double var_beta = 2.0,
                  bool use_hetero = false,
                  bool use_ratio = false,
                  bool use_paired_mean_update = false, bool use_global_mean_update = false, int joint_mean_every = 0, bool collapsed_mu_scale = false, bool use_mean_split_change = false, bool use_variance_split_change = false, int joint_variance_every = 0, int paired_variance_every = 0)
{

  if((use_paired_mean_update || use_global_mean_update) && !use_hetero) stop("Paired mean update requires heteroskedastic likelihood");
  if(use_variance_split_change && !use_hetero) stop("Variance split changes require heteroskedastic likelihood");
  if(joint_mean_every<0) stop("joint_mean_every must be nonnegative");
  if(joint_mean_every>0 && !use_hetero) stop("Joint mean refresh currently requires heteroskedastic likelihood");
  if(joint_variance_every<0 || (joint_variance_every>0 && (!use_hetero || !use_ratio))) stop("Joint variance update requires ratio likelihood and nonnegative interval");
  if(paired_variance_every<0 || (paired_variance_every>0 && (!use_hetero || !use_ratio))) stop("Paired variance update requires ratio model and nonnegative interval");
  double paired_variance_attempts=0,paired_variance_accepts=0;
  double joint_variance_attempts=0,joint_variance_accepts=0;
  double joint_mean_updates=0, joint_mean_max_leaves=0;
  double mean_split_change_attempts=0, mean_split_change_accepts=0;
  double variance_split_change_attempts=0, variance_split_change_accepts=0;
  double paired_mean_updates=0, paired_mean_skips=0;
  bool randeff = true;
  if(random_var_ix.n_elem == 1) {
    randeff = false;
  }

  if(randeff) Rcout << "Using random effects." << std::endl;

  std::ofstream treef_con;
  std::ofstream treef_mod;
  std::ofstream treef_var;
  std::ofstream treef_var_ratio;

  std::string treef_con_name = as<std::string>(treef_con_name_);
  std::string treef_mod_name = as<std::string>(treef_mod_name_);
  std::string treef_var_name = as<std::string>(treef_var_name_);
  std::string treef_var_ratio_name = as<std::string>(treef_var_ratio_name_);

  if((not treef_con_name.empty()) && (not no_output)){
    Rcout << "Saving Trees to"  << std::endl;
    Rcout << treef_con_name  << std::endl;
    Rcout << treef_mod_name  << std::endl;

    treef_con.open(treef_con_name.c_str());
    treef_mod.open(treef_mod_name.c_str());
    if(use_hetero && !treef_var_name.empty()) {
      Rcout << treef_var_name << std::endl;
      treef_var.open(treef_var_name.c_str());
    }
    if(use_ratio && !treef_var_ratio_name.empty()) {
      Rcout << treef_var_ratio_name << std::endl;
      treef_var_ratio.open(treef_var_ratio_name.c_str());
    }
  } else {  
    Rcout << "Not Saving Trees to file"  << std::endl;
  }

  RNGScope scope;
  RNG gen; //this one random number generator is used in all draws

  //double lambda = 1.0; //this one really needs to be set
  //double nu = 3.0;
  //double kfac=2.0; //original is 2.0

  Logger logger = Logger();
  char logBuff[100];

  bool log_level = false;

  logger.setLevel(log_level);
  logger.log("============================================================");
  logger.log(" Starting up BCF: ");
  logger.log("============================================================");
  if (log_level){
    logger.getVectorHead(y_, logBuff);
    Rcout << "y: " <<  logBuff << "\n";
    logger.getVectorHead(z_, logBuff);
    Rcout << "z: " <<  logBuff << "\n";
    logger.getVectorHead(w_, logBuff);
    Rcout << "w: " <<  logBuff << "\n";
  }


  logger.log("BCF is Weighted");

  // Rprintf(logBuff, "Updating Moderate Tree: %d of %d");
  // logger.log(logBuff);
  logger.log("");

  /*****************************************************************************
  /* Read, format y
  *****************************************************************************/
  std::vector<double> y; //storage for y
  double miny = INFINITY, maxy = -INFINITY;
  sinfo allys;       //sufficient stats for all of y, use to initialize the bart trees.
  double allys_y2 = 0;

  for(NumericVector::iterator it=y_.begin(); it!=y_.end(); ++it) {
    y.push_back(*it);
    if(*it<miny) miny=*it;
    if(*it>maxy) maxy=*it;
    allys.sy += *it; // sum of y
    allys_y2 += (*it)*(*it); // sum of y^2
  }
  size_t n = y.size();
  allys.n = n;

  double ybar = allys.sy/n; //sample mean
  double shat = sqrt((allys_y2-n*ybar*ybar)/(n-1)); //sample standard deviation
  /*****************************************************************************
  /* Read, format  weights
  *****************************************************************************/
  double* w = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp


  for(int j=0; j<n; j++) {
    w[j] = w_[j];
  }


  /*****************************************************************************
  /* Read, format X_con
  *****************************************************************************/
  //the n*p numbers for x are stored as the p for first obs, then p for second, and so on.
  std::vector<double> x_con;
  for(NumericVector::iterator it=x_con_.begin(); it!= x_con_.end(); ++it) {
    x_con.push_back(*it);
  }
  size_t p_con = x_con.size()/n;

  Rcout << "Using " << p_con << " control variables." << std::endl;

  //x cutpoints
  xinfo xi_con;

  xi_con.resize(p_con);
  for(int i=0; i<p_con; ++i) {
    NumericVector tmp = x_con_info_list[i];
    std::vector<double> tmp2;
    for(size_t j=0; j<tmp.size(); ++j) {
      tmp2.push_back(tmp[j]);
    }
    xi_con[i] = tmp2;
  }

  /*****************************************************************************
  /* Read, format X_mod
  *****************************************************************************/
  int ntrt = 0;
  for(size_t i=0; i<n; ++i) {
    if(z_[i]>0) ntrt += 1;
  }

  /*****************************************************************************
  /* Read, format X_var for scalar residual variance d(x)
  *****************************************************************************/
  std::vector<double> x_var;
  size_t p_var = 0;
  xinfo xi_var;
  if(use_hetero) {
    for(NumericVector::iterator it=x_var_.begin(); it!= x_var_.end(); ++it) {
      x_var.push_back(*it);
    }
    if(x_var.size() == 0) stop("x_var must be supplied when use_hetero is TRUE");
    p_var = x_var.size()/n;
    if(p_var * n != x_var.size()) stop("x_var has incompatible dimensions");
    if(ntree_var <= 0) stop("ntree_var must be positive when use_hetero is TRUE");
    if(var_lambda <= 0.0) stop("var_lambda must be positive");
    if(var_nu <= 0.0) stop("var_nu must be positive");

    Rcout << "Using " << p_var << " residual variance covariates." << std::endl;
    xi_var.resize(p_var);
    for(size_t i=0; i<p_var; ++i) {
      NumericVector tmp = x_var_info_list[i];
      std::vector<double> tmp2;
      for(size_t j=0; j<tmp.size(); ++j) tmp2.push_back(tmp[j]);
      xi_var[i] = tmp2;
    }
  }
  std::vector<double> x_mod;
  for(NumericVector::iterator it=x_mod_.begin(); it!= x_mod_.end(); ++it) {
    x_mod.push_back(*it);
  }
  size_t p_mod = x_mod.size()/n;

  Rcout << "Using " << p_mod << " potential effect moderators." << std::endl;

  //x cutpoints
  xinfo xi_mod;

  xi_mod.resize(p_mod);
  for(int i=0; i<p_mod; ++i) {
    NumericVector tmp = x_mod_info_list[i];
    std::vector<double> tmp2;
    for(size_t j=0; j<tmp.size(); ++j) {
      tmp2.push_back(tmp[j]);
    }
    xi_mod[i] = tmp2;
  }

  //  Rcout <<"\nburn,nd,number of trees: " << burn << ", " << nd << ", " << m << endl;
  //  Rcout <<"\nlambda,nu,kfac: " << lambda << ", " << nu << ", " << kfac << endl;

  /*****************************************************************************
  /* Setup the model
  *****************************************************************************/
  //--------------------------------------------------
  //trees
  std::vector<tree> t_mod(ntree_mod);
  for(size_t i=0;i<ntree_mod;i++) t_mod[i].setm(trt_init/(double)ntree_mod);

  std::vector<tree> t_con(ntree_con);
  for(size_t i=0;i<ntree_con;i++) t_con[i].setm(ybar/(double)ntree_con);

  std::vector<tree> t_var(ntree_var);
  if(use_hetero) {
    double var_leaf_init = pow(var_lambda, 1.0 / (double) ntree_var);
    for(size_t i=0;i<(size_t) ntree_var;i++) t_var[i].setm(var_leaf_init);
  }
  std::vector<tree> t_var_ratio(ntree_var);
  if(use_ratio) {
    for(size_t i=0;i<(size_t) ntree_var;i++) t_var_ratio[i].setm(1.0);
  }

  //--------------------------------------------------
  //prior parameters
  // PX scale parameter for b:
  double bscale_prec = 2;
  double bscale0 = -0.5;
  double bscale1 = 0.5;

  double mscale_prec = 1.0;
  double mscale = 1.0;
  double delta_con = 1.0;
  double delta_mod = 1.0;

  pinfo pi_mod;
  pi_mod.pbd = 1.0; //prob of birth/death move
  pi_mod.pb = .5; //prob of birth given  birth/death

  pi_mod.alpha = mod_alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
  pi_mod.beta  = mod_beta;  //2 for bart means it is harder to build big trees.
  pi_mod.tau   = mod_sd/(sqrt(delta_mod)*sqrt((double) ntree_mod)); //sigma_mu, variance on leaf parameters
  pi_mod.sigma = shat; //resid variance is \sigma^2_y/bscale^2 in the backfitting update

  pinfo pi_con;
  pi_con.pbd = 1.0; //prob of birth/death move
  pi_con.pb = .5; //prob of birth given  birth/death

  pi_con.alpha = con_alpha;
  pi_con.beta  = con_beta;
  pi_con.tau   = con_sd/(sqrt(delta_con)*sqrt((double) ntree_con)); //sigma_mu, variance on leaf parameters

  pi_con.sigma = shat/fabs(mscale); //resid variance in backfitting is \sigma^2_y/mscale^2

  double sigma = shat;

  pinfo pi_var;
  if(use_hetero) {
    pi_var.pbd = 1.0;
    pi_var.pb = 0.5;
    pi_var.alpha = var_alpha;
    pi_var.beta = var_beta;
    pi_var.tau = 1.0;
    pi_var.sigma = 1.0;
  }

  // @Peter This is where dinfo is initialized

  //--------------------------------------------------
  //dinfo for control function m(x)
//  Rcout << "ybar " << ybar << endl;
  double* allfit_con = new double[n]; //sum of fit of all trees
  for(size_t i=0;i<n;i++) allfit_con[i] = ybar;
  double* r_con = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
  dinfo di_con;
  di_con.n=n;
  di_con.p = p_con;
  di_con.x = &x_con[0];
  di_con.y = r_con; //the y for each draw will be the residual

  //--------------------------------------------------
  //dinfo for trt effect function b(x)
  double* allfit_mod = new double[n]; //sum of fit of all trees
  for(size_t i=0;i<n;i++) allfit_mod[i] = (z_[i]*bscale1 + (1-z_[i])*bscale0)*trt_init;
  double* r_mod = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
  dinfo di_mod;
  di_mod.n=n;
  di_mod.p=p_mod;
  di_mod.x = &x_mod[0];
  di_mod.y = r_mod; //the y for each draw will be the residual

  double* sigma2_fit = new double[n];
  double* sigma0_2_fit = new double[n];
  double* sigma_ratio_fit = new double[n];
  double* r_var = new double[n];
  double* r_var_ratio = new double[ntrt];
  double* ftemp_var = new double[n];
  double* ftemp_var_ratio = new double[n];
  dinfo di_var;
  dinfo di_var_ratio;
  dinfo di_var_ratio_all;
  if(use_hetero) {
    for(size_t i=0;i<n;i++) {
      sigma0_2_fit[i] = var_lambda;
      sigma_ratio_fit[i] = 1.0;
      sigma2_fit[i] = var_lambda;
    }
    di_var.n = n;
    di_var.p = p_var;
    di_var.x = &x_var[0];
    di_var.y = r_var;
    if(use_ratio) {
      di_var_ratio.n = ntrt;
      di_var_ratio.p = p_var;
      di_var_ratio.x = &x_var[0];
      di_var_ratio.y = r_var_ratio;
      di_var_ratio_all.n = n;
      di_var_ratio_all.p = p_var;
      di_var_ratio_all.x = &x_var[0];
      di_var_ratio_all.y = r_var;
    }
  } else {
    for(size_t i=0;i<n;i++) {
      sigma0_2_fit[i] = sigma * sigma;
      sigma_ratio_fit[i] = 1.0;
      sigma2_fit[i] = sigma * sigma;
    }
  }

  //--------------------------------------------------
  //setup for random effects
  size_t random_dim = random_des.n_cols;
  int nr=1;
  if(randeff) nr = n;

  arma::vec r(nr); //working residuals
  arma::vec Wtr(random_dim); // W'r

  arma::mat WtW = random_des.t()*random_des; //W'W
  arma::mat Sigma_inv_random = diagmat(1/(random_var_ix*random_var));

  // PX parameters
  arma::vec eta(random_var_ix.n_cols); //random_var_ix is num random effects by num variance components
  eta.fill(1.0);

  for(size_t k=0; k<nr; ++k) {
    r(k) = y[k] - allfit_con[k] - allfit_mod[k];
  }

  Wtr = random_des.t()*r;
  arma::vec gamma = solve(WtW/(sigma*sigma)+Sigma_inv_random, Wtr/(sigma*sigma));
  arma::vec allfit_random = random_des*gamma;
  if(!randeff) allfit_random.fill(0);

  //--------------------------------------------------
  //storage for the fits
  double* allfit = new double[n]; //yhat
  for(size_t i=0;i<n;i++) {
    allfit[i] = allfit_mod[i] + allfit_con[i];
    if(randeff) allfit[i] += allfit_random[i];
  }
  double* ftemp  = new double[n]; //fit of current tree

  NumericVector sigma_post(nd);
  NumericVector msd_post(nd);
  NumericVector bsd_post(nd);
  NumericVector b0_post(nd);
  NumericVector b1_post(nd);
  NumericMatrix m_post(nd,n);
  NumericMatrix yhat_post(nd,n);
  NumericMatrix b_post(nd,n);
  NumericMatrix sigma2_post(nd,n);
  NumericMatrix sigma0_2_post(nd,n);
  NumericMatrix sigma1_2_post(nd,n);
  NumericMatrix log_var_ratio_post(nd,n);
  arma::mat gamma_post(nd,gamma.n_elem);
  arma::mat random_var_post(nd,random_var.n_elem);

  //  NumericMatrix spred2(nd,dip.n);


  // The default output precision is of C++ is 5 or 6 dp, depending on compiler.
  // I don't have much justification for 32, but it seems like a sensible number
  int save_tree_precision = 32;

  //save stuff to tree file
  // NB: use "\n" here, not std::endl -- endl force-flushes the stream on every
  // insertion, and with num_trees x num_draws lines per file that turns into
  // millions of synchronous writes (dominates wall time as kernel/sys time).
  // Content on disk is identical either way; buffering is flushed on close().
  if(not treef_con_name.empty()){
    treef_con << std::setprecision(save_tree_precision) << xi_con << "\n"; //cutpoints
    treef_con << ntree_con << "\n";  //number of trees
    treef_con << di_con.p << "\n";  //dimension of x's
    treef_con << nd << "\n";

    treef_mod << std::setprecision(save_tree_precision) << xi_mod << "\n"; //cutpoints
    treef_mod << ntree_mod << "\n";  //number of trees
    treef_mod << di_mod.p << "\n";  //dimension of x's
    treef_mod << nd << "\n";

    if(use_hetero && !treef_var_name.empty()) {
      treef_var << std::setprecision(save_tree_precision) << xi_var << "\n";
      treef_var << ntree_var << "\n";
      treef_var << di_var.p << "\n";
      treef_var << nd << "\n";
    }
    if(use_ratio && !treef_var_ratio_name.empty()) {
      treef_var_ratio << std::setprecision(save_tree_precision) << xi_var << "\n";
      treef_var_ratio << ntree_var << "\n";
      treef_var_ratio << di_var.p << "\n";
      treef_var_ratio << nd << "\n";
    }
  }

  //*****************************************************************************
  /* MCMC
   * note: the allfit objects are all carrying the appropriate scales
   */
  //*****************************************************************************
  Rcout << "\n============================================================\nBeginning MCMC:\n============================================================\n";
  time_t tp;
  int time1 = time(&tp);

  size_t save_ctr = 0;
  bool verbose_itr = false;


  double* weight      = new double[n];
  double* weight_het  = new double[n];

  logger.setLevel(0);

  bool printTrees = false;

  for(size_t iIter=0;iIter<(nd*thin+burn);iIter++) {
    // verbose_itr = iIter>=burn;
    verbose_itr = false;

    if(verbose_sigma){
        if(iIter%status_interval==0) {
            Rcout << "iteration: " << iIter << " sigma/SD(y): "<< sigma << endl;
        }
    }

    logger.setLevel(verbose_itr);

    logger.log("==============================================");
    std::snprintf(logBuff, sizeof(logBuff), "MCMC iteration: %d of %d Start", (int) iIter + 1, nd*thin+burn);
    logger.log(logBuff);
    std::snprintf(logBuff, sizeof(logBuff), "sigma %f, mscale %f, bscale0 %f, bscale1 %f",sigma, mscale, bscale0, bscale1);
    logger.log(logBuff);
    logger.log("==============================================");
    if (verbose_itr){
      logger.getVectorHead(y, logBuff);
      Rcout << "           y: " <<  logBuff << "\n";

      logger.getVectorHead(allfit, logBuff);
      Rcout << "Current Fit : " <<  logBuff << "\n";

      logger.getVectorHead(allfit_con, logBuff);
      Rcout << "allfit_con  : " <<  logBuff << "\n";

      logger.getVectorHead(allfit_mod, logBuff);
      Rcout << "allfit_mod  : " <<  logBuff << "\n";
    }

    for (int k=0; k<n; ++k){
      double obs_var = use_hetero ? sigma2_fit[k] : sigma * sigma;
      weight[k] = w[k]*mscale*mscale/obs_var; // precision for the control-tree residual
    }

    for(size_t k=0; k<ntrt; ++k) {
      double obs_var = use_hetero ? sigma2_fit[k] : sigma * sigma;
      weight_het[k] = w[k]*bscale1*bscale1/obs_var;
    }
    for(size_t k=ntrt; k<n; ++k) {
      double obs_var = use_hetero ? sigma2_fit[k] : sigma * sigma;
      weight_het[k] = w[k]*bscale0*bscale0/obs_var;
    }

    logger.log("=====================================");
    logger.log("- Tree Processing");
    logger.log("=====================================");

    //draw trees for m(x)
    for(size_t iTreeCon=0;iTreeCon<ntree_con;iTreeCon++) {

      logger.log("==================================");
      std::snprintf(logBuff, sizeof(logBuff), "Updating Control Tree: %d of %d", (int) iTreeCon + 1, ntree_con);
      logger.log(logBuff);
      logger.log("==================================");
      logger.startContext();

      logger.log("Attempting to Print Tree Pre Update \n");
      if(verbose_itr && printTrees){
        t_con[iTreeCon].pr(xi_con);
        Rcout << "\n\n";
      }

      fit(t_con[iTreeCon], // tree& t
          xi_con, // xinfo& xi
          di_con, // dinfo& di
          ftemp); // std::vector<double>& fv


      logger.log("Attempting to Print Tree Post first call to fit \n");
      if(verbose_itr && printTrees){
        t_con[iTreeCon].pr(xi_con);
        Rcout << "\n\n";
      }

      for(size_t k=0;k<n;k++) {
        if(ftemp[k] != ftemp[k]) {
          Rcout << "control tree " << iTreeCon <<" obs "<< k<<" "<< endl;
          Rcout << t_con[iTreeCon] << endl;
          stop("nan in ftemp");
        }

        allfit[k]     = allfit[k]     -mscale*ftemp[k];
        allfit_con[k] = allfit_con[k] -mscale*ftemp[k];

        r_con[k] = (y[k]-allfit[k])/mscale;

        if(r_con[k] != r_con[k]) {
          Rcout << (y[k]-allfit[k]) << endl;
          Rcout << mscale << endl;
          Rcout << r_con[k] << endl;
          stop("NaN in resid");
        }
      }



      if(verbose_itr && printTrees){
        logger.getVectorHead(weight, logBuff);
        Rcout << "\n weight: " <<  logBuff << "\n\n";
      }
      logger.log("Starting Birth / Death Processing");
      logger.startContext();
      bd(t_con[iTreeCon], // tree& x
         xi_con, // xinfo& xi
         di_con, // dinfo& di
         weight, // phi
         pi_con, // pinfo& pi
         gen,
         logger); // RNG& gen
      logger.stopContext();

      logger.log("Attempting to Print Tree Post db \n");
      if(verbose_itr && printTrees){
        t_con[iTreeCon].pr(xi_con);
        Rcout << "\n";
      }

      if (verbose_itr && printTrees){
        logger.log("Printing Current Status of Fit");

        logger.getVectorHead(z_, logBuff);
        // logger.log(logBuff);
        Rcout << "\n          z : " <<  logBuff << "\n";

        logger.getVectorHead(y, logBuff);
        Rcout << "          y : " <<  logBuff << "\n";

        logger.getVectorHead(allfit, logBuff);
        Rcout << "Fit - Tree  : " <<  logBuff << "\n";

        logger.getVectorHead(r_con, logBuff);
        Rcout << "     r_con  : " <<  logBuff << "\n\n";

        Rcout <<" MScale: " << mscale << "\n";

        Rcout <<" bscale0 : " << bscale0 << "\n";

        Rcout <<" bscale1 : " << bscale1 << "\n\n";

      }
      logger.log("Starting To Draw Mu");
      logger.startContext();

      if(use_mean_split_change) {
        ++mean_split_change_attempts;
        mean_split_change_accepts+=mean_split_change(t_con[iTreeCon],xi_con,di_con,weight,pi_con,gen);
      }
      drmu(t_con[iTreeCon],  // tree& x
           xi_con, // xinfo& xi
           di_con, // dinfo& di
           pi_con, // pinfo& pi,
           weight,
           gen); // RNG& gen

      logger.stopContext();

      logger.log("Attempting to Print Tree Post drmu \n");
      if(verbose_itr  && printTrees){
        t_con[iTreeCon].pr(xi_con);
        Rcout << "\n";
      }

      fit(t_con[iTreeCon],
          xi_con,
          di_con,
          ftemp);

      for(size_t k=0;k<n;k++) {
        allfit[k] += mscale*ftemp[k];
        allfit_con[k] += mscale*ftemp[k];
      }

      logger.log("Attempting to Print tree Post second call to fit \n");

      if(verbose_itr && printTrees){
        t_con[iTreeCon].pr(xi_con);
        Rcout << "\n";

      }
      logger.stopContext();
    }


    for(size_t iTreeMod=0;iTreeMod<ntree_mod;iTreeMod++) {
      logger.log("==================================");
      std::snprintf(logBuff, sizeof(logBuff), "Updating Moderate Tree: %d of %d", (int) iTreeMod + 1, ntree_mod);
      logger.log(logBuff);
      logger.log("==================================");
      logger.startContext();


      logger.log("Attempting to Print Tree Pre Update \n");
      if(verbose_itr && printTrees){
        t_mod[iTreeMod].pr(xi_mod);
        Rcout << "\n";
      }

      fit(t_mod[iTreeMod],
          xi_mod,
          di_mod,
          ftemp);

      logger.log("Attempting to Print Tree Post first call to fit");
      if(verbose_itr && printTrees){
        t_mod[iTreeMod].pr(xi_mod);
        Rcout << "\n";
      }

      for(size_t k=0;k<n;k++) {
        if(ftemp[k] != ftemp[k]) {
          Rcout << "moderator tree " << iTreeMod <<" obs "<< k<<" "<< endl;
          Rcout << t_mod[iTreeMod] << endl;
          stop("nan in ftemp");
        }
        double bscale = (k<ntrt) ? bscale1 : bscale0;
        allfit[k] = allfit[k]-bscale*ftemp[k];
        allfit_mod[k] = allfit_mod[k]-bscale*ftemp[k];
        r_mod[k] = (y[k]-allfit[k])/bscale;
      }
      logger.log("Starting Birth / Death Processing");
      logger.startContext();
      bd(t_mod[iTreeMod],
         xi_mod,
         di_mod,
         weight_het,
         pi_mod,
         gen,
         logger);
      logger.stopContext();

      logger.log("Attempting to Print Tree  Post bd \n");
      if(verbose_itr && printTrees){
        t_mod[iTreeMod].pr(xi_mod);
        Rcout << "\n";
      }

      if (verbose_itr && printTrees){
        logger.log("Printing Status of Fit");

        logger.getVectorHead(z_, logBuff);
        Rcout << "\n          z : " <<  logBuff << "\n";

        logger.getVectorHead(y, logBuff);
        Rcout << "          y : " <<  logBuff << "\n";

        logger.getVectorHead(allfit, logBuff);
        Rcout << "Fit - Tree  : " <<  logBuff << "\n";

        logger.getVectorHead(r_mod, logBuff);
        Rcout << "     r_mod  : " <<  logBuff << "\n\n";

        Rcout <<" MScale: " << mscale << "\n";

        Rcout <<" bscale0 : " << bscale0 << "\n";

        Rcout <<" bscale1 : " << bscale1 << "\n\n";

      }
      logger.log("Starting To Draw Mu");
      logger.startContext();
      if(use_mean_split_change) {
        ++mean_split_change_attempts;
        mean_split_change_accepts+=mean_split_change(t_mod[iTreeMod],xi_mod,di_mod,weight_het,pi_mod,gen);
      }
      drmu(t_mod[iTreeMod],
            xi_mod,
            di_mod,
            pi_mod,
            weight_het,
            gen);
      logger.stopContext();



      logger.log("Attempting to Print Tree Post drmuhet \n");
      if(verbose_itr && printTrees){
        t_mod[iTreeMod].pr(xi_mod);
        Rcout << "\n";
      }

      fit(t_mod[iTreeMod],
          xi_mod,
          di_mod,
          ftemp);

      for(size_t k=0;k<ntrt;k++) {
        allfit[k] += bscale1*ftemp[k];
        allfit_mod[k] += bscale1*ftemp[k];
      }
      for(size_t k=ntrt;k<n;k++) {
        allfit[k] += bscale0*ftemp[k];
        allfit_mod[k] += bscale0*ftemp[k];
      }

      logger.log("Attempting to Print Tree Post second call to fit");

      if(verbose_itr && printTrees){
        t_mod[iTreeMod].pr(xi_mod);
        Rcout << "\n";
      }
      logger.stopContext();

    } // end tree lop

    logger.setLevel(verbose_itr);

    logger.log("=====================================");
    logger.log("- MCMC iteration Cleanup");
    logger.log("=====================================");

    if(use_bscale) {
      double ww0 = 0.0, ww1 = 0.;
      double rw0 = 0.0, rw1 = 0.;
      for(size_t k=0; k<n; ++k) {
        double bscale = (k<ntrt) ? bscale1 : bscale0;
        double obs_var = use_hetero ? sigma2_fit[k] : sigma * sigma;
        double scale_factor = (w[k]*allfit_mod[k]*allfit_mod[k])/(obs_var*bscale*bscale);

        if(scale_factor!=scale_factor) {
          Rcout << " scale_factor " << scale_factor << endl;
          stop("");
        }

        double randeff_contrib = randeff ? allfit_random[k] : 0.0;

        double r = (y[k] - allfit_con[k] - randeff_contrib)*bscale/allfit_mod[k];

        if(r!=r) {
          Rcout << "bscale " << k << " r " << r << " mscale " <<mscale<< " b*z " << allfit_mod[k]*z_[k] << " bscale " << bscale0 << " " <<bscale1 << endl;
          stop("");
        }
        if(k<ntrt) {
          ww1 += scale_factor;
          rw1 += r*scale_factor;
        } else {
          ww0 += scale_factor;
          rw0 += r*scale_factor;
        }
      }
      logger.log("Drawing bscale 1");
      logger.startContext();
      double bscale1_old = bscale1;
      double bscale_fc_var = 1/(ww1 + bscale_prec);
      bscale1 = bscale_fc_var*rw1 + gen.normal(0., 1.)*sqrt(bscale_fc_var);
      if(verbose_itr){

        Rcout << "Original bscale1 : " << bscale1_old << "\n";
        Rcout << "bscale_prec : " << bscale_prec << ", ww1 : " << ww1 << ", rw1 : " << rw1 << "\n";
        Rcout << "New  bscale1 : " << bscale1 << "\n\n";
      }
      logger.stopContext();


      logger.log("Drawing bscale 0");
      logger.startContext();
      double bscale0_old = bscale0;
      bscale_fc_var = 1/(ww0 + bscale_prec);
      bscale0 = bscale_fc_var*rw0 + gen.normal(0., 1.)*sqrt(bscale_fc_var);
      if(verbose_itr){
        Rcout << "Original bscale0 : " << bscale0_old << "\n";
        Rcout << "bscale_prec : " << bscale_prec << ", ww0 : " << ww0 << ", rw0 : " << rw0 << "\n";
        Rcout << "New  bscale0 : " << bscale0 << "\n\n";
      }
      logger.stopContext();

      for(size_t k=0; k<ntrt; ++k) {
        allfit_mod[k] = allfit_mod[k]*bscale1/bscale1_old;
      }
      for(size_t k=ntrt; k<n; ++k) {
        allfit_mod[k] = allfit_mod[k]*bscale0/bscale0_old;
      }

      if(!b_half_normal) {
        double ssq = 0.0;
        tree::npv bnv;
        typedef tree::npv::size_type bvsz;
        double endnode_count = 0.0;

        for(size_t iTreeMod=0;iTreeMod<ntree_mod;iTreeMod++) {
          bnv.clear();
          t_mod[iTreeMod].getbots(bnv);
          bvsz nb = bnv.size();
          for(bvsz ii = 0; ii<nb; ++ii) {
            double mm = bnv[ii]->getm(); //node parameter
            // delta_mod is an absolute Gamma precision. Its conditional
            // uses the fixed base leaf variance, excluding delta_mod itself.
            ssq += mm*mm*((double)ntree_mod)/(mod_sd*mod_sd);
            endnode_count += 1.0;
          }
        }
        delta_mod = gen.gamma(0.5*(1. + endnode_count), 1.0)/(0.5*(1 + ssq));
      }
      if(verbose_itr){
        Rcout << "Original pi_mod.tau : " <<  pi_mod.tau << "\n";
      }

      pi_mod.tau   = mod_sd/(sqrt(delta_mod)*sqrt((double) ntree_mod));

      if(verbose_itr){
        Rcout << "New pi_mod.tau : " <<  pi_mod.tau << "\n\n";
      }

    } else {
      bscale0 = -0.5;
      bscale1 =  0.5;
    }
    pi_mod.sigma = sigma;


    if(use_mscale) {
      double ww = 0.;
      double rw = 0.;
      for(size_t k=0; k<n; ++k) {
        double obs_var = use_hetero ? sigma2_fit[k] : sigma * sigma;
        double scale_factor = (w[k]*allfit_con[k]*allfit_con[k])/(obs_var*mscale*mscale);
        if(scale_factor!=scale_factor) {
          Rcout << " scale_factor " << scale_factor << endl;
          stop("");
        }

        double randeff_contrib = randeff ? allfit_random[k] : 0.0;

        double r = (y[k] - allfit_mod[k]- randeff_contrib)*mscale/allfit_con[k];
        if(r!=r) {
          Rcout << "mscale " << k << " r " << r << " mscale " <<mscale<< " b*z " << allfit_mod[k]*z_[k] << " bscale " << bscale0 << " " <<bscale1 << endl;
          stop("");
        }
        ww += scale_factor;
        rw += r*scale_factor;
      }

      logger.log("Drawing mscale");


      double mscale_old = mscale;
      double mscale_fc_var = 1/(ww + mscale_prec);
      mscale = mscale_fc_var*rw + gen.normal(0., 1.)*sqrt(mscale_fc_var);
      if(verbose_itr){
        Rcout << "Original mscale : " << mscale_old << "\n";
        Rcout << "mscale_prec : " << mscale_prec << ", ww : " << ww << ", rw : " << rw << "\n";
        Rcout << "New  mscale : " << mscale << "\n\n";
      }


      //Rcout<< mscale_fc_var << " " << rw <<" " << mscale << endl;

      for(size_t k=0; k<n; ++k) {
        allfit_con[k] = allfit_con[k]*mscale/mscale_old;
      }

      if(collapsed_mu_scale) {
        // phi = sqrt(delta_old) * theta; rho = a / sqrt(delta_old).
        // Marginal rho is Cauchy(0,1), independent of fixed Gaussian phi.
        // rho | lambda ~ N(0,1/lambda), lambda | rho ~ Gamma(1,(1+rho^2)/2).
        mscale_prec = gen.gamma(1.0, 1.0)/(0.5*(1.0 + mscale*mscale));
        pi_con.tau = con_sd/sqrt((double)ntree_con);
      } else {
      // update delta_con

      double ssq = 0.0;
      tree::npv bnv;
      typedef tree::npv::size_type bvsz;
      double endnode_count = 0.0;

      for(size_t iTreeCon=0;iTreeCon<ntree_con;iTreeCon++) {
        bnv.clear();
        t_con[iTreeCon].getbots(bnv);
        bvsz nb = bnv.size();
        for(bvsz ii = 0; ii<nb; ++ii) {
          double mm = bnv[ii]->getm(); //node parameter
          // Var(leaf | delta_con) = con_sd^2 / (ntree_con * delta_con).
          // Conditioning on leaves therefore uses their fixed base variance.
          ssq += mm*mm*((double)ntree_con)/(con_sd*con_sd);
          endnode_count += 1.0;
        }
      }

      delta_con = gen.gamma(0.5*(1. + endnode_count), 1.0)/(0.5*(1 + ssq));
      if(verbose_itr){
        logger.log("Updating pi_con.tau");
        Rcout << "Original pi_con.tau : " <<  pi_con.tau << "\n";
      }

      pi_con.tau   = con_sd/(sqrt(delta_con)*sqrt((double) ntree_con));

      if(verbose_itr){
        Rcout << "New pi_con.tau : " <<  pi_con.tau << "\n\n";
      }

      }

    } else {
      mscale = 1.0;
    }
    pi_con.sigma = sigma/fabs(mscale); //should be sigma/abs(mscale) for backfitting

    //sync allfits after scale updates, if necessary. Could do smarter backfitting updates inline
    if(use_mscale || use_bscale) {
      logger.log("Sync allfits after scale updates");

      for(size_t k=0; k<n; ++k) {
        double randeff_contrib = randeff ? allfit_random[k] : 0.0;
        allfit[k] = allfit_con[k] + allfit_mod[k] + randeff_contrib;
      }
    }

    if(randeff) {
      Rcout << "==================================\n";
      Rcout << "- Random Effects \n";
      Rcout << "==================================\n";

      //update random effects
      for(size_t k=0; k<n; ++k) {
        r(k) = y[k] - allfit_con[k] - allfit_mod[k];
        allfit[k] -= allfit_random[k];
      }

      Wtr = random_des.t()*r;

      arma::mat adj = diagmat(random_var_ix*eta);
      //    Rcout << adj << endl << endl;
      arma::mat Phi = adj*WtW*adj/(sigma*sigma) + Sigma_inv_random;
      arma::vec m = adj*Wtr/(sigma*sigma);
      //Rcout << m << Phi << endl << Sigma_inv_random;
      gamma = rmvnorm_post(m, Phi);

      //Rcout << "updated gamma";

      // Update px parameters eta

      arma::mat adj2 = diagmat(gamma)*random_var_ix;
      arma::mat Phi2 = adj2.t()*WtW*adj2/(sigma*sigma) + arma::eye(eta.size(), eta.size());
      arma::vec m2 = adj2.t()*Wtr/(sigma*sigma);
      //Rcout << m << Phi << endl << Sigma_inv_random;
      eta = rmvnorm_post(m2, Phi2);

      //Rcout << "updated eta";

      // Update variance parameters

      arma::vec ssqs   = random_var_ix.t()*(gamma % gamma);
      //Rcout << "A";
      arma::rowvec counts = sum(random_var_ix, 0);
      //Rcout << "B";
      for(size_t ii=0; ii<random_var_ix.n_cols; ++ii) {
        random_var(ii) = 1.0/gen.gamma(0.5*(random_var_df + counts(ii)), 1.0)*2.0/(random_var_df + ssqs(ii));
      }
      //Rcout << "updated vars" << endl;
      Sigma_inv_random = diagmat(1/(random_var_ix*random_var));

      allfit_random = random_des*diagmat(random_var_ix*eta)*gamma;

      //Rcout << "recom allfit vars" << endl;

      for(size_t k=0; k<n; ++k) {
        allfit[k] = allfit_con[k] + allfit_mod[k] + allfit_random(k); //+= allfit_random[k];
      }
    }

    if(joint_mean_every>0 && (iIter+1)%joint_mean_every==0) {
      size_t leaves=joint_mean_refresh(t_con,t_mod,xi_con,xi_mod,di_con,di_mod,
        y.data(),w,sigma2_fit,ntrt,mscale,bscale0,bscale1,pi_con.tau,pi_mod.tau,
        allfit_con,allfit_mod,allfit,gen);
      ++joint_mean_updates;joint_mean_max_leaves=std::max(joint_mean_max_leaves,(double)leaves);
    }
    if(use_global_mean_update)
      global_mean_refresh(t_con,t_mod,y.data(),w,sigma2_fit,n,ntrt,
        mscale,bscale0,bscale1,pi_con.tau,pi_mod.tau,
        allfit_con,allfit_mod,allfit,gen);

    if(use_paired_mean_update) {
      for(size_t j=0;j<t_mod.size();++j) {
        bool ok=paired_mean_refresh(t_con[(j+iIter)%t_con.size()],t_mod[j],
          xi_con,xi_mod,di_con,di_mod,y.data(),w,sigma2_fit,ntrt,
          mscale,bscale0,bscale1,pi_con.tau,pi_mod.tau,
          allfit_con,allfit_mod,allfit,gen);
        if(ok) ++paired_mean_updates; else ++paired_mean_skips;
      }
    }

    if(use_hetero) {
      for(size_t iTreeVar=0; iTreeVar<(size_t) ntree_var; ++iTreeVar) {
        fit(t_var[iTreeVar], xi_var, di_var, ftemp_var);
        for(size_t k=0; k<n; ++k) {
          sigma0_2_fit[k] = sigma0_2_fit[k] / ftemp_var[k];
          sigma2_fit[k] = sigma0_2_fit[k] * ((use_ratio && k < (size_t)ntrt) ? sigma_ratio_fit[k] : 1.0);
          double restemp = y[k] - allfit[k];
          // w is a likelihood precision: Var(y_k | f_k) = sigma2_k / w_k.
          // The variance-tree likelihood must use the same weighted squared
          // residual as the mean/scale conditionals and homoscedastic update.
          r_var[k] = std::max(1e-12, w[k] * restemp * restemp / sigma2_fit[k]);
        }
        varbd(t_var[iTreeVar], xi_var, di_var, pi_var, var_nu, var_lambda, gen);
        if(use_variance_split_change){
          ++variance_split_change_attempts;
          variance_split_change_accepts+=variance_split_change(t_var[iTreeVar],xi_var,di_var,pi_var,var_nu,var_lambda,gen);
        }
        vardrmu(t_var[iTreeVar], xi_var, di_var, var_nu, var_lambda, gen);
        fit(t_var[iTreeVar], xi_var, di_var, ftemp_var);
        for(size_t k=0; k<n; ++k) {
          sigma0_2_fit[k] *= ftemp_var[k];
          sigma2_fit[k] = sigma0_2_fit[k] * ((use_ratio && k < (size_t)ntrt) ? sigma_ratio_fit[k] : 1.0);
          if(!std::isfinite(sigma2_fit[k]) || sigma2_fit[k] <= 0.0) {
            stop("nonpositive or non-finite sigma2_fit in heteroscedastic variance update");
          }
        }
      }
      if(use_ratio) {
        for(size_t iTreeVarRatio=0; iTreeVarRatio<(size_t) ntree_var; ++iTreeVarRatio) {
          fit(t_var_ratio[iTreeVarRatio], xi_var, di_var_ratio_all, ftemp_var_ratio);
          for(size_t k=0; k<n; ++k) {
            sigma_ratio_fit[k] = sigma_ratio_fit[k] / ftemp_var_ratio[k];
          }
          for(size_t k=0; k<(size_t)ntrt; ++k) {
            sigma2_fit[k] = sigma0_2_fit[k] * sigma_ratio_fit[k];
            double restemp = y[k] - allfit[k];
            r_var_ratio[k] = std::max(1e-12, w[k] * restemp * restemp / sigma2_fit[k]);
          }
          varbd(t_var_ratio[iTreeVarRatio], xi_var, di_var_ratio, pi_var,
                var_nu, var_lambda, gen);
          if(use_variance_split_change){
            ++variance_split_change_attempts;
            variance_split_change_accepts+=variance_split_change(t_var_ratio[iTreeVarRatio],xi_var,di_var_ratio,pi_var,var_nu,var_lambda,gen);
          }
          vardrmu(t_var_ratio[iTreeVarRatio], xi_var, di_var_ratio,
                  var_nu, var_lambda, gen);
          fit(t_var_ratio[iTreeVarRatio], xi_var, di_var_ratio_all, ftemp_var_ratio);
          for(size_t k=0; k<n; ++k) {
            sigma_ratio_fit[k] *= ftemp_var_ratio[k];
            if(!std::isfinite(sigma_ratio_fit[k]) || sigma_ratio_fit[k] <= 0.0) {
              stop("nonpositive or non-finite sigma_ratio_fit in heteroscedastic ratio update");
            }
            sigma2_fit[k] = sigma0_2_fit[k] * ((k < (size_t)ntrt) ? sigma_ratio_fit[k] : 1.0);
            if(!std::isfinite(sigma2_fit[k]) || sigma2_fit[k] <= 0.0) {
              stop("nonpositive or non-finite sigma2_fit in heteroscedastic ratio update");
            }
          }
        }
      }
      if(joint_variance_every>0 && (iIter+1)%joint_variance_every==0) {
        ++joint_variance_attempts;
        joint_variance_accepts+=joint_variance_refresh(t_var,t_var_ratio,y.data(),allfit,w,
          sigma0_2_fit,sigma_ratio_fit,sigma2_fit,n,ntrt,var_nu,var_lambda,gen);
      }
      if(paired_variance_every>0 && (iIter+1)%paired_variance_every==0) {
        for(size_t j=0;j<t_var.size();++j){
          size_t chosen=std::min((size_t)(gen.uniform()*ntrt),(size_t)ntrt-1);
          ++paired_variance_attempts;
          paired_variance_accepts+=paired_variance_refresh(t_var[j],t_var_ratio[(j+iIter)%t_var_ratio.size()],
            xi_var,di_var.x,di_var.p,chosen,y.data(),allfit,w,sigma0_2_fit,sigma_ratio_fit,
            sigma2_fit,n,ntrt,var_nu,var_lambda,gen);
        }
      }
      double mean_sigma2 = 0.0;
      for(size_t k=0; k<n; ++k) mean_sigma2 += sigma2_fit[k];
      sigma = sqrt(mean_sigma2 / (double) n);
    } else {
      // ---------------------------------------------------------
      logger.log("Draw Sigma");
      // ---------------------------------------------------------
      double rss = 0.0;
      double restemp = 0.0;
      for(size_t k=0;k<n;k++) {
        restemp = y[k]-allfit[k];
        rss += w[k]*restemp*restemp;
      }
      sigma = sqrt((nu*lambda + rss)/gen.chi_square(nu+n));
    }
    pi_con.sigma = sigma/fabs(mscale);
    pi_mod.sigma = sigma; // Is this another copy paste Error?

    if( ((iIter>=burn) & (iIter % thin==0)) )  {
      if(not treef_con_name.empty()){
        for(size_t j=0;j<ntree_con;j++) treef_con << std::setprecision(save_tree_precision) << t_con[j] << "\n"; // save trees
        for(size_t j=0;j<ntree_mod;j++) treef_mod << std::setprecision(save_tree_precision) << t_mod[j] << "\n"; // save trees
        if(use_hetero && !treef_var_name.empty()) {
          for(size_t j=0;j<(size_t) ntree_var;j++) treef_var << std::setprecision(save_tree_precision) << t_var[j] << "\n";
        }
        if(use_ratio && !treef_var_ratio_name.empty()) {
          for(size_t j=0;j<(size_t) ntree_var;j++) treef_var_ratio << std::setprecision(save_tree_precision) << t_var_ratio[j] << "\n";
        }
      }

      msd_post(save_ctr) = mscale;
      bsd_post(save_ctr) = bscale1-bscale0;
      b0_post(save_ctr)  = bscale0;
      b1_post(save_ctr)  = bscale1;


      gamma_post.row(save_ctr) = (diagmat(random_var_ix*eta)*gamma).t();
      random_var_post.row(save_ctr) = (sqrt( eta % eta % random_var)).t();

      sigma_post(save_ctr) = sigma;
      for(size_t k=0;k<n;k++) {
        m_post(save_ctr, k) = allfit_con[k];
        yhat_post(save_ctr, k) = allfit[k];
        sigma2_post(save_ctr, k) = use_hetero ? sigma2_fit[k] : sigma * sigma;
        sigma0_2_post(save_ctr, k) = use_hetero ? sigma0_2_fit[k] : sigma * sigma;
        sigma1_2_post(save_ctr, k) = use_hetero ? sigma0_2_fit[k] * (use_ratio ? sigma_ratio_fit[k] : 1.0) : sigma * sigma;
        log_var_ratio_post(save_ctr, k) = use_ratio ? std::log(sigma_ratio_fit[k]) : 0.0;
      }
      for(size_t k=0;k<n;k++) {
        double bscale = (k<ntrt) ? bscale1 : bscale0;
        b_post(save_ctr, k) = (bscale1-bscale0)*allfit_mod[k]/bscale;
      }
      //}
      save_ctr += 1;
    }
    logger.log("==============================================");
    std::snprintf(logBuff, sizeof(logBuff), "MCMC iteration: %d of %d End", (int) iIter + 1, nd*thin+burn);
    logger.log(logBuff);
    std::snprintf(logBuff, sizeof(logBuff), "sigma %f, mscale %f, bscale0 %f, bscale1 %f",sigma, mscale, bscale0, bscale1);
    logger.log(logBuff);
    logger.log("==============================================");
    if (verbose_itr){
      logger.getVectorHead(y, logBuff);
      Rcout << "           y: " <<  logBuff << "\n";

      logger.getVectorHead(allfit, logBuff);
      Rcout << "Current Fit : " <<  logBuff << "\n";

      logger.getVectorHead(allfit_con, logBuff);
      Rcout << "allfit_con  : " <<  logBuff << "\n";

      logger.getVectorHead(allfit_mod, logBuff);
      Rcout << "allfit_mod  : " <<  logBuff << "\n";
    }

  } // end MCMC Loop

  int time2 = time(&tp);
  Rcout << "\n============================================================\n MCMC Complete \n============================================================\n";

  Rcout << "time for loop: " << time2 - time1 << endl;

  t_mod.clear(); t_con.clear(); t_var.clear(); t_var_ratio.clear();
  delete[] allfit;
  delete[] allfit_mod;
  delete[] allfit_con;
  delete[] sigma2_fit;
  delete[] sigma0_2_fit;
  delete[] sigma_ratio_fit;
  delete[] r_var;
  delete[] r_var_ratio;
  delete[] ftemp_var;
  delete[] ftemp_var_ratio;
  delete[] r_mod;
  delete[] r_con;
  delete[] ftemp;

  if(not treef_con_name.empty()){
    treef_con.close();
    treef_mod.close();
    if(use_hetero && !treef_var_name.empty()) treef_var.close();
    if(use_ratio && !treef_var_ratio_name.empty()) treef_var_ratio.close();
  }

  return(List::create(_["yhat_post"] = yhat_post, _["m_post"] = m_post, _["b_post"] = b_post,
                      _["sigma2_post"] = sigma2_post,
                      _["sigma0_2_post"] = sigma0_2_post,
                      _["sigma1_2_post"] = sigma1_2_post,
                      _["log_var_ratio_post"] = log_var_ratio_post,
                      _["sigma"] = sigma_post, _["msd"] = msd_post, _["bsd"] = bsd_post, _["b0"] = b0_post, _["b1"] = b1_post,
                      _["gamma"] = gamma_post, _["random_var_post"] = random_var_post,
                      _["collapsed_mu_scale"] = collapsed_mu_scale,
                      _["paired_variance_attempts"] = paired_variance_attempts,
                      _["paired_variance_accepts"] = paired_variance_accepts,
                      _["joint_variance_attempts"] = joint_variance_attempts,
                      _["joint_variance_accepts"] = joint_variance_accepts,
                      _["variance_split_change_attempts"] = variance_split_change_attempts,
                      _["variance_split_change_accepts"] = variance_split_change_accepts,
                      _["mean_split_change_attempts"] = mean_split_change_attempts,
                      _["mean_split_change_accepts"] = mean_split_change_accepts,
                      _["joint_mean_updates"] = joint_mean_updates,
                      _["joint_mean_max_leaves"] = joint_mean_max_leaves,
                      _["global_mean_update"] = use_global_mean_update,
                      _["paired_mean_update"] = use_paired_mean_update,
                      _["paired_mean_updates"] = paired_mean_updates,
                      _["paired_mean_skips"] = paired_mean_skips
  ));
}
