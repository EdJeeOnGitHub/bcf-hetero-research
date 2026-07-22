#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "funs.h"
#include "rng.h"
#include "tree.h"

using namespace Rcpp;

namespace {

xinfo sbc_xinfo(const List& cutpoints)
{
  xinfo xi(cutpoints.size());
  for(R_xlen_t j = 0; j < cutpoints.size(); ++j) {
    NumericVector cp = cutpoints[j];
    xi[j] = as<std::vector<double> >(cp);
  }
  return xi;
}

std::vector<double> sbc_row_major(const NumericMatrix& x)
{
  std::vector<double> out(x.nrow() * x.ncol());
  for(int i = 0; i < x.nrow(); ++i) {
    for(int j = 0; j < x.ncol(); ++j) out[i * x.ncol() + j] = x(i, j);
  }
  return out;
}

bool sbc_valid_leaf_counts(tree& t, xinfo& xi,
                           const std::vector<double>& x,
                           size_t p, const std::vector<size_t>& rows,
                           size_t min_leaf)
{
  std::map<size_t, size_t> counts;
  tree::npv leaves;
  t.getbots(leaves);
  for(tree::npv::const_iterator it = leaves.begin(); it != leaves.end(); ++it) {
    counts[(*it)->nid()] = 0;
  }
  for(std::vector<size_t>::const_iterator it = rows.begin(); it != rows.end(); ++it) {
    double* xx = const_cast<double*>(&x[(*it) * p]);
    counts[t.bn(xx, xi)->nid()] += 1;
  }
  for(std::map<size_t, size_t>::const_iterator it = counts.begin(); it != counts.end(); ++it) {
    if(it->second < min_leaf) return false;
  }
  return true;
}

void sbc_draw_topology(tree& t, xinfo& xi, double alpha, double beta, RNG& gen)
{
  std::vector<size_t> pending(1, 1);
  while(!pending.empty()) {
    size_t nid = pending.back();
    pending.pop_back();
    tree::tree_p node = t.getptr(nid);
    if(node == 0 || !cansplit(node, xi)) continue;
    double p_grow = alpha / std::pow(1.0 + node->depth(), beta);
    if(gen.uniform() >= p_grow) continue;

    std::vector<size_t> goodvars;
    getgoodvars(node, xi, goodvars);
    size_t v = goodvars[static_cast<size_t>(std::floor(gen.uniform() * goodvars.size()))];
    int L = 0;
    int U = static_cast<int>(xi[v].size()) - 1;
    node->rg(v, &L, &U);
    size_t c = static_cast<size_t>(L + std::floor(gen.uniform() * (U - L + 1)));
    if(!t.birth(nid, v, c, 0.0, 0.0)) stop("failed to grow an SBC prior tree");
    pending.push_back(2 * nid + 1);
    pending.push_back(2 * nid);
  }
}

tree sbc_draw_tree(xinfo& xi, const std::vector<double>& x, size_t p,
                   const std::vector<size_t>& support_rows, size_t min_leaf,
                   double alpha, double beta, double leaf_df,
                   double leaf_scale, bool variance_leaf, int max_attempts,
                   RNG& gen, std::vector<double>& leaf_values)
{
  for(int attempt = 0; attempt < max_attempts; ++attempt) {
    tree candidate;
    sbc_draw_topology(candidate, xi, alpha, beta, gen);
    if(!sbc_valid_leaf_counts(candidate, xi, x, p, support_rows, min_leaf)) continue;

    tree::npv leaves;
    candidate.getbots(leaves);
    for(tree::npv::iterator it = leaves.begin(); it != leaves.end(); ++it) {
      double value = variance_leaf
        ? leaf_df * leaf_scale / gen.chi_square(leaf_df)
        : gen.normal(0.0, leaf_scale);
      (*it)->setm(value);
      leaf_values.push_back(value);
    }
    return candidate;
  }
  stop("could not draw a prior tree satisfying the minimum leaf-size support");
}

std::vector<double> sbc_predict_tree(tree& t, xinfo& xi,
                                     const std::vector<double>& x,
                                     size_t n, size_t p)
{
  std::vector<double> out(n);
  for(size_t i = 0; i < n; ++i) {
    double* xx = const_cast<double*>(&x[i * p]);
    out[i] = t.bn(xx, xi)->getm();
  }
  return out;
}

List sbc_draw_forest(const NumericMatrix& x_matrix, const List& cutpoints,
                     int num_trees, double alpha, double beta,
                     double leaf_df, double leaf_scale, bool variance_leaf,
                     const std::vector<size_t>& support_rows, int min_leaf,
                     int max_attempts, RNG& gen)
{
  xinfo xi = sbc_xinfo(cutpoints);
  std::vector<double> x = sbc_row_major(x_matrix);
  size_t n = x_matrix.nrow();
  size_t p = x_matrix.ncol();
  NumericVector fitted(n, variance_leaf ? 1.0 : 0.0);
  IntegerVector tree_sizes(num_trees);
  IntegerVector tree_depths(num_trees);
  std::vector<double> leaf_values;

  for(int j = 0; j < num_trees; ++j) {
    tree t = sbc_draw_tree(xi, x, p, support_rows, min_leaf, alpha, beta,
                           leaf_df, leaf_scale, variance_leaf, max_attempts,
                           gen, leaf_values);
    std::vector<double> pred = sbc_predict_tree(t, xi, x, n, p);
    for(size_t i = 0; i < n; ++i) {
      if(variance_leaf) fitted[i] *= pred[i]; else fitted[i] += pred[i];
    }
    tree_sizes[j] = static_cast<int>(t.treesize());
    tree::cnpv nodes;
    t.getnodes(nodes);
    size_t max_depth = 0;
    for(tree::cnpv::const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
      max_depth = std::max(max_depth, (*it)->depth());
    }
    tree_depths[j] = static_cast<int>(max_depth);
  }

  return List::create(
    _["fitted"] = fitted,
    _["tree_sizes"] = tree_sizes,
    _["tree_depths"] = tree_depths,
    _["leaf_values"] = wrap(leaf_values)
  );
}

} // namespace

// Draw a joint BCF prior realization for simulation-based calibration.
// Internal test helper; hyperparameters must already be on the fitting scale.
// [[Rcpp::export]]
List bcfSbcPriorDraw(NumericMatrix x_control, NumericMatrix x_moderate,
                     NumericMatrix x_variance, NumericVector z,
                     List cutpoints_control, List cutpoints_moderate,
                     List cutpoints_variance, std::string variance_model,
                     int ntree_control, int ntree_moderate, int ntree_variance,
                     double sd_control, double sd_moderate,
                     double base_control, double power_control,
                     double base_moderate, double power_moderate,
                     double variance_base, double variance_power,
                     double sigma_nu, double sigma_lambda,
                     double variance_nu_tree, double variance_lambda_tree,
                     int min_leaf = 5, int max_attempts = 100000)
{
  int n = x_control.nrow();
  if(n <= 0 || x_moderate.nrow() != n || x_variance.nrow() != n || z.size() != n) {
    stop("SBC prior inputs must have the same positive number of observations");
  }
  if(variance_model != "homoscedastic" && variance_model != "shared" &&
     variance_model != "ratio") stop("unknown SBC variance model");
  if(ntree_control <= 0 || ntree_moderate <= 0 || min_leaf <= 0) {
    stop("SBC tree counts and min_leaf must be positive");
  }
  if(variance_model != "homoscedastic" && ntree_variance <= 0) {
    stop("heteroscedastic SBC requires positive variance-tree count");
  }

  std::vector<size_t> all_rows(n);
  std::vector<size_t> treated_rows;
  for(int i = 0; i < n; ++i) {
    all_rows[i] = static_cast<size_t>(i);
    if(z[i] == 1.0) treated_rows.push_back(static_cast<size_t>(i));
  }
  if(variance_model == "ratio" && treated_rows.size() < static_cast<size_t>(2 * min_leaf)) {
    stop("ratio SBC needs at least 2 * min_leaf treated observations");
  }

  RNG gen;
  List control = sbc_draw_forest(
    x_control, cutpoints_control, ntree_control, base_control, power_control,
    0.0, sd_control / std::sqrt(static_cast<double>(ntree_control)), false,
    all_rows, min_leaf, max_attempts, gen
  );
  List moderate = sbc_draw_forest(
    x_moderate, cutpoints_moderate, ntree_moderate, base_moderate, power_moderate,
    0.0, sd_moderate / std::sqrt(static_cast<double>(ntree_moderate)), false,
    all_rows, min_leaf, max_attempts, gen
  );

  NumericVector mu = control["fitted"];
  NumericVector tau = moderate["fitted"];
  NumericVector sigma0_2(n);
  NumericVector sigma1_2(n);
  NumericVector log_var_ratio(n, 0.0);
  List baseline;
  List ratio;

  if(variance_model == "homoscedastic") {
    double sigma2 = sigma_nu * sigma_lambda / gen.chi_square(sigma_nu);
    std::fill(sigma0_2.begin(), sigma0_2.end(), sigma2);
    std::fill(sigma1_2.begin(), sigma1_2.end(), sigma2);
  } else {
    baseline = sbc_draw_forest(
      x_variance, cutpoints_variance, ntree_variance, variance_base,
      variance_power, variance_nu_tree, variance_lambda_tree, true,
      all_rows, min_leaf, max_attempts, gen
    );
    sigma0_2 = as<NumericVector>(baseline["fitted"]);
    if(variance_model == "ratio") {
      ratio = sbc_draw_forest(
        x_variance, cutpoints_variance, ntree_variance, variance_base,
        variance_power, variance_nu_tree, variance_lambda_tree, true,
        treated_rows, min_leaf, max_attempts, gen
      );
      NumericVector ratio_fit = ratio["fitted"];
      for(int i = 0; i < n; ++i) {
        sigma1_2[i] = sigma0_2[i] * ratio_fit[i];
        log_var_ratio[i] = std::log(ratio_fit[i]);
      }
    } else {
      sigma1_2 = clone(sigma0_2);
    }
  }

  NumericVector sigma2(n);
  NumericVector mean(n);
  for(int i = 0; i < n; ++i) {
    sigma2[i] = z[i] == 1.0 ? sigma1_2[i] : sigma0_2[i];
    mean[i] = mu[i] + z[i] * tau[i];
  }

  return List::create(
    _["mu"] = mu,
    _["tau"] = tau,
    _["mean"] = mean,
    _["sigma2"] = sigma2,
    _["sigma0_2"] = sigma0_2,
    _["sigma1_2"] = sigma1_2,
    _["log_var_ratio"] = log_var_ratio,
    _["control"] = control,
    _["moderate"] = moderate,
    _["variance"] = baseline,
    _["variance_ratio"] = ratio
  );
}
