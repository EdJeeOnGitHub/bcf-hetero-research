#include "varfuns.h"
#include "funs.h"

#include <cmath>
#include <map>

double varlh(double S, double n, double nu, double lambda)
{
  double half = 0.5 * (n + nu);
  double half_nu = 0.5 * nu;
  return -0.5 * n * log(M_PI)
       + half_nu * log(nu * lambda)
       - lgamma(half_nu)
       + lgamma(half)
       - half * log(S + nu * lambda);
}

void vargetsuffBirth(tree& x, tree::tree_cp nx, size_t v, size_t c,
                     xinfo& xi, dinfo& di,
                     size_t& nl, double& Sl,
                     size_t& nr, double& Sr)
{
  nl = 0; nr = 0; Sl = 0.0; Sr = 0.0;
  for(size_t i=0; i<di.n; ++i) {
    double* xx = di.x + i * di.p;
    if(nx == x.bn(xx, xi)) {
      if(xx[v] < xi[v][c]) {
        nl += 1;
        Sl += di.y[i];
      } else {
        nr += 1;
        Sr += di.y[i];
      }
    }
  }
}

void vargetsuffDeath(tree& x, tree::tree_cp l, tree::tree_cp r,
                     xinfo& xi, dinfo& di,
                     size_t& nl, double& Sl,
                     size_t& nr, double& Sr)
{
  nl = 0; nr = 0; Sl = 0.0; Sr = 0.0;
  for(size_t i=0; i<di.n; ++i) {
    double* xx = di.x + i * di.p;
    tree::tree_cp bn = x.bn(xx, xi);
    if(bn == l) {
      nl += 1;
      Sl += di.y[i];
    } else if(bn == r) {
      nr += 1;
      Sr += di.y[i];
    }
  }
}

double vardrawnodemu(double S, double n, double nu, double lambda, RNG& gen)
{
  return (nu * lambda + S) / gen.chi_square(nu + n);
}

void varallsuff(tree& x, xinfo& xi, dinfo& di,
                tree::npv& bnv,
                std::vector<size_t>& nv_leaf,
                std::vector<double>& Sv)
{
  bnv.clear();
  x.getbots(bnv);

  typedef tree::npv::size_type bvsz;
  bvsz nb = bnv.size();
  nv_leaf.assign(nb, 0);
  Sv.assign(nb, 0.0);

  std::map<tree::tree_cp, size_t> bnmap;
  for(bvsz i=0; i<nb; ++i) bnmap[bnv[i]] = i;

  for(size_t i=0; i<di.n; ++i) {
    double* xx = di.x + i * di.p;
    size_t ni = bnmap[x.bn(xx, xi)];
    nv_leaf[ni] += 1;
    Sv[ni] += di.y[i];
  }
}

void vardrmu(tree& t, xinfo& xi, dinfo& di,
             double nu, double lambda, RNG& gen)
{
  tree::npv bnv;
  std::vector<size_t> nv_leaf;
  std::vector<double> Sv;
  varallsuff(t, xi, di, bnv, nv_leaf, Sv);
  for(tree::npv::size_type i=0; i<bnv.size(); ++i) {
    bnv[i]->setm(vardrawnodemu(Sv[i], static_cast<double>(nv_leaf[i]), nu, lambda, gen));
  }
}

bool varbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi,
           double nu, double lambda, RNG& gen)
{
  tree::npv goodbots;
  double PBx = getpb(x, xi, pi, goodbots);

  if(gen.uniform() < PBx) {
    size_t ni = floor(gen.uniform() * goodbots.size());
    tree::tree_p nx = goodbots[ni];

    std::vector<size_t> goodvars;
    getgoodvars(nx, xi, goodvars);
    size_t vi = floor(gen.uniform() * goodvars.size());
    size_t v = goodvars[vi];

    int L = 0, U = xi[v].size() - 1;
    nx->rg(v, &L, &U);
    size_t c = L + floor(gen.uniform() * (U - L + 1));

    double Pbotx = 1.0 / goodbots.size();
    size_t dnx = nx->depth();
    double PGnx = pi.alpha / pow(1.0 + dnx, pi.beta);

    double PGly, PGry;
    if(goodvars.size() > 1) {
      PGly = pi.alpha / pow(1.0 + dnx + 1.0, pi.beta);
      PGry = PGly;
    } else {
      PGly = (((int)c - 1) < L) ? 0.0 : pi.alpha / pow(1.0 + dnx + 1.0, pi.beta);
      PGry = (U < (int)(c + 1)) ? 0.0 : pi.alpha / pow(1.0 + dnx + 1.0, pi.beta);
    }

    double PDy;
    if(goodbots.size() > 1) {
      PDy = 1.0 - pi.pb;
    } else {
      PDy = ((PGry == 0) && (PGly == 0)) ? 1.0 : 1.0 - pi.pb;
    }

    size_t nnogs = x.nnogs();
    tree::tree_cp nxp = nx->getp();
    double Pnogy;
    if(nxp == 0) {
      Pnogy = 1.0;
    } else {
      Pnogy = nxp->isnog() ? 1.0 / nnogs : 1.0 / (nnogs + 1.0);
    }

    size_t nl, nr;
    double Sl, Sr;
    vargetsuffBirth(x, nx, v, c, xi, di, nl, Sl, nr, Sr);

    double alpha = 0.0;
    if((nl >= 5) && (nr >= 5)) {
      double n_l = static_cast<double>(nl);
      double n_r = static_cast<double>(nr);
      double lhl = varlh(Sl, n_l, nu, lambda);
      double lhr = varlh(Sr, n_r, nu, lambda);
      double lht = varlh(Sl + Sr, n_l + n_r, nu, lambda);
      double alpha1 = (PGnx * (1.0 - PGly) * (1.0 - PGry) * PDy * Pnogy) /
                      ((1.0 - PGnx) * PBx * Pbotx);
      alpha = std::min(1.0, alpha1 * exp(lhl + lhr - lht));
    }

    if(gen.uniform() < alpha) {
      double mul = vardrawnodemu(Sl, static_cast<double>(nl), nu, lambda, gen);
      double mur = vardrawnodemu(Sr, static_cast<double>(nr), nu, lambda, gen);
      x.birth(nx->nid(), v, c, mul, mur);
      return true;
    }
    return false;
  }

  tree::npv nognds;
  x.getnogs(nognds);
  size_t ni = floor(gen.uniform() * nognds.size());
  tree::tree_p nx = nognds[ni];

  size_t dny = nx->depth();
  double PGny = pi.alpha / pow(1.0 + dny, pi.beta);
  double PGlx = pgrow(nx->getl(), xi, pi);
  double PGrx = pgrow(nx->getr(), xi, pi);
  double PBy = (!(nx->getp())) ? 1.0 : pi.pb;

  int ngood = goodbots.size();
  if(cansplit(nx->getl(), xi)) --ngood;
  if(cansplit(nx->getr(), xi)) --ngood;
  ++ngood;
  double Pboty = 1.0 / ngood;
  double PDx = 1.0 - PBx;
  double Pnogx = 1.0 / nognds.size();

  size_t nl, nr;
  double Sl, Sr;
  vargetsuffDeath(x, nx->getl(), nx->getr(), xi, di, nl, Sl, nr, Sr);

  double n_l = static_cast<double>(nl);
  double n_r = static_cast<double>(nr);
  double lhl = varlh(Sl, n_l, nu, lambda);
  double lhr = varlh(Sr, n_r, nu, lambda);
  double lht = varlh(Sl + Sr, n_l + n_r, nu, lambda);

  double alpha1 = ((1.0 - PGny) * PBy * Pboty) /
                  (PGny * (1.0 - PGlx) * (1.0 - PGrx) * PDx * Pnogx);
  double alpha = std::min(1.0, alpha1 * exp(lht - lhl - lhr));

  if(gen.uniform() < alpha) {
    double mu = vardrawnodemu(Sl + Sr, n_l + n_r, nu, lambda, gen);
    x.death(nx->nid(), mu);
    return true;
  }
  return false;
}
