/*
 * Scalar product-of-trees variance updates for heteroscedastic BCF.
 *
 * This is the scalar residual-variance subset of the bayesm.HART varbart
 * machinery: each variance tree has positive inverse-chi-square leaves and
 * the full variance function is the product of tree predictions.
 */
#ifndef GUARD_varfuns_h
#define GUARD_varfuns_h

#include "tree.h"
#include "info.h"
#include "rng.h"

double varlh(double S, double n, double nu, double lambda);

void vargetsuffBirth(tree& x, tree::tree_cp nx, size_t v, size_t c,
                     xinfo& xi, dinfo& di,
                     size_t& nl, double& Sl,
                     size_t& nr, double& Sr);

void vargetsuffDeath(tree& x, tree::tree_cp l, tree::tree_cp r,
                     xinfo& xi, dinfo& di,
                     size_t& nl, double& Sl,
                     size_t& nr, double& Sr);

double vardrawnodemu(double S, double n, double nu, double lambda, RNG& gen);

void varallsuff(tree& x, xinfo& xi, dinfo& di,
                tree::npv& bnv,
                std::vector<size_t>& nv_leaf,
                std::vector<double>& Sv);

void vardrmu(tree& t, xinfo& xi, dinfo& di,
             double nu, double lambda, RNG& gen);

bool varbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi,
           double nu, double lambda, RNG& gen);

#endif
