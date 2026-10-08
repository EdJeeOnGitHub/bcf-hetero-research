// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include "varfuns.h"
#include "funs.h"
// Test fixture links to the installed package's actual MH kernel.
// [[Rcpp::export]]
Rcpp::List variance_change_reference(Rcpp::NumericMatrix X,Rcpp::NumericVector y,
  Rcpp::NumericVector phi,Rcpp::List cuts,int draws,int ancestor=-1,double nu=8,double lambda=.6) {
  int n=X.nrow(),p=X.ncol();
  std::vector<double> xx(n*p);xinfo xi(p);
  for(int j=0;j<p;++j){Rcpp::NumericVector c=cuts[j];xi[j]=Rcpp::as<std::vector<double>>(c);
    for(int i=0;i<n;++i)xx[i*p+j]=X(i,j);}
  dinfo di;di.n=n;di.p=p;di.x=xx.data();std::vector<double> stat(n);for(int i=0;i<n;++i)stat[i]=phi[i]*y[i]*y[i];di.y=stat.data();
  pinfo pi;pi.alpha=.95;pi.beta=2;pi.tau=1;pi.sigma=1;
  tree t;tree::tree_p node=&t;
  if(ancestor>=0){t.birth(1,0,ancestor,1,1);node=t.getl();}
  t.birth(node->nid(),0,0,1,1);
  std::vector<int> offsets(p+1,0);
  for(int j=0;j<p;++j)offsets[j+1]=offsets[j]+xi[j].size();
  Rcpp::IntegerVector counts(offsets[p]);int accepted=0;RNG gen;
  for(int i=0;i<draws+10000;++i){
    accepted+=variance_split_change(t,xi,di,pi,nu,lambda,gen);
    if(i>=10000)++counts[offsets[node->getv()]+node->getc()];
  }
  return Rcpp::List::create(Rcpp::Named("counts")=counts,Rcpp::Named("accepted")=accepted);
}
