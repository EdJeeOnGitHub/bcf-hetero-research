// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include "bd.h"
#include "funs.h"
// Test fixture links to the installed package's actual MH kernel.
// [[Rcpp::export]]
Rcpp::List split_change_reference(Rcpp::NumericMatrix X,Rcpp::NumericVector y,
  Rcpp::NumericVector phi,Rcpp::List cuts,int draws,int ancestor=-1,double tau=.5) {
  int n=X.nrow(),p=X.ncol();
  std::vector<double> xx(n*p);xinfo xi(p);
  for(int j=0;j<p;++j){Rcpp::NumericVector c=cuts[j];xi[j]=Rcpp::as<std::vector<double>>(c);
    for(int i=0;i<n;++i)xx[i*p+j]=X(i,j);}
  dinfo di;di.n=n;di.p=p;di.x=xx.data();di.y=y.begin();
  pinfo pi;pi.alpha=.95;pi.beta=2;pi.tau=tau;pi.sigma=1;
  tree t;tree::tree_p node=&t;
  if(ancestor>=0){t.birth(1,0,ancestor,0,0);node=t.getl();}
  t.birth(node->nid(),0,0,0,0);
  std::vector<int> offsets(p+1,0);
  for(int j=0;j<p;++j)offsets[j+1]=offsets[j]+xi[j].size();
  Rcpp::IntegerVector counts(offsets[p]);int accepted=0;RNG gen;
  for(int i=0;i<draws+10000;++i){
    accepted+=mean_split_change(t,xi,di,phi.begin(),pi,gen);
    if(i>=10000)++counts[offsets[node->getv()]+node->getc()];
  }
  return Rcpp::List::create(Rcpp::Named("counts")=counts,Rcpp::Named("accepted")=accepted);
}
