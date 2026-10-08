// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "joint_variance.h"
// [[Rcpp::export]]
Rcpp::List joint_variance_reference(Rcpp::NumericVector yy,Rcpp::NumericVector ww,
 int nt,int draws,double nu,double lambda) {
 size_t n=yy.size();std::vector<tree> g(2),r(2);
 g[0].setm(.7);g[1].setm(1.3);r[0].setm(1.2);r[1].setm(.6);
 std::vector<double> mean(n,0),v0(n,.91),vr(n,.72),obs(n);
 for(size_t k=0;k<n;++k)obs[k]=v0[k]*(k<(size_t)nt?vr[k]:1);
 Rcpp::NumericMatrix out(draws,2);RNG gen;int accepted=0;
 for(int i=0;i<draws+10000;++i){
  accepted+=joint_variance_refresh(g,r,yy.begin(),mean.data(),ww.begin(),v0.data(),vr.data(),obs.data(),n,nt,nu,lambda,gen);
  double pg=g[0].getm()*g[1].getm(),pr=r[0].getm()*r[1].getm();
  for(size_t k=0;k<n;++k)if(std::abs(obs[k]-pg*(k<(size_t)nt?pr:1))>1e-8 || std::abs(v0[k]-pg)>1e-8 ||std::abs(vr[k]-pr)>1e-8)Rcpp::stop("Cache/leaf inconsistency");
  if(i>=10000){out(i-10000,0)=log(pg/.91);out(i-10000,1)=log(pr/.72);}
 }
 return Rcpp::List::create(Rcpp::Named("draws")=out,Rcpp::Named("accepted")=accepted);
}
