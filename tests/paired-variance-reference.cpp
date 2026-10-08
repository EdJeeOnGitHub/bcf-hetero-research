// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "paired_variance.h"
// [[Rcpp::export]]
Rcpp::List paired_variance_reference(Rcpp::NumericMatrix X,Rcpp::NumericVector yy,
 Rcpp::NumericVector ww,int nt,int chosen,int draws,double nu,double lambda){
 size_t n=X.nrow();std::vector<double> xx(n*2),mean(n,0),v0(n),vr(n),obs(n);
 for(size_t k=0;k<n;++k){xx[k*2]=X(k,0);xx[k*2+1]=X(k,1);}
 xinfo xi(2);xi[0].push_back(.5);xi[1].push_back(.5);
 tree g,r;g.birth(1,0,0,.7,1.3);r.birth(1,1,0,1.2,.6);
 tree::tree_p gl=g.getptr(g.bn(xx.data()+chosen*2,xi)->nid()),rl=r.getptr(r.bn(xx.data()+chosen*2,xi)->nid());
 double initialg=gl->getm(),initialr=rl->getm();
 for(size_t k=0;k<n;++k){v0[k]=g.bn(xx.data()+k*2,xi)->getm();
  vr[k]=r.bn(xx.data()+k*2,xi)->getm();obs[k]=v0[k]*(k<(size_t)nt?vr[k]:1);}
 Rcpp::NumericMatrix out(draws,2);RNG gen;int accepted=0;
 for(int i=0;i<draws+10000;++i){
  accepted+=paired_variance_refresh(g,r,xi,xx.data(),2,chosen,yy.begin(),mean.data(),ww.begin(),
   v0.data(),vr.data(),obs.data(),n,nt,nu,lambda,gen);
  for(size_t k=0;k<n;++k){double gg=g.bn(xx.data()+k*2,xi)->getm(),rr=r.bn(xx.data()+k*2,xi)->getm();
   if(std::abs(v0[k]-gg)>1e-8 ||std::abs(vr[k]-rr)>1e-8 ||std::abs(obs[k]-gg*(k<(size_t)nt?rr:1))>1e-8)Rcpp::stop("Local cache/leaf mismatch");}
  if(i>=10000){out(i-10000,0)=log(gl->getm()/initialg);out(i-10000,1)=log(rl->getm()/initialr);}
 }
 return Rcpp::List::create(Rcpp::Named("draws")=out,Rcpp::Named("accepted")=accepted);
}
