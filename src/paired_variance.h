#ifndef GUARD_paired_variance_h
#define GUARD_paired_variance_h
#include "tree.h"
#include "rng.h"
#include <vector>
#include <cmath>
// Selection and proposal covariance depend only on covariates and partitions.
// The two inverse-gamma prior ratios include the log-coordinate Jacobians.
inline bool paired_variance_refresh(tree& g,tree& r,xinfo& xi,
 const double* xx,size_t p,size_t chosen,const double* y,const double* mean,
 const double* w,double* v0,double* vr,double* observed,size_t n,size_t nt,
 double nu,double lambda,RNG& gen) {
 tree::tree_p gl=g.getptr(g.bn(const_cast<double*>(xx+chosen*p),xi)->nid()),
              rl=r.getptr(r.bn(const_cast<double*>(xx+chosen*p),xi)->nid());
 std::vector<unsigned char> ing(n),inr(n);
 double ng=0,nr=0,ninter=0,sg=0,sr=0,sinter=0;
 for(size_t k=0;k<n;++k){
  ing[k]=g.bn(const_cast<double*>(xx+k*p),xi)==gl;
  inr[k]=r.bn(const_cast<double*>(xx+k*p),xi)==rl;
  bool tr=k<nt,gg=ing[k],rr=tr&&inr[k];
  ng+=gg;nr+=rr;ninter+=(gg&&rr);
  double resid=y[k]-mean[k],rss=w[k]*resid*resid/observed[k];
  if(gg&&rr)sinter+=rss;else if(gg)sg+=rss;else if(rr)sr+=rss;
 }
 const double a=nu/2,b=nu*lambda/2,
 h11=ng/2+a,h12=ninter/2,h22=nr/2+a,
 l11=sqrt(h11),l21=h12/l11,l22=sqrt(h22-l21*l21),
 dr=1.7*gen.normal()/l22,dg=(1.7*gen.normal()-l21*dr)/l11;
 const double loga=-a*(dg+dr)+b/gl->getm()*(-expm1(-dg))+
 b/rl->getm()*(-expm1(-dr))-(ng*dg+nr*dr)/2+
 (sg*(-expm1(-dg))+sr*(-expm1(-dr))+sinter*(-expm1(-(dg+dr))))/2;
 if(!std::isfinite(loga)||log(gen.uniform())>=loga)return false;
 const double fg=exp(dg),fr=exp(dr);
 if(!std::isfinite(gl->getm()*fg)||!std::isfinite(rl->getm()*fr)||
  gl->getm()*fg<=0||rl->getm()*fr<=0)return false;
 for(size_t k=0;k<n;++k){
  double gg=v0[k]*(ing[k]?fg:1),rr=vr[k]*(inr[k]?fr:1),obs=gg*(k<nt?rr:1);
  if(!std::isfinite(gg)||!std::isfinite(rr)||!std::isfinite(obs)||gg<=0||rr<=0||obs<=0)return false;
 }
 gl->setm(gl->getm()*fg);rl->setm(rl->getm()*fr);
 for(size_t k=0;k<n;++k){if(ing[k])v0[k]*=fg;if(inr[k])vr[k]*=fr;
  observed[k]=v0[k]*(k<nt?vr[k]:1);}
 return true;
}
#endif
