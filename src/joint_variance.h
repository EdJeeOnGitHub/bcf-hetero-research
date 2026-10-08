#ifndef GUARD_joint_variance_h
#define GUARD_joint_variance_h
#include "tree.h"
#include "rng.h"
#include <cmath>
#include <vector>
// Symmetric proposal in log-leaf coordinates. The prior ratio includes
// their Jacobian. Proposal precision depends on fixed partition sizes only.
inline bool joint_variance_refresh(std::vector<tree>& shared,
 std::vector<tree>& ratio, const double* y,const double* mean,const double* w,
 double* v0,double* vr,double* observed,size_t n,size_t nt,
 double nu,double lambda,RNG& gen) {
 if(shared.empty() || ratio.empty() || nt==0 || nt>=n) return false;
 tree::npv g,r,bots;
 for(size_t j=0;j<shared.size();++j){bots.clear();shared[j].getbots(bots);g.insert(g.end(),bots.begin(),bots.end());}
 for(size_t j=0;j<ratio.size();++j){bots.clear();ratio[j].getbots(bots);r.insert(r.end(),bots.begin(),bots.end());}
 const double mg=shared.size(),mr=ratio.size(),a=nu/2,b=nu*lambda/2;
 double invg=0,invr=0;
 for(size_t j=0;j<g.size();++j)invg+=1/g[j]->getm();
 for(size_t j=0;j<r.size();++j)invr+=1/r[j]->getm();
 const double h11=n/2.0+a*g.size()/(mg*mg),h12=nt/2.0,
 h22=nt/2.0+a*r.size()/(mr*mr);
 const double l11=sqrt(h11),l21=h12/l11,l22=sqrt(h22-l21*l21);
 const double dr=1.7*gen.normal()/l22,dg=(1.7*gen.normal()-l21*dr)/l11;
 double rc=0,rt=0;
 for(size_t k=0;k<n;++k){double residual=y[k]-mean[k];
  double rss=w[k]*residual*residual/observed[k];if(k<nt)rt+=rss;else rc+=rss;}
 double loga=-a*g.size()*dg/mg+b*invg*(-expm1(-dg/mg))
             -a*r.size()*dr/mr+b*invr*(-expm1(-dr/mr))
             -(n-nt)*dg/2.0-nt*(dg+dr)/2.0
             +rc*(-expm1(-dg))/2.0+rt*(-expm1(-(dg+dr)))/2.0;
 if(!std::isfinite(loga) || log(gen.uniform())>=loga)return false;
 const double fg=exp(dg/mg),fr=exp(dr/mr),vg=exp(dg),vrr=exp(dr);
 for(size_t j=0;j<g.size();++j)if(!std::isfinite(g[j]->getm()*fg)||g[j]->getm()*fg<=0)return false;
 for(size_t j=0;j<r.size();++j)if(!std::isfinite(r[j]->getm()*fr)||r[j]->getm()*fr<=0)return false;
 for(size_t k=0;k<n;++k)if(!std::isfinite(v0[k]*vg)||!std::isfinite(vr[k]*vrr)||
  !std::isfinite(observed[k]*exp(dg+(k<nt?dr:0)))||
  v0[k]*vg<=0||vr[k]*vrr<=0||observed[k]*exp(dg+(k<nt?dr:0))<=0)return false;
 for(size_t j=0;j<g.size();++j)g[j]->setm(g[j]->getm()*fg);
 for(size_t j=0;j<r.size();++j)r[j]->setm(r[j]->getm()*fr);
 for(size_t k=0;k<n;++k){v0[k]*=vg;vr[k]*=vrr;observed[k]=v0[k]*(k<nt?vr[k]:1.0);}
 return true;
}
#endif
