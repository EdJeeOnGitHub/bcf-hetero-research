#ifndef BCF_PAIRED_MEAN_H
#define BCF_PAIRED_MEAN_H
#include <map>
// Conditional Gaussian refresh of two trees at fixed partitions and scales.
inline bool paired_mean_refresh(tree& con, tree& mod, xinfo& xc, xinfo& xm,
 dinfo& dc, dinfo& dm, const double* y, const double* w,
 const double* variance, size_t ntrt, double a, double b0, double b1,
 double tc, double tm, double* fc, double* fm, double* fitted, RNG& gen) {
 tree::npv lc,lm; con.getbots(lc);mod.getbots(lm);
 const size_t nc=lc.size(),nm=lm.size(),n=dc.n;
 std::map<const tree*,size_t> ic,im;
 for(size_t j=0;j<nc;++j)ic[lc[j]]=j;
 for(size_t j=0;j<nm;++j)im[lm[j]]=nc+j;
 arma::mat q=arma::eye(nc+nm,nc+nm);
 arma::vec rhs(nc+nm,arma::fill::zeros);
 std::vector<size_t> ci(n),mi(n);
 for(size_t k=0;k<n;++k){
  ci[k]=ic.at(con.bn(dc.x+k*dc.p,xc));
  mi[k]=im.at(mod.bn(dm.x+k*dm.p,xm));
  double b=k<ntrt?b1:b0,pc=w[k]/variance[k];
  if(!std::isfinite(pc)||pc<=0||!std::isfinite(tc)||tc<=0||!std::isfinite(tm)||tm<=0)
   Rcpp::stop("Invalid paired mean conditional precision");
  double hc=a*tc,hm=b*tm;
  double r=y[k]-fitted[k]+a*lc[ci[k]]->getm()+b*lm[mi[k]-nc]->getm();
  q(ci[k],ci[k])+=pc*hc*hc;q(mi[k],mi[k])+=pc*hm*hm;
  q(ci[k],mi[k])+=pc*hc*hm;q(mi[k],ci[k])+=pc*hc*hm;
  rhs(ci[k])+=pc*hc*r;rhs(mi[k])+=pc*hm*r;
 }
 arma::mat upper;
 if(!arma::chol(upper,q))return false;
 arma::vec noise(nc+nm);
 for(size_t j=0;j<nc+nm;++j)noise(j)=gen.normal();
 arma::vec draw=arma::solve(arma::trimatu(upper),
  arma::solve(arma::trimatl(upper.t()),rhs)+noise);
 if(!draw.is_finite())Rcpp::stop("Nonfinite paired mean draw");
 for(size_t k=0;k<n;++k){
  double b=k<ntrt?b1:b0;
  double delta_c=a*(tc*draw(ci[k])-lc[ci[k]]->getm());
  double delta_m=b*(tm*draw(mi[k])-lm[mi[k]-nc]->getm());
  fc[k]+=delta_c;fm[k]+=delta_m;fitted[k]+=delta_c+delta_m;
 }
 for(size_t j=0;j<nc;++j)lc[j]->setm(tc*draw(j));
 for(size_t j=0;j<nm;++j)lm[j]->setm(tm*draw(nc+j));
 return true;
}

// Refresh the two ensemble intercept directions, holding leaf contrasts fixed.
inline void global_mean_refresh(std::vector<tree>& con,std::vector<tree>& mod,
 const double* y,const double* w,const double* variance,size_t n,size_t ntrt,
 double a,double b0,double b1,double tc,double tm,
 double* fc,double* fm,double* fitted,RNG& gen) {
 double mc=con.size(),mm=mod.size();
 arma::mat q(2,2,arma::fill::zeros);arma::vec rhs(2,arma::fill::zeros);
 for(size_t j=0;j<con.size();++j){
  tree::npv leaves;con[j].getbots(leaves);
  for(size_t l=0;l<leaves.size();++l){
   q(0,0)+=1/(mc*mc*tc*tc);rhs(0)-=leaves[l]->getm()/(mc*tc*tc);
  }
 }
 for(size_t j=0;j<mod.size();++j){
  tree::npv leaves;mod[j].getbots(leaves);
  for(size_t l=0;l<leaves.size();++l){
   q(1,1)+=1/(mm*mm*tm*tm);rhs(1)-=leaves[l]->getm()/(mm*tm*tm);
  }
 }
 for(size_t k=0;k<n;++k){
  double b=k<ntrt?b1:b0,p=w[k]/variance[k],r=y[k]-fitted[k];
  q(0,0)+=p*a*a;q(1,1)+=p*b*b;q(0,1)+=p*a*b;q(1,0)+=p*a*b;
  rhs(0)+=p*a*r;rhs(1)+=p*b*r;
 }
 arma::mat upper;if(!arma::chol(upper,q))Rcpp::stop("Global mean precision is not positive definite");
 arma::vec noise(2);noise(0)=gen.normal();noise(1)=gen.normal();
 arma::vec shift=arma::solve(arma::trimatu(upper),arma::solve(arma::trimatl(upper.t()),rhs)+noise);
 if(!shift.is_finite())Rcpp::stop("Nonfinite global mean shift");
 for(size_t j=0;j<con.size();++j){tree::npv leaves;con[j].getbots(leaves);
  for(size_t l=0;l<leaves.size();++l)leaves[l]->setm(leaves[l]->getm()+shift(0)/mc);}
 for(size_t j=0;j<mod.size();++j){tree::npv leaves;mod[j].getbots(leaves);
  for(size_t l=0;l<leaves.size();++l)leaves[l]->setm(leaves[l]->getm()+shift(1)/mm);}
 for(size_t k=0;k<n;++k){double b=k<ntrt?b1:b0;
  fc[k]+=a*shift(0);fm[k]+=b*shift(1);fitted[k]+=a*shift(0)+b*shift(1);}
}


// Refresh all mean coefficients conditionally on partitions, scales and variances.
inline size_t joint_mean_refresh(std::vector<tree>& con,std::vector<tree>& mod,
 xinfo& xc,xinfo& xm,dinfo& dc,dinfo& dm,
 const double* y,const double* w,const double* variance,size_t ntrt,
 double a,double b0,double b1,double tc,double tm,
 double* fc,double* fm,double* fitted,RNG& gen) {
 std::vector<tree*> leaves;
 std::vector<std::map<const tree*,size_t> > indices(con.size()+mod.size());
 for(size_t j=0;j<con.size()+mod.size();++j){
  tree::npv nodes;
  if(j<con.size())con[j].getbots(nodes);else mod[j-con.size()].getbots(nodes);
  for(size_t l=0;l<nodes.size();++l){indices[j][nodes[l]]=leaves.size();leaves.push_back(nodes[l]);}
 }
 size_t nc=0;for(size_t j=0;j<con.size();++j)nc+=indices[j].size();
 const size_t n=dc.n,L=leaves.size();
 arma::mat h(n,L,arma::fill::zeros);
 arma::vec residual(n),precision(n);
 for(size_t k=0;k<n;++k){
  double b=k<ntrt?b1:b0;
  precision(k)=w[k]/variance[k];
  if(!std::isfinite(precision(k))||precision(k)<=0)Rcpp::stop("Invalid joint mean precision");
  residual(k)=y[k]-fitted[k]+fc[k]+fm[k];
  for(size_t j=0;j<con.size();++j)h(k,indices[j].at(con[j].bn(dc.x+k*dc.p,xc)))=a*tc;
  for(size_t j=0;j<mod.size();++j)h(k,indices[con.size()+j].at(mod[j].bn(dm.x+k*dm.p,xm)))=b*tm;
 }
 arma::mat weighted=h;weighted.each_col()%=arma::sqrt(precision);
 arma::mat q=weighted.t()*weighted;q.diag()+=1;
 arma::vec rhs=h.t()*(precision%residual),noise(L);
 arma::mat upper;if(!arma::chol(upper,q))Rcpp::stop("Joint mean precision is not positive definite");
 for(size_t j=0;j<L;++j)noise(j)=gen.normal();
 arma::vec draw=arma::solve(arma::trimatu(upper),arma::solve(arma::trimatl(upper.t()),rhs)+noise);
 if(!draw.is_finite())Rcpp::stop("Nonfinite joint mean draw");
 arma::vec new_c=h.cols(0,nc-1)*draw.rows(0,nc-1);
 arma::vec new_m=h.cols(nc,L-1)*draw.rows(nc,L-1);
 for(size_t k=0;k<n;++k){fitted[k]+=new_c(k)+new_m(k)-fc[k]-fm[k];fc[k]=new_c(k);fm[k]=new_m(k);}
 for(size_t j=0;j<L;++j)leaves[j]->setm(draw(j)*(j<nc?tc:tm));
 return L;
}
#endif
