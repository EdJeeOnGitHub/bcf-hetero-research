library(bcf)
args=commandArgs(TRUE);stopifnot(length(args)<=1L)
if(length(args)) {
 src=normalizePath(args[1],mustWork=TRUE)
} else {
 candidates=c('.', '..', '../00_pkg_src/bcf')
 roots=candidates[file.exists(file.path(candidates,'src','varfuns.h'))]
 stopifnot(length(roots)>0L)
 src=normalizePath(roots[1],mustWork=TRUE)
}
Sys.setenv(PKG_CPPFLAGS=paste('-I',shQuote(file.path(src,'src'))),
 PKG_LIBS=shQuote(system.file('libs','bcf.so',package='bcf')))
# R CMD check sets a relative startup file; child compilation changes directory.
Sys.unsetenv("R_TESTS")
Rcpp::sourceCpp(file.path(src,'tests/joint-variance-reference.cpp'))
set.seed(891)
# Independent numerical posterior over the two log-scale shifts. Each forest
# has two unequal stump leaves; the contrasts stay fixed under this kernel.
y=c(.6,-1.4,.8,1.1,-.5,.2);nt=4L;nu=5;lambda=.6
for(w in list(rep(1,6),c(.25,3,1,.5,4,2))){
 grid=seq(-8,8,length.out=801)
 a=nu/2;b=nu*lambda/2
 lp=function(eta,base)sum(-a*(log(base)+eta/2)-b/(base*exp(eta/2)))
 pg=vapply(grid,lp,numeric(1),base=c(.7,1.3))
 pr=vapply(grid,lp,numeric(1),base=c(1.2,.6))
 mat=outer(pg,pr,'+')
 for(k in seq_along(y)){
  lv=outer(log(.91)+grid,if(k<=nt)log(.72)+grid else rep(0,length(grid)),'+')
  mat=mat-.5*lv-.5*w[k]*y[k]^2*exp(-lv)
 }
 mass=exp(mat-max(mat));mass=mass/sum(mass)
 stopifnot(sum(mass[c(1,nrow(mass)),])<1e-8,sum(mass[,c(1,ncol(mass))])<1e-8)
 marginal=list(rowSums(mass),colSums(mass))
 refmean=vapply(marginal,function(p)sum(grid*p),numeric(1))
 refsecond=vapply(marginal,function(p)sum(grid^2*p),numeric(1))
 refq=sapply(marginal,function(p)approx(cumsum(p),grid,xout=c(.1,.5,.9),ties='ordered')$y)
 fit=joint_variance_reference(y,w,nt,200000L,nu,lambda)
 stopifnot(fit$accepted>10000,
  max(abs(colMeans(fit$draws)-refmean))<.025,
  max(abs(colMeans(fit$draws^2)-refsecond))<.06,
  max(abs(apply(fit$draws,2,quantile,probs=c(.1,.5,.9))-refq))<.04)
 cat('JOINT_VARIANCE_NUMERICAL_REFERENCE_PASS',paste(w,collapse=','),'\n')
}
# Full ratio model checks that serialization/prediction sees the updated trees.
set.seed(91);n=60L;x=matrix(rnorm(n*2),n);z=rep(0:1,each=n/2)
td=file.path(tempdir(),'joint-variance-replay');dir.create(td)
f=bcf(y=rnorm(n)+z,z=z,x_control=x,x_moderate=x,x_variance=x,pihat=rep(.5,n),
 w=rep(1,n),nburn=10,nsim=10,n_chains=1,n_threads=1,ntree_control=5,ntree_moderate=3,
 vartree=list(num_trees=3,nu=5,lambda=1,numcut=10),variance_model='ratio',
 use_muscale=FALSE,use_tauscale=FALSE,joint_variance_every=1L,
 save_tree_directory=td,verbose=FALSE)
stopifnot(f$joint_variance_attempts==20,f$joint_variance_accepts>0)
p=predict(f,x_predict_control=x,x_predict_moderate=x,x_predict_variance=x,
 pi_pred=rep(.5,n),z_pred=z,save_tree_directory=td,n_cores=1,verbose=FALSE)
stopifnot(max(abs(p$sigma0_2-f$sigma0_2))<1e-8,
 max(abs(p$sigma1_2-f$sigma1_2))<1e-8,max(abs(p$tau-f$tau))<1e-8)
cat('JOINT_VARIANCE_FULL_REPLAY_PASS\n')
