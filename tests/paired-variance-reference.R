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
Rcpp::sourceCpp(file.path(src,'tests/paired-variance-reference.cpp'))
set.seed(192);n=16L;nt=8L
X=cbind(rep(c(0,0,1,1),4),rep(c(0,1,0,1),4));y=rnorm(n)
nu=5;lambda=.6;grid=seq(-8,8,length.out=801)
for(chosen in c(1L,2L))for(w in list(rep(1,n),rep(c(.25,3,1,.5),4))){
 ing=X[,1]==X[chosen,1];inr=X[,2]==X[chosen,2]
 baseg=ifelse(X[,1]<.5,.7,1.3);baser=ifelse(X[,2]<.5,1.2,.6)
 a=nu/2;b=nu*lambda/2
 prior=function(eta,base)-a*(log(base)+eta)-b/(base*exp(eta))
 mat=outer(prior(grid,baseg[chosen]),prior(grid,baser[chosen]),'+')
 for(k in seq_len(n)){
  lg=log(baseg[k])+if(ing[k])grid else rep(0,length(grid))
  lr=if(k<=nt)log(baser[k])+if(inr[k])grid else rep(0,length(grid)) else rep(0,length(grid))
  lv=outer(lg,lr,'+');mat=mat-.5*lv-.5*w[k]*y[k]^2*exp(-lv)
 }
 mass=exp(mat-max(mat));mass=mass/sum(mass)
 stopifnot(sum(mass[c(1,nrow(mass)),])<1e-8,sum(mass[,c(1,ncol(mass))])<1e-8)
 marg=list(rowSums(mass),colSums(mass))
 m1=vapply(marg,function(p)sum(grid*p),numeric(1))
 m2=vapply(marg,function(p)sum(grid^2*p),numeric(1))
 q=sapply(marg,function(p)approx(cumsum(p),grid,xout=c(.1,.5,.9),ties='ordered')$y)
 f=paired_variance_reference(X,y,w,nt,chosen-1L,200000L,nu,lambda)
 stopifnot(f$accepted>10000,max(abs(colMeans(f$draws)-m1))<.025,
  max(abs(colMeans(f$draws^2)-m2))<.06,
  max(abs(apply(f$draws,2,quantile,probs=c(.1,.5,.9))-q))<.04)
 cat('PAIRED_VARIANCE_NUMERICAL_REFERENCE_PASS',chosen,paste(w,collapse=','),'\n')
}
# Actual forest integration with both global and local moves enabled.
set.seed(93);n=60L;x=matrix(rnorm(n*2),n);z=rep(0:1,each=n/2)
td=file.path(tempdir(),'paired-variance-replay');dir.create(td)
f=bcf(y=rnorm(n)+z,z=z,x_control=x,x_moderate=x,x_variance=x,pihat=rep(.5,n),
 w=rep(1,n),nburn=10,nsim=10,n_chains=1,n_threads=1,ntree_control=5,ntree_moderate=3,
 vartree=list(num_trees=3,nu=5,lambda=1,numcut=10),variance_model='ratio',
 use_muscale=FALSE,use_tauscale=FALSE,joint_variance_every=1L,paired_variance_every=1L,
 save_tree_directory=td,verbose=FALSE)
stopifnot(f$paired_variance_attempts==60,f$paired_variance_accepts>0,f$joint_variance_attempts==20)
p=predict(f,x_predict_control=x,x_predict_moderate=x,x_predict_variance=x,
 pi_pred=rep(.5,n),z_pred=z,save_tree_directory=td,n_cores=1,verbose=FALSE)
stopifnot(max(abs(p$sigma0_2-f$sigma0_2))<1e-8,
 max(abs(p$sigma1_2-f$sigma1_2))<1e-8,max(abs(p$tau-f$tau))<1e-8)
cat('PAIRED_VARIANCE_FULL_REPLAY_PASS\n')
