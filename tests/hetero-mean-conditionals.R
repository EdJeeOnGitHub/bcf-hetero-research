# Analytical checks of active mean updates with known, unequal arm variances.
# A concentrated variance prior pins sigma0^2=2 and sigma1^2=4. Constant
# predictors pin every tree to a stump, giving exact Gaussian / horseshoe
# posterior references independent of the sampler implementation.
suppressPackageStartupMessages({library(bcf);library(coda)})
stopifnot(as.character(packageVersion('bcf'))==Sys.getenv('BCF_EXPECT_VERSION','2.0.2.9011'))
paired=identical(Sys.getenv('BCF_PAIRED_MEAN_UPDATE'),'true')
set.seed(831);n=60L;z=rep(0:1,each=n/2);w=rep(c(.25,1,4),length.out=n)
v=ifelse(z==0,2,4);y=1.1+.6*z+rnorm(n)*sqrt(v/w)
x=matrix(0,n,1);root=tempfile('mean-conditionals-');dir.create(root)
fit_case=function(name,hc,moderate_sd){
 p=file.path(root,name);dir.create(p)
 set.seed(992)
 bcf(y,z,x,x,rep(.5,n),w=w,standardize=FALSE,
  collapsed_mu_scale=hc && identical(Sys.getenv("BCF_COLLAPSED_MU_SCALE"),"true"),
  mean_split_change=identical(Sys.getenv("BCF_MEAN_SPLIT_CHANGE"),"true"),
  joint_mean_every=as.integer(Sys.getenv("BCF_JOINT_MEAN_EVERY","0")),
  global_mean_update=identical(Sys.getenv("BCF_GLOBAL_MEAN_UPDATE"),"true"), paired_mean_update=paired, sd_control=1,sd_moderate=moderate_sd,use_muscale=hc,use_tauscale=FALSE,
  ntree_control=3,ntree_moderate=2,nburn=1000,nsim=6000,
  n_chains=1,n_threads=1,x_variance=x,variance_model='ratio',
  vartree=list(num_trees=2,nu=1e9,lambda=2,numcut=10,sparse=FALSE),
  save_tree_directory=p,log_file=file.path(p,'fit.log'),verbose=FALSE)
}
g=fit_case('gaussian',FALSE,.7)
stopifnot(abs(mean(g$sigma0_2)-2)<.001,abs(mean(g$sigma1_2)-4)<.002)
design=cbind(1,z-.5);precision=w/v
vc=solve(diag(c(1,1/.7^2))+crossprod(design,design*precision))
mc=as.numeric(vc%*%crossprod(design,y*precision))
projection=rbind(c(1,-.5),c(0,1))
reference_mean=as.numeric(projection%*%mc)
reference_cov=projection%*%vc%*%t(projection)
draws=cbind(mu0=g$mu[,1],tau=g$tau[,1])
neff=as.numeric(effectiveSize(draws))
mean_limit=6*sqrt(diag(reference_cov)/neff)+.001
cov_limit=6*sqrt((outer(diag(reference_cov),diag(reference_cov))+reference_cov^2)/min(neff))+.001
print(data.frame(component=colnames(draws),reference=reference_mean,
 empirical=colMeans(draws),ess=neff,tolerance=mean_limit))
stopifnot(all(abs(colMeans(draws)-reference_mean)<mean_limit),
 all(abs(cov(draws)-reference_cov)<cov_limit))
rm(g);gc(FALSE)
h=fit_case('half-cauchy',TRUE,1e-8)
p=sum(precision);likelihood_variance=1/p;ybar=sum(precision*y)/p
# lambda~half-Cauchy(1), mu|lambda~N(0,lambda^2). Under lambda=tan(t),
# its measure is constant on (0,pi/2), so the constant cancels below.
conditional=function(t){
 prior_variance=tan(t)^2;total=prior_variance+likelihood_variance
 list(weight=dnorm(ybar,0,sqrt(total)),
  mean=prior_variance/total*ybar,
  variance=prior_variance/total*likelihood_variance)
}
quad=function(f)integrate(f,0,pi/2,rel.tol=1e-9,subdivisions=500)$value
normalizer=quad(function(t)conditional(t)$weight)
ref_mean=quad(function(t){a=conditional(t);a$weight*a$mean})/normalizer
ref_var=quad(function(t){a=conditional(t);a$weight*(a$variance+(a$mean-ref_mean)^2)})/normalizer
ref_fourth=quad(function(t){a=conditional(t);b=a$mean-ref_mean
 a$weight*(3*a$variance^2+6*a$variance*b^2+b^4)})/normalizer
u=h$mu[,1];ess_mean=as.numeric(effectiveSize(u))
ess_variance=as.numeric(effectiveSize((u-ref_mean)^2))
mean_limit=6*sqrt(ref_var/ess_mean)+.001
var_limit=6*sqrt((ref_fourth-ref_var^2)/ess_variance)+.001
print(data.frame(reference_mean=ref_mean,empirical_mean=mean(u),
 reference_variance=ref_var,empirical_variance=var(u),ess_mean=ess_mean,
 mean_tolerance=mean_limit,variance_tolerance=var_limit))
stopifnot(abs(mean(u)-ref_mean)<mean_limit,abs(var(u)-ref_var)<var_limit)
cat('ACTIVE_HETERO_MEAN_CONDITIONALS_PASS\n')
