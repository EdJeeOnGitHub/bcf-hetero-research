library(bcf)
stopifnot(as.character(packageVersion('bcf'))=='2.0.2.9011')
set.seed(422);n=60L;x=matrix(rnorm(n*2),n);z=rep(0:1,each=n/2);y=rnorm(n)+z
for(hc in c(FALSE,TRUE)) {
 td=tempfile('shared-sampler-');dir.create(td)
 f=bcf(y=y,z=z,x_control=x,x_moderate=x,x_variance=x,pihat=rep(.5,n),w=rep(1,n),
 nburn=20,nsim=30,n_chains=1,n_threads=1,ntree_control=5,ntree_moderate=3,
 vartree=list(num_trees=3,nu=5,lambda=1,numcut=10),variance_model='ratio',
 use_muscale=hc,use_tauscale=FALSE,joint_mean_every=5L,variance_split_change=TRUE,
 joint_variance_every=1L,paired_variance_every=1L,save_tree_directory=td,verbose=FALSE)
 stopifnot(f$joint_mean_updates==10L,f$joint_variance_attempts==50L,
 f$paired_variance_attempts==150L,f$joint_variance_accepts>0,f$paired_variance_accepts>0)
 p=predict(f,x_predict_control=x,x_predict_moderate=x,x_predict_variance=x,
 pi_pred=rep(.5,n),z_pred=z,save_tree_directory=td,n_cores=1,verbose=FALSE)
 stopifnot(max(abs(p$sigma0_2-f$sigma0_2))<1e-8,
 max(abs(p$sigma1_2-f$sigma1_2))<1e-8,max(abs(p$tau-f$tau))<1e-8)
 cat('SHARED_SAMPLER_REPLAY_PASS mu_half_cauchy=',hc,'\n')
}
