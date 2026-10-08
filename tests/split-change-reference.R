library(bcf);library(RcppParallel)
RcppParallel::setThreadOptions(numThreads=1L)
args=commandArgs(TRUE);stopifnot(length(args)<=1L)
if(length(args)) {
 src=normalizePath(args[1],mustWork=TRUE)
} else {
 candidates=c('.', '..', '../00_pkg_src/bcf')
 roots=candidates[file.exists(file.path(candidates,'src','varfuns.h'))]
 stopifnot(length(roots)>0L)
 src=normalizePath(roots[1],mustWork=TRUE)
}
Sys.setenv(PKG_CPPFLAGS=paste('-I',shQuote(file.path(src,'src')),sep=''),
 PKG_LIBS=shQuote(system.file('libs','bcf.so',package='bcf')))
# R CMD check sets a relative startup file; child compilation changes directory.
Sys.unsetenv("R_TESTS")
Rcpp::sourceCpp(file.path(src,'tests/split-change-reference.cpp'))
set.seed(42);n=30L
X=cbind(seq_len(n),sample(seq_len(n)))
y=.3*sin(seq_len(n))+rnorm(n,sd=.5);phi=rep(c(.3,1,4),length.out=n)
check=function(X,cuts,ancestor=-1L){
 p=ncol(X);tau=.5;depth=if(ancestor<0)0 else 1
 region=if(ancestor<0)rep(TRUE,n)else X[,1]<cuts[[1]][ancestor+1L]
 allowed=lapply(seq_len(p),function(v)if(v==1 && ancestor>=0)seq_len(ancestor)else seq_along(cuts[[v]]))
 mass=numeric(sum(lengths(cuts)));offset=c(0,cumsum(lengths(cuts)))
 for(v in seq_len(p))for(c in allowed[[v]]){
  left=region & X[,v]<cuts[[v]][c];right=region & !left
  if(sum(left)<5 || sum(right)<5)next
  # Independent Gaussian marginal via the full covariance of observations.
  Z=cbind(left[region],right[region])*1
  V=diag(1/phi[region])+tau^2*tcrossprod(Z)
  likelihood=-.5*(as.numeric(determinant(V,logarithm=TRUE)$modulus)+
                     sum(y[region]*solve(V,y[region])))
  terminal=function(side){
   available=any(vapply(seq_len(p),function(j){
    a=allowed[[j]]
    if(j==v)a=if(side=='left')a[a<c]else a[a>c]
    length(a)>0
   },logical(1)))
   if(available)1-.95/(depth+2)^2 else 1
  }
  mass[offset[v]+c]=exp(likelihood)*terminal('left')*terminal('right')/
                         (p*length(allowed[[v]]))
 }
 expected=mass/sum(mass)
 draw=split_change_reference(X,y,phi,cuts,120000L,ancestor,tau)
 observed=draw$counts/sum(draw$counts)
 stopifnot(draw$accepted>1000,max(abs(expected-observed))<.015)
 cat('SPLIT_CHANGE_EXACT_REFERENCE_PASS',p,ancestor,max(abs(expected-observed)),'\n')
}
check(X[,1,drop=FALSE],list(c(6.5,12.5,18.5,24.5)))
check(X,list(c(6.5,12.5,18.5,24.5),c(10.5,20.5)))
check(X,list(c(6.5,12.5,18.5,24.5),c(10.5,20.5)),ancestor=3L)
