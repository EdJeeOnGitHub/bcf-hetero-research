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
Rcpp::sourceCpp(file.path(src,'tests/variance-change-reference.cpp'))
set.seed(42);n=30L
X=cbind(seq_len(n),sample(seq_len(n)))
y=.3*sin(seq_len(n))+rnorm(n,sd=.5);phi=rep(c(.3,1,4),length.out=n)
check=function(X,cuts,ancestor=-1L){
 p=ncol(X);nu=8;lambda=.6;depth=if(ancestor<0)0 else 1
 region=if(ancestor<0)rep(TRUE,n)else X[,1]<cuts[[1]][ancestor+1L]
 allowed=lapply(seq_len(p),function(v)if(v==1 && ancestor>=0)seq_len(ancestor)else seq_along(cuts[[v]]))
 mass=numeric(sum(lengths(cuts)));offset=c(0,cumsum(lengths(cuts)))
 for(v in seq_len(p))for(c in allowed[[v]]){
  left=region & X[,v]<cuts[[v]][c];right=region & !left
  if(sum(left)<5 || sum(right)<5)next
  # Independently integrate the Gaussian likelihood against an IG prior.
  marginal=function(ix){
   a=nu/2;b=nu*lambda/2;yy=y[ix];pp=phi[ix]
   mode=(b+sum(pp*yy^2)/2)/(a+length(yy)/2)
   logf=function(t){
    v=exp(t);out=rep(-Inf,length(t));valid=is.finite(v)&v>0
    out[valid]=vapply(which(valid),function(k){
     sum(dnorm(yy,0,sqrt(v[k]/pp),log=TRUE))+a*log(b)-lgamma(a)-
       (a+1)*t[k]-b/v[k]+t[k]
    },numeric(1));out
   }
   peak=logf(log(mode))
   peak+log(integrate(function(t)exp(logf(t)-peak),-Inf,Inf,
                      rel.tol=1e-9,abs.tol=0)$value)
  }
  likelihood=marginal(left)+marginal(right)
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
 draw=variance_change_reference(X,y,phi,cuts,120000L,ancestor,nu,lambda)
 observed=draw$counts/sum(draw$counts)
 stopifnot(draw$accepted>1000,max(abs(expected-observed))<.015)
 cat('VARIANCE_CHANGE_EXACT_REFERENCE_PASS',p,ancestor,max(abs(expected-observed)),'\n')
}
check(X[,1,drop=FALSE],list(c(6.5,12.5,18.5,24.5)))
check(X,list(c(6.5,12.5,18.5,24.5),c(10.5,20.5)))
check(X,list(c(6.5,12.5,18.5,24.5),c(10.5,20.5)),ancestor=3L)
