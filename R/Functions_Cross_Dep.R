#==================================================================================================>
#' Cross-sectional dependence in models
#' 
#' This routine computes all models considered for cross-sectional dependence using endogenously
#' generated weight matrices. Currently "lmcd" for linear regression, "logitcd" (logit model), 
#' "poiscd" (poisson model), "probitcd" (probit model), "nbcd" (negative binomial) are supported.
#' Note: the data should first be sorted by unit ID \code{Pid}, then by time ID \code{Tid}
#' 
#' @param pars vector of parameter values at which to evaluate the function
#' @param Y outcome variable
#' @param X matrix of covariates or design matrix
#' @param Xm matrix of other control variables
#' @param Xi matrix of variable(s) that drive cross-sectional dependence
#' @param Tid time ID
#' @param Pid unit ID
#' @param fun a function that parameterises cross-sectional dependence. it can be user-written. It must
#' take exactly three inputs: xi,xj, and k. see \link{polyexp}, \link{polyexp_mul}, \link{inter_xij}
#' for example.
#' @param k order of polynomial approximation for \code{fun()}; this can be set as a dummy for non-polynomial
#' forms of \code{fun}
#' @param nt number of terms in approximation; should be of length equal to output of \code{fun()}
#' @param utid unique time IDs to be used in the estimation; should leave at least 1 period for lagged Xi
#' @param modclass the class of model. See description above for classes supported.
#' @param rval character string naming the output to return; "crit" - value of criterion, 
#' "critobs" - observation-specific criterion values, squared residuals for 
#' "lmcd", negative log-likelihoods for MLE methods, etc., or "Xweight" vector X*W where W collects weight matrices
#' through time.
#' @param eta extra parameter to be passed, eg. \eqn{\theta} for negative binomial,\eqn{\tau} for quantile
#' regression.
#' @param dWzero logical. Should the diagonal elements of the weight matrix be set to zero? 
#' Defaults to \code{FALSE}. Choose option \code{TRUE} if own drivers \code{Xi} should affect the outcome
#' 
#' @return criterion value or vector of observation-specific criterion values if \code{crit_obs=TRUE}
#' 
#' @examples 
#' pars=c(1.0,0.5,0.8); N = 5; Tp = 6; fnp<- function(x,y,k){-(0.5*y^4+(x-y)^4)^.25} #dummy k
#' datpois=gdat_cd(pars=pars,N=N,Tp=Tp,seed=2,fun=fnp,eta = 200,modclass="poiscd") #poisson data 
#' k=1; lp=k*(k+1)/2; startp = rep(0.2,(lp+2))
#' RpacSPD::ncd_gen(pars=startp,Y=datpois$Y,X=datpois$X,Xm=NULL,Xi=datpois$X,Tid=datpois$tpID,
#' Pid=datpois$psID,fun=fnp,k=k,nt=lp,utid=c(2:Tp),modclass="poiscd") #return function value
#' RpacSPD::ncd_gen(pars=startp,Y=datpois$Y,X=datpois$X,Xm=NULL,Xi=datpois$X,Tid=datpois$tpID,
#' Pid=datpois$psID,fun=fnp,k=k,nt=lp,utid=c(2:Tp),modclass="poiscd",rval="Xw") #return Xw
#' 
#' @export

ncd_gen<- function(pars,Y,X,Xm=NULL,Xi=X,Tid,Pid,fun,k=1,nt,utid,modclass="lmcd",rval="crit",eta=0.5,dWzero=FALSE)
{
  alf = pars[1]; rho = pars[2]; 
  
  if(!is.null(Xm)){
    Xm=as.matrix(Xm)
    npxm = ncol(Xm)
    bets = pars[3:(2+npxm)]; zp = 2+npxm
  }else{
    zp = 2
  }
  
  deli = pars[-c(1:zp)];
  if(nt!=length(deli)){
    cat("The model class",modclass,"requires",I(nt+zp),"parameters","\n")
    stop("Length of parameter vector incorrect.")
  }
  
  N = length(unique(Pid)) #extract N
  ID = which(Tid>=min(utid))
  Tp = max(utid) #extract T
  XX = rep(NA,length(Y)) # to store weighted covariates
  Xi = as.matrix(Xi)
  # initialise weight matrix
  W = matrix(NA,N,N)
  for (t in utid) {
    for (i in 1:N) {
      xit1 = Xi[((i-1)*Tp + t-1),] #extract xi at t-1
      for (j in 1:N) {
        if(i!=j){ #for the diagonal term
          xjt1 = Xi[((j-1)*Tp + t-1),] #extract xj at t-1
          W[i,j] = exp(sum(c(fun(xit1,xjt1,k))*deli))
        }else{ #if diagonal, act on diagonal of W
          if(!dWzero){ 
            xjt1 = Xi[((j-1)*Tp + t-1),] #extract xj at t-1
            W[i,j] = exp(sum(c(fun(xit1,xjt1,k))*deli))
          }else{
            W[i,j] = 0.0  
          }
        }  #end if
      } #end for j
      #row-normalise
      W[i,] = W[i,]/sum(W[i,])
      id=(i-1)*Tp + t
      idintv<- ((1:N)-1)*Tp + t
      XX[id] = sum(X[idintv]*W[i,])
    }
  }
  if(!is.null(Xm)){
    lc = alf+rho*XX[ID] + Xm[ID,]%*%matrix(bets,ncol = 1)
  }else{
    lc = alf+rho*XX[ID]
  }
  if(rval=="crit"){
    val = crit(c_crit(Y[ID],lc,eta=eta,modclass = modclass))
  }else if(rval=="critobs"){
    val = crit_obs(c_crit(Y[ID],lc,eta=eta,modclass = modclass))
  }else if(rval=="Xw"){
    val = XX
  }else{
    stop(paste("The return value type",rval,"is not supported."))
  }
  return(val)
}

#==================================================================================================>
#' Huber-White Sandwich Estimator
#' 
#' A Sandwich Huber-White variance-covariance estimator.
#' 
#' @param hessmat the hessian matrix of dimension kxk evaluated at the converged parameter values
#' @param estfmat empirical estimating function
#' 
#' @return varcov a sandwich variance-covariance matrix
#' @export

swvcovHR<- function(hessmat, estfmat)
{
  estfmat = as.matrix(estfmat)
  if(nrow(estfmat)<ncol(estfmat)){estfmat=t(estfmat)} #ensure it is n x k, n>=k
  n = nrow(estfmat); k=ncol(estfmat)
  if(k!=ncol(hessmat)){"The hessian matrix and estfmat are non-conformable."}
  meat <- matrix(0,k,k)
  
  for (i in 1:n) {
    zk = matrix(estfmat[i,],nrow = 1)
    meat<- meat + t(zk)%*%zk #take sum of outer products
  } #end for i
  meat = (1/n)*meat
  bread = solve((1/n)*hessmat)
  
  varcov = bread%*%meat%*%bread
  varcov/(n-k) #adjust for degrees of freedom
}


#==================================================================================================>
#' Cross-sectional dependence in regression
#' 
#' A routine that takes starting values and model specifications as in \code{ncd_gen()} and returns
#' regression results.
#' 
#' @param startp vector of starting values
#' @param ... arguments to be passed to \code{ncd_gen()} \strong{except} argument \code{rval}. 
#' Argument names must match exactly.
#' @param rvcov Logical. Should the variance-covariance matrix be returned?
#' 
#' @return A list
#' \itemize{
#'   \item coefs vector of coefficients
#'   \item stde vector of standard errors
#'   \item tstat vector of t-statistics
#'   \item pval vector of p-values
#'   \item varcov variance-covariance matrix if \code{rvcov} is set to \code{TRUE}
#'   \item Wstat a Wald chi-square statistic
#'   \item pvwald Prob>Wstat 
#' }
#' 
#' @seealso ncd_gen
#' 
#' @examples 
#' pars=c(1.0,0.5,0.8); N = 5; Tp = 6; fnp<- function(x,y,k){-(0.5*y^4+(x-y)^4)^.25} #dummy k
#' datpois=gdat_cd(pars=pars,N=N,Tp=Tp,seed=2,fun=fnp,eta = 200,modclass="poiscd") #poisson data 
#' k=1; lp=k*(k+1)/2; startp = rep(0.2,(lp+2)) # fun() is known
#' zg1=RpacSPD::reg_cd(startp=startp,Y=datpois$Y,X=datpois$X,Xm=NULL,Xi=datpois$X,Tid=datpois$tpID,
#' Pid=datpois$psID,fun=fnp,k=k,nt=lp,utid=c(2:Tp),modclass="poiscd") #return function value
#' BIC(zg1) #compute BIC of fitted model
#'
#' k=4; lp=k*(k+1)/2; startp = rep(0.2,(lp+2)) # fun() is polynomial approximated
#' zg4=RpacSPD::reg_cd(startp=startp,Y=datpois$Y,X=datpois$X,Xm=NULL,Xi=datpois$X,Tid=datpois$tpID,
#' Pid=datpois$psID,fun=polyexp,k=k,nt=lp,utid=c(2:Tp),modclass="poiscd") #return function value
#' BIC(zg4) #compute BIC of fitted model
#' @export

reg_cd<- function(startp,...,rvcov=FALSE){
  argz = as.list(match.call())
  fn<- function(pars,rval){ncd_gen(pars = pars,...,rval=rval)}
  optobj<- stats::optim(par = startp,fn=fn,hessian = TRUE,rval="crit") #include possibility of multistart
  estfmat=pracma::jacobian(fn,optobj$par,rval="critobs")
  vcovfn=swvcovHR(hessmat=optobj$hessian, estfmat=estfmat)
  stde<- sqrt(diag(vcovfn)) # Compute standard errors
  tstat<- optobj$par/stde # Compute t statistics
  n = max(dim(estfmat)); k = min(dim(estfmat))
  df = n - k #obtain degree of freedom
  pval = 2*(1-stats::pt(abs(tstat),df))
  Wstat = matrix(optobj$par,nrow = 1)%*%solve(vcovfn)%*%matrix(optobj$par,ncol = 1) #compute wald chi-square statistic
  pvwald<- 1-stats::pchisq(Wstat,df)
  if(!rvcov){
    ans = list(coefs=optobj$par, stde=stde, tstat=tstat, pval=pval,fval=optobj$value,df=df,Wstat=Wstat,pvwald=pvwald)
  }else{
    ans = list(coefs=optobj$par, stde=stde, tstat=tstat, pval=pval,fval=optobj$value,df=df,varcov=vcovfn,Wstat=Wstat,pvwald=pvwald)
  }
  class(ans)<- c(argz$modclass,"regcd")
  ans
}

#==================================================================================================>
#' Cross-sectional dependence in regression
#' 
#' A routine like \code{reg_cd()} that partially optimises with respect \strong{\eqn{\delta}} and 
#' uses internal R routines to optimise with respect to \strong{\eqn{\beta}}. This is particularly
#' helpful in high dimensional settings.
#' 
#' @param startdel vector of starting values for \strong{\eqn{\delta}}
#' @param Y outcome variable
#' @param X matrix of covariates or design matrix
#' @param Xm matrix of other control variables
#' @param ... other arguments to be passed to \link{ncd_gen} \strong{except} arguments listed here and \code{rval}.
#' Argument names must match exactly.
#' @param modclass the class of model. See description above for classes supported.
#' @param rvcov Logical. Should the variance-covariance matrix be returned?
#' 
#' @return A list
#' \itemize{
#'   \item coefs vector of coefficients
#'   \item stde vector of standard errors
#'   \item tstat vector of t-statistics
#'   \item pval vector of p-values
#'   \item varcov variance-covariance matrix if \code{rvcov} is set to \code{TRUE}
#'   \item Wstat a Wald chi-square statistic
#'   \item pvwald Prob>Wstat 
#' }
#' 
#' @seealso \link{ncd_gen}, \link{reg_cd}
#' 
#' @examples 
#' pars = c(1.0,0.5,0.8); pars2=pars = c(1.0,0.5,0.8,0.1,-0.1); N = 10; Tp = 16 
#' fnp<- function(x,y,k) {-(0.5*y^4 + (x-y)^4)^.25} # a dummy k
#' datpois = gdat_cd(pars=pars,N=N,Tp=Tp,seed=2,fun=fnp,eta = 200,modclass="poiscd") 
#' datpois2 = gdat_cd(pars=pars2,N=N,Tp=Tp,ncXm=2,seed=2,fun=fnp,eta = 200,modclass="poiscd")
#' k=1; lp=k*(k+1)/2; startp = rep(0.2,lp); # fun() is known
#' zg1=RpacSPD::reg_cdir(startdel=startp,Y=datpois$Y,X=datpois$X,Xm=NULL,Xi=datpois$X,Tid=datpois$tpID,
#' Pid=datpois$psID,fun=fnp,k=k,nt=lp,utid=c(2:Tp),modclass="poiscd",rvcov=TRUE) #return function value
#' BIC(zg1) #compute BIC of fitted model
#'
#' k=4; lp=k*(k+1)/2; startp = rep(0,lp); # fun() is polynomial approximated
#' zg4=RpacSPD::reg_cdir(startdel=startp,Y=datpois2$Y,X=datpois2$X,Xm=datpois2[c("X1","X2")],
#' Xi=datpois2$X,Tid=datpois2$tpID,Pid=datpois2$psID,fun=polyexp,k=k,nt=lp,utid=c(2:Tp),
#' modclass="poiscd",rvcov=TRUE)
#' BIC(zg4) #compute BIC of fitted model
#' 
#' @export

reg_cdir<- function(startdel,Y,X,Xm,...,modclass="lmcd",rvcov=FALSE){
  argz = as.list(match.call())
  fncd<- function(pars,rval,eta=0.5){ncd_gen(pars = pars,Y=Y,X=X,Xm=Xm,...,modclass=modclass,rval=rval,eta=eta)}
  # assign class and pass to regir
  obdat<- list(Y=Y,Xmat=cbind(X,Xm))
  class(obdat)<- modclass 
  bets=stats::coef(regir(obdat)) #as starting values for beta
  
  fncdir<- function(del,rval="crit"){
  Xw = fncd(pars=c(bets,del),rval="Xw") #only del counts in computing Xw, bets is a dummy    
  obdat<- list(Y=Y,Xmat=cbind(Xw,Xm))
  class(obdat)<- modclass
  if(rval=="crit"){
  val = -stats::logLik(regir(obdat)) #return negative log-likelihood
  }else if(rval=="coefs"){
    if(modclass!="nbcd"){
    val = stats::coef(regir(obdat)) #return beta
    }else{
      robj = regir(obdat); th = robj$theta
      parz = stats::coef(robj) #return beta  
      val = list(pars=parz,theta = th)
    }
  }
  val
  }
  
  optobj<- stats::optim(par = startdel,fn=fncdir) #include possibility of multistart
  
  if(modclass!="nbcd"){
    bets = fncdir(del=optobj$par, rval="coefs") #extract optimal beta
    pars = c(bets,optobj$par) #merge vectors beta and delta  
    estfmat=pracma::jacobian(fncd,pars,rval = "critobs")
    hessmat=pracma::hessian(fncd,pars,rval = "crit")
  }else{
    zg = fncdir(del=optobj$par, rval="coefs") #extract optimal beta
    bets = zg$pars; th = zg$theta
    pars = c(bets,optobj$par) #merge vectors beta and delta  
    estfmat=pracma::jacobian(fncd,pars,rval = "critobs",eta=th)
    hessmat=pracma::hessian(fncd,pars,rval = "crit",eta=th)
  }
  
  vcovfn=swvcovHR(hessmat=hessmat, estfmat=estfmat)
  stde<- sqrt(diag(vcovfn)) # Compute standard errors
  tstat<- pars/stde # Compute t statistics
  n = max(dim(estfmat)); k = min(dim(estfmat))
  df = n - k #obtain degree of freedom
  pval = 2*(1-stats::pt(abs(tstat),df))
  Wstat = matrix(pars,nrow = 1)%*%solve(vcovfn)%*%matrix(pars,ncol = 1) #compute wald chi-square statistic
  pvwald<- 1-stats::pchisq(Wstat,df)
  if(!rvcov){
    ans = list(coefs=pars, stde=stde, tstat=tstat, pval=pval,fval=optobj$value,df=df,Wstat=Wstat,pvwald=pvwald)
  }else{
    ans = list(coefs=pars, stde=stde, tstat=tstat, pval=pval,fval=optobj$value,df=df,varcov=vcovfn,Wstat=Wstat,pvwald=pvwald)
  }
  class(ans)<- c(argz$modclass,"regcd")
  ans
}


#==================================================================================================>
#' Generate Data
#' 
#' Generate panel data with cross-sectional dependence for classes of models viz. Generate data for model 
#' with a single weight delta
#' 
#' @param pars vector of parameters to generate the data.
#' @param N number of units in the panel
#' @param Tp number of time periods
#' @param ncXm number of other control variables to include
#' @param seed seed for random number generation
#' @param fun function that models cross-sectional dependence function \eqn{\phi(.,.)}
#' @param sder standard deviation for the error term. this option is applicable to the class "lmcd"
#' @param eta additional parameter, like \eqn{\theta} in the negative binomial class "nbcd" and 
#' seed for simulating the outcome for MLE methods. see \link{c_out_dep}
#' @param modclass the class of model; see \link{ncd_gen} for models supported
#' 
#' @return A data frame of Y,X,psID,tpID, where psID and tpID are unit and time IDs, Xm of ncXm columns
#'  if ncXm is specified
#' 
#' @seealso \link{c_out_dep}
#' 
#' @examples 
#' pars = c(1.0,0.5,0.8); pars2=pars = c(1.0,0.5,0.8,0.1,-0.1); N = 10; Tp = 16 
#' fnp<- function(x,y) {-(0.5*y^4 + (x-y)^4)^.25} #poisson model
#' datpois = gdat_cd(pars=pars,N=N,Tp=Tp,seed=2,fun=fnp,eta = 200,modclass="poiscd") 
#' datpois2 = gdat_cd(pars=pars2,N=N,Tp=Tp,ncXm=2,seed=2,fun=fnp,eta = 200,modclass="poiscd")
#' summary(datpois); summary(datpois2);#NA's in first period for Y due to lag
#' @export

gdat_cd<- function(pars,N,Tp,ncXm=0,seed,fun,sder=NULL,eta=NULL,modclass="lmcd"){
  if(!is.null(seed)){set.seed(seed)}; 
  if(modclass=="lmcd"){eta = stats::rnorm(N*Tp,sd=sder) }
  X = stats::rnorm(N*Tp,1,1)
  Y = rep(NA,N*Tp) # vector of outcome variable
  
  alf = pars[1]; rho = pars[2]; delv=pars[3];
  if(ncXm>=1){
    bets = pars[-c(1:3)]
    if(length(bets)!=ncXm){
      stop(paste("The number of parameters must match;",I(ncXm+3),"parameters required."))
      }
    Xm = matrix(stats::rnorm(ncXm*N*Tp,0,1),ncol = ncXm)
    alf = alf + Xm%*%matrix(bets,ncol = 1) # pass as alf + linear combination of Xm and bets
    alf = c(alf)
  }

  # individual and time identifiers:
  psID = sort(rep(1:N,Tp))
  tpID = rep(1:Tp,N)
  
  for (t in 2:Tp) {
    # construct weight matrix at t
    W = matrix(NA,N,N)
    
    for (i in 1:N) {
      xit1 = X[(i-1)*Tp + t-1] #extract xi at t-1
      for (j in 1:N) {
        xjt1 = X[(j-1)*Tp + t-1] #extract xj at t-1
        W[i,j] = exp(fun(xit1,xjt1)*delv)
      }
      #row-normalise
      W[i,] = W[i,]/sum(W[i,])
      # generate outcome at time t
      id=(i-1)*Tp + t
      idintv<- (c(1:N)-1)*Tp + t
      if(!(ncXm>=1)){
        lc = alf+rho*sum(X[idintv]*W[i,])  
      }else{
        lc = alf[id]+rho*sum(X[idintv]*W[i,])
      }
      
      if(modclass=="lmcd"){
        Y[id]=out_dep(c_out_dep(lc,eta[id],modclass=modclass))
      }else{
        Y[id]=out_dep(c_out_dep(lc,eta,modclass=modclass))
      }
    }
  }
  if(!(ncXm>=1)){
  dat<- data.frame(Y,X,psID,tpID)
  }else{
    dat<- data.frame(Y,X,Xm,psID,tpID)
  }
  return(dat)
}

#=============================================================================================>
#' Bivariate Polynomial Approximation
#' 
#' Take a k-order polynomial approximation of a bivariate \eqn{\phi(xi,xj)}. The expansion excludes 
#' terms that are solely in the first variable term because coefficients on them are not identified.
#' 
#' @param xi first variable term
#' @param xj second variable term
#' @param k order of polynomial approximation
#' 
#' @return vector of length \eqn{k(k+1)/2} 
#' 
#' @examples 
#' polyexp(1.2,3.1,3)
#' polyexp(1,0,1) # an order one expansion useful for dummy variables
#' @export

polyexp<- function(xi,xj,k){
  nt = k*(k+1)/2 #number of terms in the expansion
  Z = rep(NA,nt)
  cnt = 0
  for (i in 0:(k-1)) {
    for (j in 1:(k-i)) {
      cnt = cnt+1
      Z[cnt] = (xi^i)*(xj^j)
    }
  }
  Z
}
#=============================================================================================>
#' Multivariate Polynomial Approximation
#' 
#' Take a k-order polynomial approximation of a multivariate \eqn{\phi(xi,xj)} where both xi and xj 
#' are vectors of length l. This expansion involves l pairwise polynomial expansions using 
#' \code{polyexp()}. The expansion excludes terms that are solely in the first variable term 
#' because coefficients on them are not identified.
#' 
#' @param xi vector of first variable terms
#' @param xj vector of second variable terms
#' @param k order of polynomial approximation
#' 
#' @return vector of length \eqn{l*k(k+1)/2} 
#' 
#' @examples
#' set.seed(2); xi=stats::runif(3); xj=stats::runif(3); polyexp_mul(xi,xj,2)
#' 
#' @export

polyexp_mul<- function(xi,xj,k){
  xi = c(xi); xj = c(xj)
  nz = length(xi) 
  if(nz!=length(xj)){
    stop(paste("The length of xi and xj must be equal:",nz,"!=",length(xj)))
  }
  Z = c(0)
  for (i in 1:nz) {
    Z=c(Z,polyexp(xi[i],xj[i],k))
  }
  Z[-1] 
}

#=============================================================================================>
#' Binary Interaction
#' 
#' Interact several factors. This function allows for binary \eqn{\{0,1\}} interaction between
#' observable factors xi and xj. xi and xj should be of the same length. It creates a 1 if i 
#' and j belong to the same category (xi=xj) and zero otherwise. This function also allows 
#' several factors i.e. xi and xj as vectors
#' 
#' @param xi vector of first variable terms
#' @param xj vector of second variable terms
#' @param k a dummy input as length of xi; does not need to be specified.
#' 
#' @return vector of binaries \eqn{\{0,1\}} the same length as xi
#' 
#' @examples
#' xi=c(0,1,0); xj=c(0,1,1); inter_xij(xi,xj,2)
#' 
#' @export

inter_xij<- function(xi,xj,k=length(xi)){
  xi = c(xi); xj = c(xj)
  if(k!=length(c(xi))){k = length(xi)}
  Z = c(0)
  for (i in 1:k) {
    Z[i] = ifelse(xi[i]==xj[i],1,0)
  }
  Z
}



