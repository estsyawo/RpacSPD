#' Cross-sectional dependence in models
#' 
#' This routine computes all models considered for cross-sectional dependence using endogenously
#' generated weight matrices. Currently "lm" for linear regression, "logit" (logit model), 
#' "poisson" (poisson model), "probit" (probit model), "negbin" (negative binomial) are supported.
#' Note: the data should first be sorted by unit ID \code{Pid}, then by time ID \code{Tid}
#' 
#' @param pars vector of parameter values at which to evaluate the function
#' @param Y outcome variable
#' @param X matrix of covariates or design matrix
#' @param Xm matrix of other control variables
#' @param Xi matrix of variable(s) that drive cross-sectional dependence
#' @param Tid time ID
#' @param Pid unit ID
#' @param fun a function that parameterises the form of dependence. it can be user-written
#' @param k order of polynomial approximation for \code{fun()}
#' @param nt number of terms in approximation; should be of length equal to output of \code{fun()}
#' @param utid unique time IDs to be used in the estimation; should leave at least 1 period for lagged Xi
#' @param modclass the class of model. See description above for classes supported.
#' @param crit_obs a logical for vector of observation-specific criterion values, squared residuals for 
#' "lm", negative log-likelihoods for MLE methods, etc.
#' 
#' @return criterion value or vector of observation-specific criterion values if \code{crit_obs=TRUE}
#' 
#' @export

ncd_gen<- function(pars,Y,X,Xm=NULL,Xi=X,Tid,Pid,fun=fn4,k=1,nt,utid,modclass="lm",crit_obs=FALSE)
{
  if(!is.null(Xm)){npxm = ncol(as.matrix(Xm))}else{npxm=0}
  
  alf = pars[1]; rho = pars[2]; bets = pars[3:(2+npxm)]
  deli = pars[-c(1:(2+npxm))];
  if(nt!=length(deli)){stop("Length of parameter vector incorrect.")}
  
  N = length(unique(Pid)) #extract N
  ID = which(Tid>=min(utid))
  Tp = max(utid) #extract T
  XX = rep(NA,length(Y)) # to store weighted covariates
  Xi = as.matrix(Xi)
  
  for (t in utid) {
    # construct weight matrix at t
    W = matrix(NA,N,N)
    for (i in 1:N) {
      xit1 = Xi[((i-1)*Tp + t-1),] #extract xi at t-1
      for (j in 1:N) {
        xjt1 = Xi[((j-1)*Tp + t-1),] #extract xj at t-1
        W[i,j] = exp(sum(c(fun(xit1,xjt1,k))*deli))
      }
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
  if(!crit_obs){
    val = crit(c_crit(Y[ID],lc,eta=eta,klass = modclass))
  }else{
    val = crit_obs(c_crit(Y[ID],lc,eta=eta,klass = modclass))
  }
  return(val)
}
