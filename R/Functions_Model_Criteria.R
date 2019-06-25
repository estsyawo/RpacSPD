#' Model criterion function
#' 
#' A generic S3 function for evaluating the criteria of models, eg. sum of squares for linear 
#' regression, negative sum of the log of individual likelihoods for MLE, etc.
#' 
#' @param object the object to be passed to the concrete class constructor \code{c_crit()}
#' 
#' @export

crit<- function(object) UseMethod("crit")

#==================================================================================================>
#' Model criterion function
#' 
#' A generic S3 function for evaluating the criteria of models, eg. squares of residuals for linear 
#' regression, negative log of individual likelihoods for MLE, etc.
#' 
#' @param object the object to be passed to the concrete class constructor \code{c_crit()}
#' 
#' @export
crit_obs<- function(object) UseMethod("crit_obs")

##==================================================================================================>
#' Concrete class constructor
#' 
#' A function for constructing functions for concrete classes of models for the \code{crit()} family of
#'  of functions.
#' 
#' @param Y vector of the outcome variable
#' @param lc linear combination \eqn{lc = X'\beta}
#' @param eta for models needing extra parameter besides beta in the linear combination \code{lc}. 
#' eta can be a vector
#' 
#' @param modclass the class of model. Currently "lmcd" for linear regression, "logitcd" (logit model), 
#' "poiscd" (poisson model), "probitcd" (probit model), "nbcd" (negative binomial) are supported.
#' 
#' @return object an object list with class attribute modclass.
#' 
#' @export

c_crit<- function(Y, lc, eta=NULL, modclass="lmcd"){ #assign class of object
  if(is.null(eta)){ 
    object<- list(Y,lc)
    class(object)<- modclass
    names(object)<- c("Y","lc")  
  }else{ 
    object<- list(Y,lc,eta)
    class(object)<- modclass
    names(object)<- c("Y","lc","eta")  
  }
  object
}

#==================================================================================================>
#' Implemention - linear regression
#' 
#' This function evaluates the criterion for the linear regression class of model. 
#' 
#' @param object Object with class attribute "lmcd" obtained from the function \code{c_crit()}
#' @return sum of squared errors
#' 
#' @export

crit.lmcd<- function(object){sum((object$Y-object$lc)^2)} #for least squares

#==================================================================================================>
#' Implemention - logit model
#' 
#' This function evaluates the criterion for the logit class of model. 
#' 
#' @param object Object with class attribute "logitcd" obtained from the function \code{c_crit()}
#' @return negative log likelihood
#' 
#' @export

crit.logitcd<- function(object) {-sum(object$Y*object$lc - log(1+exp(object$lc))) } #return neg. log-like for logit

#==================================================================================================>
#' Implemention - poisson model
#' 
#' This function evaluates the criterion for the poisson class of model. 
#' 
#' @param object Object with class attribute "poiscd" obtained from the function \code{c_crit()}
#' @return negative log likelihood
#' 
#' @export

crit.poiscd<- function(object) {-sum(object$Y*object$lc - exp(object$lc))} #return neg. log-like for poisson

#==================================================================================================>
#' Implemention - probit model
#' 
#' This function evaluates the criterion for the probit class of model. 
#' 
#' @param object Object with class attribute "probitcd" obtained from the function \code{c_crit()}
#' @return negative log likelihood
#' 
#' @export

crit.probitcd<- function(object){
  plc = stats::pnorm(object$lc)
  lk = -sum(object$Y*log(plc)+(1-object$Y)*log(1-plc))
  lk
} 

#==================================================================================================>
#' Implemention - negative binomial model
#' 
#' This function evaluates the criterion for the negative binomial class of model. 
#' 
#' @param object Object with class attribute "nbcd" obtained from the function \code{c_crit()}
#' @return negative log likelihood
#' 
#' @export

crit.nbcd<- function(object){
  if(is.null(object$eta)){
    stop("Error: The negative binomial class requires an additional parameter eta")
  } #end if
  theta = object$eta; Y = object$Y; mu = exp(object$lc); 
  lv = lgamma(theta+Y)-lgamma(theta)-lgamma(Y+1)+theta*log(theta)+Y*log(mu+(Y == 0))-(theta+Y)*log(theta+mu)
  -sum(lv)
} #return the neg. log-like for negative binomial

#==================================================================================================>
#' Implemention - linear regression model
#' 
#' This function evaluates to a vector of squared residuals for the linear regression class of model. 
#' 
#' @param object Object with class attribute "lmcd" obtained from the function \code{c_crit()}
#' @return vector of squared residuals
#' 
#' @export

crit_obs.lmcd<- function(object){c((object$Y-object$lc)^2)} #for least squares

#==================================================================================================>
#' Implemention - logit model
#' 
#' This function evaluates to a vector of log-likelihoods for the logit class of model. 
#' 
#' @param object Object with class attribute "logitcd" obtained from the function \code{c_crit()}
#' @return vector of negative log likelihoods
#' 
#' @export

crit_obs.logitcd<- function(object) {-c(object$Y*object$lc - log(1+exp(object$lc))) } #return neg. log-like for logit

#==================================================================================================>
#' Implemention - poisson model
#' 
#' This function evaluates to a vector of log-likelihoods for the poisson class of model. 
#' 
#' @param object Object with class attribute "poiscd" obtained from the function \code{c_crit()}
#' @return vector of negative log likelihoods
#' 
#' @export

crit_obs.poiscd<- function(object) {-c(object$Y*object$lc - exp(object$lc))} #return neg. log-like for poisson

#==================================================================================================>
#' Implemention - probit model
#' 
#' This function evaluates to a vector of log-likelihoods for the probit class of model. 
#' 
#' @param object Object with class attribute "probitcd" obtained from the function \code{c_crit()}
#' @return vector of negative log likelihoods
#' 
#' @export

crit_obs.probitcd<- function(object){
  plc = stats::pnorm(object$lc)
  lk = -c(object$Y*log(plc)+(1-object$Y)*log(1-plc))
  lk
}

#==================================================================================================>
#' Implemention - negative binomial model
#' 
#' This function evaluates to a vector of log-likelihoods for the negative binomial class of model. 
#' 
#' @param object Object with class attribute "nbcd" obtained from the function \code{c_crit()}
#' @return vector of negative log likelihoods
#' 
#' @export

crit_obs.nbcd<- function(object){
  if(is.null(object$eta)){
    stop("Error: This negative binomial class requires an additional parameter eta")
  } #end if
  theta = object$eta; Y = object$Y; mu = exp(object$lc); 
  lv = lgamma(theta+Y)-lgamma(theta)-lgamma(Y+1)+theta*log(theta)+Y*log(mu+(Y == 0))-(theta+Y)*log(theta+mu)
  #lv = eta2*log(eta2/dn) + object$Y*log(mu/dn) + log(gamma(object$Y+eta2)/gamma(eta2))
  -c(lv)
}

#=============================================================================================>
#' Randomly generate outcome Y
#' 
#' A generic function for randomly generating outcome variables for model classes.
#' 
#' @param object Object of class attribute obtained from \code{c_out_dep()}
#'  
#' @export

out_dep<- function(object) UseMethod("out_dep")

#=============================================================================================>
#' Object class constructor
#' 
#' This function constructs an object class for the generic function \code{out_dep()}
#' 
#' @param lc linear combination `lc = X'beta`
#' @param eta for models needing extra parameter besides beta in the linear combination \code{lc}. 
#' eta can be a vector. For "lmcd", \code{eta} is the vector of errors; must of the same length as Y
#' For "nbcd", \code{eta} must be a vector of two elements; the first one is parameter \code{th}
#' and the second is the seed for random generation. For other classes, eta is the seed for random 
#' generation.
#' 
#' @param modclass the class of model. Currently "lmcd" for linear regression, "logitcd" (logit model), 
#' "poiscd" (poisson model), "probitcd" (probit model), "nbcd" (negative binomial) are supported.
#' 
#' @return object an object list with class attribute modclass.
#' 
#' @export

c_out_dep<- function(lc,eta=NULL,modclass="lmcd"){ #assign class object
  if(is.null(eta)){ 
    object<- lc
    class(object)<- modclass
    names(object)<- c("lc")  
  }else{ 
    object<- list(lc,eta)
    class(object)<- modclass
    names(object)<- c("lc","eta")  
  }
  object
}

#=============================================================================================>
#' Implemention - linear regression model
#' 
#' Generate outcome variable for linear regression class of model
#' 
#' @param object object of class "lmcd"
#' 
#' @return Y vector of outcome variable

out_dep.lmcd<- function(object){
  if(length(object$lc)==length(object$eta)){
    Y=object$lc + object$eta #eta component is the error term. must be of length lc
  }else{
    stop("lc and error terms eta must be of the same length.")
  }
  Y
}

#=============================================================================================>
#' Implemention - logit model
#' 
#' Generate outcome variable for logit regression class of model
#' 
#' @param object object of class "logitcd"
#' 
#' @return Y vector of outcome variable
#' @export

out_dep.logitcd<- function(object){
  probs= stats::plogis(object$lc)
  N = length(probs); Y = rep(0,N)
  set.seed(object$eta) #seed is passed as eta
  for (i in 1:N) {
    Y[i]=sample(c(1,0),1,replace = TRUE,prob = c(probs[i],(1-probs[i])))
  }
  Y
}

#=============================================================================================>
#' Implemention - probit model
#' 
#' Generate outcome variable for the probit class of model
#' 
#' @param object object of class "probitcd"
#' 
#' @return Y vector of outcome variable
#' @export

out_dep.probitcd<- function(object){
    probs= stats::pnorm(object$lc)
    N = length(probs); Y = rep(0,N)
    set.seed(object$eta) #seed is passed as eta
    for (i in 1:N) {
      Y[i]=sample(c(1,0),1,replace = TRUE,prob = c(probs[i],(1-probs[i])))
    }
    Y
} 

#=============================================================================================>
#' Implemention - poisson model
#' 
#' Generate outcome variable for the poisson class of model
#' 
#' @param object object of class "poiscd"
#' 
#' @return Y vector of outcome variable
#' @export

   
  out_dep.poiscd<- function(object){
    mu= exp(object$lc)
    N = length(mu); Y = rep(0,N)
    set.seed(object$eta) #seed is passed as eta
    for (i in 1:N) {
      Y[i]=stats::rpois(1,lambda = mu[i])
    }
    Y
  } 

  #=============================================================================================>
  #' Implemention - negative binomial model
  #' 
  #' Generate outcome variable for the negative binomial class of model
  #' 
  #' @param object object of class "nbcd"
  #' 
  #' @return Y vector of outcome variable
  #' @export
  
  out_dep.nbcd<- function(object){
    mu= exp(object$lc)
    N = length(mu); Y = rep(0,N)
    theta = object$eta[1] #eta must be a vector of two elements: theta and seed
    seed = object$eta[2]; 
    set.seed(seed) #seed is passed as eta
    for (i in 1:N) {
      Y[i]=MASS::rnegbin(1,mu = mu[i],theta = theta)
    }
    Y
  } 
  
#=============================================================================================>
#' Model criterion function
#' 
#' A generic S3 function the Schwarz Information Criterion
#' 
#' @param object the object to be passed to the concrete class constructor \code{BIC}
#' 
#' @export
  
BIC<- function(object) UseMethod("BIC")

#=============================================================================================>
#' Schwarz Information Criterion
#' 
#' Compute the Schwarz Information Criterion for the general class of models
#' 
#' @param object object returned as output from \link{reg_cd}
#' 
#' @return BIC value
#' @export

BIC.regcd<- function(object){
  k = length(object$coefs)
  2*object$fval + (log(object$df+k)*k)
}
#=============================================================================================>
#' Schwarz Information Criterion - lmcd class
#' 
#' A class specific implementation of the Schwarz Information Criterion for "lmcd" class
#' 
#' @param object object returned as output from \link{reg_cd}
#' 
#' @return BIC value
#' @export

BIC.lmcd<- function(object){
  k = length(object$coefs)
  n = object$df+k
  n*log(object$fval/n) + (log(n)*k)
}

#=============================================================================================>
#' Model criterion function
#' 
#' A generic S3 function the Schwarz Information Criterion
#' 
#' @param object the object to be passed to the concrete class constructor \code{BIC}
#' 
#' @export

AIC<- function(object) UseMethod("AIC")

#=============================================================================================>
#' Schwarz Information Criterion
#' 
#' Compute the Schwarz Information Criterion for the general class of models
#' 
#' @param object object returned as output from \link{reg_cd}
#' 
#' @return BIC value
#' @export

AIC.regcd<- function(object){
  k = length(object$coefs)
  2*object$fval + 2*k
}
#=============================================================================================>
#' Schwarz Information Criterion - lmcd class
#' 
#' A class specific implementation of the Schwarz Information Criterion for "lmcd" class
#' 
#' @param object object returned as output from \link{reg_cd}
#' 
#' @return BIC value
#' @export

AIC.lmcd<- function(object){
  k = length(object$coefs)
  n = object$df+k
  n*log(object$fval/n) + 2*k
}
#=============================================================================================>
#' Model criterion function
#' 
#' A generic S3 function the Schwarz Information Criterion
#' 
#' @param object the object to be passed to the concrete class constructor \code{BIC}
#' 
#' @export

HQIC<- function(object) UseMethod("BIC")

#=============================================================================================>
#' Schwarz Information Criterion
#' 
#' Compute the Schwarz Information Criterion for the general class of models
#' 
#' @param object object returned as output from \link{reg_cd}
#' 
#' @return BIC value
#' @export

HQIC.regcd<- function(object){
  k = length(object$coefs)
  2*object$fval + 2*k*log(log(object$df+k))
}
#=============================================================================================>
#' Schwarz Information Criterion - lmcd class
#' 
#' A class specific implementation of the Schwarz Information Criterion for "lmcd" class
#' 
#' @param object object returned as output from \link{reg_cd}
#' 
#' @return BIC value
#' @export

HQIC.lmcd<- function(object){
  k = length(object$coefs)
  n = object$df+k
  n*log(object$fval/n) + 2*k*log(log(n))
}
#=============================================================================================>
#' Model criterion function
#' 
#' A generic S3 function as wrapper for internal R routines for classes of models implemented
#' in this package. See details \link{c_crit} for the list of classes supported.
#' 
#' @param object the object to be passed to the concrete class constructor \code{regir}
#' @param ... additional paramters to be passed to the internal routine
#' 
#' @export

regir<- function(object,...) UseMethod("regir")

#=============================================================================================>
#' Regression - lmcd class
#' 
#' A linear regression implementation for the "lmcd" class. It uses \code{\link[stats]{lm}} ##\code{lm()} from the stats package.
#' 
#' @param object a list of Y - outcome variable and Xmat - design matrix of class "lmcd"
#' @param ... additional parameters to be passed to \code{\link[stats]{lm}}
#' 
#' @return fitted model object
#' @export

regir.lmcd<- function(object,...){
  dat = data.frame(object$Y,object$Xmat); names(dat)[1]="Y"
  stats::lm(object$Y~.,data = dat,...)
}

#=============================================================================================>
#' Regression - logitcd class
#' 
#' A logit regression implementation for the "logitcd" class. It uses \code{\link[stats]{glm}}
#' 
#' @param object a list of Y - outcome variable and Xmat - design matrix of class "logitcd"
#' @param ... additional parameters to be passed to \code{\link[stats]{glm}}
#' 
#' @return fitted model object
#' @export

regir.logitcd<- function(object,...){
  fam = stats::binomial(link = "logit")
  dat = data.frame(object$Y,object$Xmat); names(dat)[1]="Y"
  stats::glm(object$Y~.,family = fam, data = dat,...)
}

#=============================================================================================>
#' Regression - poiscd class
#' 
#' A poisson regression implementation for the "poiscd" class. It uses \code{\link[stats]{glm}}
#' 
#' @param object a list of Y - outcome variable and Xmat - design matrix of class "poiscd"
#' @param ... additional parameters to be passed to \code{\link[stats]{glm}}
#' 
#' @return fitted model object
#' @export

regir.poiscd<- function(object,...){
  fam = stats::poisson(link = "log")
  dat = data.frame(object$Y,object$Xmat); names(dat)[1]="Y"
  stats::glm(Y~.,family = fam,data = dat,...)
}

#=============================================================================================>
#' Regression - probitcd class
#' 
#' A probit regression implementation for the "probitcd" class. It uses \code{\link[stats]{glm}}
#' 
#' @param object a list of Y - outcome variable and Xmat - design matrix of class "probitcd"
#' @param ... additional parameters to be passed to \code{\link[stats]{glm}}
#' 
#' @return fitted model object
#' @export

regir.probitcd<- function(object,...){
  fam = stats::binomial(link = "probit")
  dat = data.frame(object$Y,object$Xmat); names(dat)[1]="Y"
  stats::glm(object$Y~.,family = fam, data = dat,...)
}

#=============================================================================================>
#' Regression - nbcd class
#' 
#' A negative binomial regression implementation for the "nbcd" class. It uses \code{\link[MASS]{glm.nb}}
#' 
#' @param object a list of Y - outcome variable and Xmat - design matrix of class "nbcd"
#' @param ... additional parameters to be passed to \code{\link[MASS]{glm.nb}}
#' 
#' @return fitted model object
#' @export

regir.nbcd<- function(object,...){
  dat = data.frame(object$Y,object$Xmat); names(dat)[1]="Y"
  MASS::glm.nb(object$Y~.,data = dat,...,maxit = 75)
}

#=============================================================================================>
#' Regression - qreg class
#' 
#' A quantile regression implementation for the "qreg" class. It uses \code{\link[quantreg]{rq}}
#' 
#' @param object a list of Y - outcome variable and Xmat - design matrix of class "qreg"
#' @param ... additional parameters to be passed to \code{\link[quantreg]{rq}}
#' 
#' @return fitted model object
#' @export

regir.qreg<- function(object,...){
  dat = data.frame(object$Y,object$Xmat); names(dat)[1]="Y"
  quantreg::rq(object$Y~.,data = dat,...)
}
