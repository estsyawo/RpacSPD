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

#==================================================================================================>
#' Concrete class constructor
#' 
#' A function for constructing functions for concrete classes of models for the \code{crit} family of
#' of functions.
#' 
#' @param Y vector of the outcome variable
#' @param lc linear combination \code{lc = X'beta}
#' @param eta for models needing extra parameter besides beta in the linear combination \code{lc}. 
#' eta can be a vector
#' 
#' @param klass the class of model. Currently "lm" for linear regression, "logit" (logit model), 
#' "poisson" (poisson model), "probit" (probit model), "negbin" (negative binomial) are supported.
#' 
#' @return object an object list with class attribute klass.
#' 
#' @export

c_crit<- function(Y,lc,eta=NULL,klass="lm"){ #assign class of object
  if(is.null(eta)){ 
    object<- list(Y,lc)
    class(object)<- klass
    names(object)<- c("Y","lc")  
  }else{ 
    object<- list(Y,lc,eta)
    class(object)<- klass
    names(object)<- c("Y","lc","eta")  
  }
  object
}

#==================================================================================================>
#' Implemention - linear regression
#' 
#' This function evaluates the criterion for the linear regression class of model. 
#' 
#' @param object Object with class attribute "lm" obtained from the function \code{c_crit()}
#' @return sum of squared errors
#' 
#' @export

crit.lm<- function(object){sum((object$Y-object$lc)^2)} #for least squares

#==================================================================================================>
#' Implemention - logit model
#' 
#' This function evaluates the criterion for the logit class of model. 
#' 
#' @param object Object with class attribute "logit" obtained from the function \code{c_crit()}
#' @return negative log likelihood
#' 
#' @export

crit.logit<- function(object) {-sum(object$Y*object$lc - log(1+exp(object$lc))) } #return neg. log-like for logit

#==================================================================================================>
#' Implemention - poisson model
#' 
#' This function evaluates the criterion for the poisson class of model. 
#' 
#' @param object Object with class attribute "pois" obtained from the function \code{c_crit()}
#' @return negative log likelihood
#' 
#' @export

crit.pois<- function(object) {-sum(object$Y*object$lc - exp(object$lc))} #return neg. log-like for poisson

#==================================================================================================>
#' Implemention - probit model
#' 
#' This function evaluates the criterion for the probit class of model. 
#' 
#' @param object Object with class attribute "probit" obtained from the function \code{c_crit()}
#' @return negative log likelihood
#' 
#' @export

crit.probit<- function(object){
  plc = pnorm(object$lc)
  lk = -sum(object$Y*log(plc)+(1-object$Y)*log(1-plc))
  lk
} 

#==================================================================================================>
#' Implemention - negative binomial model
#' 
#' This function evaluates the criterion for the negative binomial class of model. 
#' 
#' @param object Object with class attribute "negbin" obtained from the function \code{c_crit()}
#' @return negative log likelihood
#' 
#' @export

crit.negbin<- function(object){
  if(is.null(object$eta)){
    stop("Error: This negative binomial class requires an additional parameter eta")
  } #end if
  eta2=1/(object$eta^2)
  mu = exp(object$lc); dn = eta2+mu
  lv = eta2*log(eta2/dn) + object$Y*log(mu/dn) + log(gamma(object$Y+eta2)/gamma(eta2))
  -sum(lv)
} #return the neg. log-like for negative binomial

#==================================================================================================>
#' Implemention - linear regression model
#' 
#' This function evaluates to a vector of squared residuals for the linear regression class of model. 
#' 
#' @param object Object with class attribute "lm" obtained from the function \code{c_crit()}
#' @return vector of squared residuals
#' 
#' @export

crit_obs.lm<- function(object){c((object$Y-object$lc)^2)} #for least squares

#==================================================================================================>
#' Implemention - logit model
#' 
#' This function evaluates to a vector of log-likelihoods for the logit class of model. 
#' 
#' @param object Object with class attribute "logit" obtained from the function \code{c_crit()}
#' @return vector of negative log likelihoods
#' 
#' @export

crit_obs.logit<- function(object) {-c(object$Y*object$lc - log(1+exp(object$lc))) } #return neg. log-like for logit

#==================================================================================================>
#' Implemention - poisson model
#' 
#' This function evaluates to a vector of log-likelihoods for the poisson class of model. 
#' 
#' @param object Object with class attribute "pois" obtained from the function \code{c_crit()}
#' @return vector of negative log likelihoods
#' 
#' @export

crit_obs.pois<- function(object) {-c(object$Y*object$lc - exp(object$lc))} #return neg. log-like for poisson

#==================================================================================================>
#' Implemention - probit model
#' 
#' This function evaluates to a vector of log-likelihoods for the probit class of model. 
#' 
#' @param object Object with class attribute "probit" obtained from the function \code{c_crit()}
#' @return vector of negative log likelihoods
#' 
#' @export

crit_obs.probit<- function(object){
  plc = pnorm(object$lc)
  lk = -c(object$Y*log(plc)+(1-object$Y)*log(1-plc))
  lk
}

#==================================================================================================>
#' Implemention - negative binomial model
#' 
#' This function evaluates to a vector of log-likelihoods for the negative binomial class of model. 
#' 
#' @param object Object with class attribute "negbin" obtained from the function \code{c_crit()}
#' @return vector of negative log likelihoods
#' 
#' @export

crit_obs.negbin<- function(object){
  if(is.null(object$eta)){
    stop("Error: This negative binomial class requires an additional parameter eta")
  } #end if
  eta2=1/(object$eta^2)
  mu = exp(object$lc); dn = eta2+mu
  lv = eta2*log(eta2/dn) + object$Y*log(mu/dn) + log(gamma(object$Y+eta2)/gamma(eta2))
  -c(lv)
}


#=============================================================================================>
#' Randomly generate outcome Y
#' 
#' A generic function for randomly generating outcome variables for model classes.
#' 
#'  @param object Object of class attribute obtained from \code{c_out_dep()}
#'  
#'  @export

out_dep<- function(object) UseMethod("out_dep")

#=============================================================================================>
#' Object class constructor
#' 
#' This function constructs an object class for the generic function \code{out_dep()}
#' 
#' @param Y vector of the outcome variable
#' @param lc linear combination \code{lc = X'beta}
#' @param eta for models needing extra parameter besides beta in the linear combination \code{lc}. 
#' eta can be a vector. For "lm", \code{eta} is the vector of errors; must of the same length as Y
#' For "negbin", \code{eta} must be a vector of two elements; the first one is parameter \code{th}
#' and the second is the seed for random generation. For other classes, eta is the seed for random 
#' generation.
#' 
#' @param klass the class of model. Currently "lm" for linear regression, "logit" (logit model), 
#' "poisson" (poisson model), "probit" (probit model), "negbin" (negative binomial) are supported.
#' 
#' @return object an object list with class attribute klass.
#' 
#' @export

c_out_dep<- function(lc,eta=NULL,klass="lm"){ #assign class object
  if(is.null(eta)){ 
    object<- lc
    class(object)<- klass
    names(object)<- c("lc")  
  }else{ 
    object<- list(lc,eta)
    class(object)<- klass
    names(object)<- c("lc","eta")  
  }
  object
}

#=============================================================================================>
#' Implemention - linear regression model
#' 
#' Generate outcome variable for linear regression class of model
#' 
#' @param object object of class "lm"
#' 
#' @return Y vector of outcome variable

out_dep.lm<- function(object){
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
#' @param object object of class "logit"
#' 
#' @return Y vector of outcome variable
#' @export

out_dep.logit<- function(object){
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
#' @param object object of class "probit"
#' 
#' @return Y vector of outcome variable
#' @export

out_dep.probit<- function(object){
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
#' @param object object of class "pois"
#' 
#' @return Y vector of outcome variable
#' @export

   
  out_dep.pois<- function(object){
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
  #' @param object object of class "negbin"
  #' 
  #' @return Y vector of outcome variable
  #' @export
  
  out_dep.negbin<- function(object){
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

