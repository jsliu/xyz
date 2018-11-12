unscale<-function(X) {
  mu_X<-attr(X,"scaled:center");sigma_X<-attr(X,"scaled:scale")
  return(t(apply(X, 1, function(r)r*mu_X + attr(X, 'scaled:center'))))
}

standardize_result<-function(result,X,Y,standardize,standardize_response) {
  L<-length(result[[1]])
  main_effects<-result[[1]];
  beta_main<-result[[2]];
  intr_effects<-result[[3]];
  beta_intr<-result[[4]];
  intercept<-result[[6]];

  if(standardize) {
    mu_X <- colMeans(X); sigma_X <- apply(X,2,sd)

    for(l in 1:L) {
      for(i in 1:length(main_effects[[l]])) {
        beta_main[[l]][i]<-beta_main[[l]][i]/sigma_X[main_effects[[l]][i]]
        intercept[l]<-intercept[l]-beta_main[[l]][i]*mu_X[main_effects[[l]][i]]
      }
      for(i in 1:length(beta_intr[[l]])) {
        variable<-X[,intr_effects[[l]][1,i]]*X[,intr_effects[[l]][2,i]]
        mu_XX<-mean(variable)
        sigma_XX<-stats::sd(variable)
        beta_intr[[l]][i]<-beta_intr[[l]][i]/sigma_XX
        intercept[l]<-intercept[l]-beta_intr[[l]][i]*mu_XX
      }
    }
  }

  if(standardize_response) {
    mu_Y<-attr(Y,"scaled:center"); sigma_Y<-attr(Y,"scaled:scale")
    for(l in 1:L) {
      beta_main[[l]]<-beta_main[[l]]*sigma_Y
      beta_intr[[l]]<-beta_intr[[l]]*sigma_Y
      intercept[l]<-sigma_Y*intercept[l]+mu_Y
    }
  }

  result[[2]] <- beta_main
  result[[4]] <- beta_intr
  result[[6]] <- intercept

  return(result)
}

#' Elasticnet with interactions (glmnet)
#' @param X A matrix.
#' @param Y A vector.
#' @param L An integer indicating how many projection steps are performed.
#' @param standardize A boolean indicating if X should be scaled and centered.
#' @param standardize_response A boolean indicating if Y should be scaled and centered.
#' @param n_lambda A natural number indicating how long the path of lambdas should be.
#' @param lambdas A vector of decreasing real numbers containing user specified values of lambda.
#' @param alpha  A real number between 0 and 1 (the elastic net parameter)
#' @return \code{N} strongest interactions (of type \code{type}) between \code{X} and \code{Y} after \code{L} projections.
#' @references G. Thanei, N. Meinshausen and R. Shah (2016). The xyz algorithm for fast interaction search in high-dimensional data. <https://arxiv.org/pdf/1610.05108v1.pdf>
#' @examples
#' n<-300
#' p<-1000
#' #build matrix of predictors
#' X<-matrix(rnorm(n*p),n,p)
#' #build a main effect and an interaction into Y
#' Y<-4*X[,1]*X[,2]-5*X[,4]+rnorm(n)
#' result<-xyz_regression(X,Y,n_lambda=10,alpha=0.9,L=10)
#' #print the result
#' print(result)
#' #plot the result
#' plot(result)
#' @export
#' @useDynLib xyz2
#' @importFrom Rcpp sourceCpp
xyz_regression<-function(X,Y,weights=NULL,lambdas=NULL,n_lambda=10,alpha=0.9,L=10,standardize_response=FALSE) {
  L<-round(L)
  if(L < 1) {
    stop("Number of runs has to be at least 1.")
  }
  if(L > 100) {
    warning("You choose very high number of runs.")
  }
  if(n_lambda < 1) {
    stop("Length of lambda sequence cannot be smaller than 1.")
  }

  #do translating checks here
  if(!is.matrix(X)) {
    X<-as.matrix(X)
    if(!is.matrix(X)) {
      stop("X cannot be coerced to a matrix.")
    }
  }

  n<-dim(X)[1]
  p<-dim(X)[2]
  if(is.null(weights)) {
    weights <- rep(1.0/n,n)
  }
  if (alpha < 0 || alpha > 1)
    stop("alpha has to be between 0 and 1.")

  if(!is.vector(Y)) {
    stop("Y has to be a vector.")
  }
  if(length(Y) != n) {
    stop("Y and X have to have the same number of rows")
  }
  if(length(weights) != n){
    stop("weights and Y have to have the same number of rows")
  }
  if(n < 10) {
    stop(paste("You have ",n," samples. The number of samples should at least be 10.",sep=""))
  }

  if(standardize_response) {
    Y<-scale(Y)
  }
  if(is.null(lambdas)) {
    lambdas<-rep(-1,n_lambda)
  }

  wts <- weights/sum(weights)

  max_main_effects<-100
  max_interaction_effects<-20
  standardize <- TRUE
  result<-gaussiglmnet(X, Y, wts, lambdas, alpha, standardize, max_main_effects, max_interaction_effects, 2, L)
  L<-length(result[[1]])
  for(i in 1:L) {
    result[[1]][[i]]<-result[[1]][[i]]+1
    result[[3]][[i]]<-result[[3]][[i]]+1
    result$df[i] <-length(result[[2]][[i]])+length(result[[4]][[i]])
  }

  result<-standardize_result(result,X,Y,standardize,standardize_response)
  class(result)<-"xyz_regression_result"
  return(result)
  stop("You reached the end of the function and it won't return anything. This is not good and u should be ashamed.")
}

#' @export
summary.xyz_regression_result<-function(object,...) {
  x<-object
  l_main_1<-length(x[[1]])
  l_intr<-dim(x[[3]][[l_main_1]])[2]
  l_main<-length(x[[1]][[l_main_1]])
  if(l_main==1 & x[[2]][[l_main_1]][1]==0) {
    l_main = 0
  }
  if(l_intr==1 & x[[4]][[l_main_1]][1]==0) {
    l_intr = 0
  }
  output<-paste("Discovered main effects: ",l_main," Discovered interaction effects: ",l_intr,sep="")
  return(output)
}

#' @export
print.xyz_regression_result<-function(x,whichlambda=-1,...) {
  lambda_sequence<-round(x[[5]],digits=5)
  cat("Lambda sequence:\n")
  cat(paste("lambda",1:length(lambda_sequence),"=",lambda_sequence,"\n",sep=""))
  l_main_1<-length(x[[1]])
  l_intr<-dim(x[[3]][[l_main_1]])[2]
  l_main<-length(x[[1]][[l_main_1]])
  if(l_main==1 & x[[2]][[l_main_1]][1]==0) {
    l_main = 0

  }
  if(l_intr==1 & x[[4]][[l_main_1]][1]==0) {
    l_intr = 0
  }
  cat(paste("Discovered main effects: ",l_main," Discovered interaction effects: ",l_intr,"\n",sep=""))
  l<-whichlambda
  cat("Model parameters:\n")
  cat("intercept: ",x[[6]][1],"\n",sep="")
  if(l_main + l_intr > 0) {
    if(l<1) {
      l<-length(lambda_sequence)
    } else {
      l <- round(l)
      if(l < 1 | l > length(lambda_sequence)) {
        l<-length(lambda_sequence)
      }
    }
    cat(paste("Printing effects for lambda",l,"=",lambda_sequence[l],"\n",sep=""))
  }
  if(l_main > 0) {
    cat("Main effects:\n")
    ord<-order(abs(x[[2]][[l]]),decreasing=TRUE)
    x[[1]][[l]]<-x[[1]][[l]][ord]
    x[[2]][[l]]<-x[[2]][[l]][ord]
    for(i in 1:length(x[[1]][[l]])) {
      cat("Main effect: ",x[[1]][[l]][i]," coefficient: ",x[[2]][[l]][i],"\n",sep="")
    }
  }
  if(l_intr > 0) {
    cat("Interaction effect:\n")
    ord<-order(abs(x[[4]][[l]]),decreasing=TRUE)
    x[[3]][[l]]<-matrix(x[[3]][[l]][,ord],2,length(x[[4]][[l]]))
    x[[4]][[l]]<-x[[4]][[l]][ord]
    for(i in 1:length(x[[4]][[l]])) {
      cat("Interaction effect: (",x[[3]][[l]][1,i],",",x[[3]][[l]][2,i],") coefficient: ",x[[4]][[l]][i],"\n",sep = "")
    }
  }
}

#' @importFrom graphics matplot legend
#' @export
plot.xyz_regression_result<-function(x,xvar=c("norm","lambda")) {
  if( length(x) != 7 ) {
    stop("fit has parts missing")
  }
  main_effects<-x[[1]]
  main_coefs<-x[[2]]
  intr_effects<-x[[3]]
  intr_coefs<-x[[4]]
  lambdas<-x[[5]]
  df <- x[[7]]
  path_length<-length(lambdas)

  l<-length(x[[1]])

  #find all main and intr effects
  all_main <- unique(unlist(main_effects))
  all_intr <- matrix(unlist(intr_effects),ncol=2,byrow=TRUE)
  all_intr <- all_intr[!duplicated(all_intr),]

  nr_main<-length(all_main)
  nr_intr<-length(all_intr)/2

  main_path<-list()
  intr_path<-list()
  for(i in 1:nr_main) {
    temp<-rep(0,path_length)
    for(j in 1:path_length) {
      temp_pos<- -1
      for(k in 1:length(main_effects[[j]])) {
        if(main_effects[[j]][k] == all_main[i]) {
          temp_pos<-k
        }
      }
      if(temp_pos>0) {
        temp[j] <- main_coefs[[j]][temp_pos]
      } else {
        temp[j] <-0
      }
    }
    main_path[[i]]<-temp
  }
  if(nr_intr>1) {
    for(i in 1:nr_intr) {
      temp<-rep(0,path_length)
      for(j in 1:path_length) {
        temp_pos<- -1
        intr_effects[[j]]<-matrix(intr_effects[[j]],2,length(intr_coefs[[j]]))
        for(k in 1:dim(intr_effects[[j]])[2]) {

          if(intr_effects[[j]][1,k] == all_intr[i,1] && intr_effects[[j]][2,k] == all_intr[i,2]) {
            temp_pos<-k
          }
          if(intr_effects[[j]][2,k] == all_intr[i,1] && intr_effects[[j]][1,k] == all_intr[i,2]) {
            temp_pos<-k
          }
        }
        if(temp_pos>0) {
          temp[j] <- intr_coefs[[j]][temp_pos]
        } else {
          temp[j] <-0
        }
      }
      intr_path[[i]]<-temp
    }
  } else {
    intr_path[[1]]<- unlist(intr_coefs)
  }
  all_paths<-matrix(c(unlist(main_path),unlist(intr_path)),nrow=nr_main+nr_intr,byrow=TRUE)

  beta_norm<-colSums(abs(all_paths))

  xvar <- match.arg(xvar)
  switch(xvar,
        "norm"={xaxis=beta_norm; xlab="L1 Norm"; approx.f=1},
        "lambda"={xaxis=log(lambdas); xlab="Log Lambda"; approx.f=0}
        )
  yaxis <- t(all_paths)
  colors_plot<-c(rep("blue",nr_main),rep("red",nr_intr))
  matplot(xaxis,yaxis,type="l",col=colors_plot,lty=1,lwd=2,ylab="coefficients",xlab=xlab)
  atdf <- pretty(xaxis)
  prettydf <- approx(x=xaxis,y=df,xout=atdf,rule=2,method="constant",f=approx.f)$y
  axis(3,at=atdf,labels=prettydf,tcl=NA)
  legend("topleft",legend=c("main effects","interaction effects"),lty=c(1,1),col=c("blue","red"))
}

#' @export
predict.xyz_regression_result<-function(object,newdata,...) {
  l<-length(object[[1]])
  Y_pred_matrix<-matrix(0,dim(newdata)[1],l)
  for(i in 1:l) {
    main_effects<-c(object[[1]][[i]],1)
    beta_main<-c(object[[2]][[i]],0)
    intr_effects<-cbind(object[[3]][[i]],c(1,1))
    beta_intr<-c(object[[4]][[i]],0)
    intercept<-object[[6]][i]
    Y_pred_matrix[,i]<-intercept+newdata[,main_effects]%*%beta_main+(newdata[,intr_effects[1,]]*newdata[,intr_effects[2,]])%*%beta_intr
  }
  return(Y_pred_matrix)
}

#' @export
ic.xyz_regression<-function(x,y,crit=c("bic","aic","aicc","hqc"),...)
{
  if (is.matrix(x) == FALSE) {
    x <- as.matrix(x)
  }
  if (is.vector(y) == FALSE) {
    y <- as.vector(y)
  }
  crit<-match.arg(crit)
  n<-length(y)
  model <- xyz_regression(X = x, Y = y, ...)
  lambda <- model[[5]]
  df <- sapply(model[[1]],length) + sapply(model[[3]],length)

  yhat <- predict(model,x,...)
  residuals <- (y - yhat)
  mse <- colMeans(residuals^2)
  sse <- colSums(residuals^2)
  nvar <- df + 1
  bic <- n*log(mse)+nvar*log(n)
  aic <- n*log(mse)+2*nvar
  aicc <- aic+(2*nvar*(nvar+1))/(n-nvar-1)
  hqc <- n*log(mse)+2*nvar*log(log(n))
  sst <- (n-1)*var(y)
  r2 <- 1 - (sse/sst)
  adjr2 <- (1 - (1 - r2) * (n - 1)/(nrow(x) - nvar - 1))

  crit <- switch(crit,bic=bic,aic=aic,aicc=aicc,hqc=hqc)

  selected <- best.model <- max(which(crit == min(crit)))

  ic <- c(bic=bic[selected],aic=aic[selected],aicc=aicc[selected],hqc=hqc[selected])

  result <- list(best.model=selected,ic=ic,lambda=lambda[[selected]], nvar=nvar[selected],
              glmnet=model,residuals=residuals[,selected],fitted.values=yhat[,selected],
              ic.range=crit, call = match.call())

  class(result) <- "ic.xyz_regression"
  return(result)
}

#' @export
coef.ic.xyz_regression <- function(object)
{
  selected <- object$best.model
  model <- object$glmnet
  main_effects <- c(model[[1]][[selected]],1)
  beta_main <- c(model[[2]][[selected]],0)
  intr_effects <- cbind(model[[3]][[selected]],c(1,1))
  beta_intr <- c(model[[4]][[selected]],0)
  intercept <- model[[6]][selected]

  return(list(main_effects=main_effects,beta_main=beta_main,intr_effects=intr_effects,beta_intr=beta_intr,intercept=intercept))
}

#' @export
predict.ic.xyz_regression <- function(object, newdata, ...)
{
  coeffs <- coef(object)
  main_effects <- coeffs$main_effects
  beta_main <- coeffs$beta_main
  intr_effects <- coeffs$intr_effects
  beta_intr <- coeffs$beta_intr
  intercept <- coeffs$intercept
  Y_pred <- intercept+newdata[,main_effects]%*%beta_main+(newdata[,intr_effects[1,]]*newdata[,intr_effects[2,]])%*%beta_intr

  return(as.vector(Y_pred))
}

#' @export
plot.ic.xyz_regression <- function(object,...)
{
  n <- object$glmnet$df
  llambda=log(object$glmnet[[5]])
  ic=object$ic.range
  ylab=names(object$ic[which(object$ic==ic[which.min(ic)])])

  plot(llambda, ic, xlab = "log(Lambda)", ylab = substr(ylab, 1, 3),...)
  abline(v = llambda[which.min(ic)], lty = 2)
  axis(3, at = llambda, labels = n)
}

#' #' @export
#' covariates <- function(x,y,wei,standardize)
#' {
#'   cov.mat <- absolute_covariates(x,y,wei,standardize)
#'   return(cov.mat)
#' }
#'
#' #' @export
#' createLambda <- function(max.cov,alpha,n_lambda,factor.inv)
#' {
#'   lambdas <- create_lambda_sequence(max.cov,alpha,n_lambda,factor.inv)
#'   return(lambdas)
#' }
