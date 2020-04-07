genpar <- function(X, Y, rho, scale=TRUE, delta.EBIC=0.5, delta.scalefree = 1, maxit.weight=100, maxit.glasso = 10, n.test = 10, nfolds.glmnet=10, msg = TRUE){
    if(class(X)!="matrix") stop("\"X\" must be a matrix.")
    if(class(Y)!="numeric") stop("\"Y\" must be a numeric vector.")
    if(nrow(X)!=length(Y)) stop("The number of rows in \"X\" must be equal to length of \"Y\".")
    if(class(rho)!="numeric") stop("\"rho\" must be a numeric vector.")
    if(sum(rho>0)!=length(rho>0)) stop("All of the elements of \"rho\" must be positive values.")
    if(class(scale)!="logical") stop("\"scale\" must be logical.")
    if(class(delta.EBIC)!="numeric") stop("\"delta.EBIC\" must be numeric.")
    if(delta.EBIC<=0) stop("\"delta.EBIC\" must be positive")
    if(class(delta.scalefree)!="numeric") stop("\"delta.scalefree\" must be numeric.")
    if(delta.scalefree<=0) stop("\"delta.scalefree\" must be positive")
    if(class(maxit.weight)!="numeric") stop("\"maxit.weight\" must be numeric.")
    if(maxit.weight<=0) stop("\"maxit.weight\" must be a positive integer.")
    if(class(maxit.glasso)!="numeric") stop("\"maxit.glasso\" must be numeric.")
    if(maxit.glasso<=0) stop("\"maxit.glasso\" must be a positive integer.")
    if(class(n.test)!="numeric") stop("\"n.test\" must be numeric.")
    if(n.test<=0) stop("\"n.test\" must be a positive integer.")
    if(class(nfolds.glmnet)!="numeric") stop("\"nfolds.glmnet\" must be numeric.")
    if(nfolds.glmnet<=0) stop("\"nfolds.glmnet\" must be a positive integer.")
    if(class(msg)!="logical") stop("\"msg\" must be logical.")


    n <- nrow(X)
    p <- ncol(X)
    if(scale) X <- scale(X)
    covx <- t(X)%*%X/(n-1)
    
    nrho <- length(rho)
    BIC <- EBIC <- niter <- rep(NA, nrho)
    convergence <- rep(TRUE, nrho)
    omega2.list <- sigma2.list <- vector(nrho, mode="list")
    omega3.list <- sigma3.list <- vector(nrho, mode="list")
    
    for(irho in 1:nrho){
        rho0 <- rho[irho]
        tol <- 1
        count <- 1
        tolall <- NULL
        condition <- TRUE
        wmat <- matrix(1,p,p)
        diag(wmat) <- 0
        
        if(msg) cat(paste("network estimation for rho = ", rho[irho], ": computing", sep = ""))
        while(condition){
            #fit glasso
            fit <- glasso(covx, rho=rho0*wmat, maxit=maxit.glasso) 
            
            #calculate weights
            if(count > 1) omegaprevious <- omega
            wmat <- matrix(1,p,p)
            diag(wmat) <- 0
            omega <- Matrix(fit$wi)
            sigma <- Matrix(fit$w)
            w <- 1/(Matrix::colSums(abs(omega)) - abs(Matrix::diag(omega)) + delta.scalefree)
            wmat <- sqrt(tcrossprod(w))
            diag(wmat) <- 0
            
            #calculate model selection criteria
            if(count==1){
                loglik = -n/2*(-log(Matrix::det(omega)) + sum(Matrix::diag(omega%*%covx)))
                BIC[irho] <- -2*(loglik) + (sum(omega!=0)-p)/2*log(n)
                EBIC[irho] <- -2*(loglik) + (sum(omega!=0)-p)/2*log(n) + 4 * (sum(omega!=0)+p)/2 * delta.EBIC * log(p)
            }
            
            if(count>1){
                tol <- sum(abs(omegaprevious-omega))
            }else{
                tol <- 1
            }
            tolall[count] <- tol
            
            condition <- tol > 1e-5 && count <= maxit.weight
          
            if(count == 1){
                omega2.list[[irho]] <- omega
                sigma2.list[[irho]] <- sigma
            }            
            if(!condition){
                omega3.list[[irho]] <- omega
                sigma3.list[[irho]] <- sigma
                if(count == maxit.weight) convergence[irho] <- FALSE
                niter[irho] <- count
            }
            count <- count+1
            if (msg) cat(".")
        }
        if(msg) cat(paste("\n", sep = ""))
    }
    
    
    EBICbest <- which.min(EBIC)
    Sigma.all <- list(Sigma1 = covx, Sigma2 = sigma2.list[[EBICbest]], Sigma3 = sigma3.list[[EBICbest]])
    
    # determine true beta and variance of error
    num.test <- sample(1:n, n.test)
    X.train <- X[-num.test, ]
    X.test <- X[num.test, ]
    Y.train <- Y[-num.test]
    Y.test <- Y[num.test]
    
    res <- cv.glmnet(X.train, Y.train, nfolds=nfolds.glmnet)
    Y.pred <- predict(res, newx = X.test, s = "lambda.min")
    #sd.error <- sqrt(sum((Y.pred-Y.test)^2)/10)
    sd.error <- sqrt(sum((Y.pred-Y.test)^2)/num.test)
    
    # prepare for generating data (determine mu and A)
    mu <- colMeans(X)
    #mu <- apply(X, 2, mean)
    
    # if Sigmas are not full rank, we add a tiny number to diagonal element
    A.all <- try(lapply(Sigma.all, chol), silent = TRUE)
    if(class(A.all) == "try-error"){
        Sigma.all <- lapply(Sigma.all, function(x){x <- x + x[1, 1] * 1e-5 * diag(nrow(x))})
        A.all <- lapply(Sigma.all, chol)
    }

    ans <- list(Sigma=Sigma.all, A=A.all, res=res, beta=coef(res), sd.error=sd.error, omega2=omega2.list, omega3=omega3.list, BIC=BIC, EBIC=EBIC, rho=rho, p=p, mu=mu)
    ans$call <- match.call()
    class(ans) <- "genpar"
    ans
}