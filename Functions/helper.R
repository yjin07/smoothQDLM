library(Rcpp)
library(quantreg)
library(hqreg)
source("Functions/ADMM.R")            
source("Functions/utils.R")           
sourceCpp("Functions/utils.cpp") 


# ---------------------------------------------
# Fit quantile regression with Unimodal penalty
# ---------------------------------------------
# y : n * 1 response 
# X : list of length num.expo. 
# Z : n * p2 covariates
# lam1.vec : tuning parameters for Unimodal
# lam2.vec : tuning parameters for Smoothness
# ---------------------------------------------
cv.q.unimodal <- function(y, X, Z, yval, Xval, Zval, lam1.vec, lam2.vec, tau0, sigma=1e-4, verbose = 0) {
    n <- length(y); nval <- length(yval)
    p2 <- ncol(Z)
    num.expo <- length(X)
    p.expo <- rep(0, num.expo)
    for (i in 1 : num.expo) {
        p.expo[i] <- ncol(X[[i]])
    }
    p.aggre <- sum(p.expo)

    nlam1 <- length(lam1.vec); nlam2 <- length(lam2.vec)
    hat.beta.mat <- array(0, dim = c(p.aggre, nlam1 + 1, nlam2))
    hat.gamma.mat <- array(0, dim = c(p2, nlam1 + 1, nlam2))
    m.list <- array(0, dim = c(num.expo, nlam1, nlam2))
    cvErr <- matrix(0, nrow = nlam1, ncol = nlam2)

    # Position matrix
    pos.Mat <- matrix(nrow = num.expo, ncol = 2)
    start_pos <- 1
    for(i in 1:num.expo) {
        end_pos <- start_pos + p.expo[i] - 1
        pos.Mat[i, ] <- c(start_pos, end_pos)
        start_pos <- end_pos + 1
    }

    D.list <- list()
    for (jjj in 1 : num.expo) {
        D.list[[jjj]] <- Dmat(p.expo[jjj])
    }
    D <- diagonal_matrix_bind(D.list)

    Xcat <- NULL
    Xvalcat <- NULL
    for (jjj in 1 : num.expo) {
        Xcat <- cbind(Xcat, X[[jjj]])
        Xvalcat <- cbind(Xvalcat, Xval[[jjj]])
    }

    # ! ------------------------------------------
    # ! Store frequently used variables ahead
    # ! ------------------------------------------
    Z.tilde <- rbind(Z, matrix(0, nrow = nrow(D), ncol = p2))
    y.tilde <- c(y, rep(0, nrow(D)))
    XZcat <- cbind(Xcat, Z)
    ZtZZt <- solve(crossprod(Z.tilde), t(Z.tilde))
    # ! ------------------------------------------


    rqfit0 <- rq(y ~ Xcat + Z - 1, tau = tau0)
    hat.beta.qr <- rqfit0$coef[1 : ncol(Xcat)]
    hat.gamma.qr <- rqfit0$coef[(ncol(Xcat) + 1) : (ncol(Xcat) + p2)]

    for (iter2 in 1 : nlam2) {
        if (verbose > 0) {
            cat("\n\n========================================\n")
            cat(iter2, "th lam2:", lam2.vec[iter2], "\n")
            cat("========================================\n")
        }
        lam2 <- lam2.vec[iter2]

        hat.beta.mat[, 1, iter2] <- hat.beta.qr
        hat.gamma.mat[, 1, iter2] <- hat.gamma.qr

        # ! ------------------------------------------
        # ! Store frequently used variables ahead
        # ! ------------------------------------------
        Xcat.tilde <- rbind(Xcat, sqrt(2 * lam2 / sigma) * D)
        ZtZZtX <- ZtZZt %*% Xcat.tilde
        PzX.tilde <- Z.tilde %*% ZtZZtX
        Xcat.bar <- Xcat.tilde - PzX.tilde
        eta <- max(eigen(crossprod(Xcat.bar))$values)
        # ! ------------------------------------------

        for (i in 1 : nlam1) {
            lam1 <- lam1.vec[i]
            if (verbose > 1) {
                cat("\n+++++++++++++++++++++++++++++\n")
                cat(i, "th lam1:", lam1, "\n")
                cat("+++++++++++++++++++++++++++++\n")
            }

            hat.beta <- hat.beta.mat[, i, iter2]
            hat.gamma <- hat.gamma.mat[, i, iter2]
            hat.rr <- y - Xcat %*% hat.beta - Z %*% hat.gamma
            m.expo <- rep(0, num.expo)
            for (jjj in 1 : num.expo) {
                obj.mm <- rep(0, p.expo[jjj] + 1)
                for (mm in 0 : p.expo[jjj]) {
                    obj.mm[mm + 1] <- eval_unimodal(hat.beta[pos.Mat[jjj, 1] : pos.Mat[jjj, 2]], mm)
                }
                m.expo[jjj] <- which.min(obj.mm) - 1
            }
            if (verbose > 2) {
                cat("maximal loc: ", which.max(hat.beta), "\n")
                cat("m start:", m.expo, '\n')
            }
            obj.pre <- Inf
            obj <- eval_obj_uni(y, Xcat, Z, hat.beta, hat.gamma, tau0, m.expo, lam1, lam2, D, pos.Mat)
            num.iter <- 0
            while (abs(obj.pre - obj) > 1e-6 && num.iter < 1e3) {
                tmp <- runADMM.uni(y, Xcat, Z, 
                                y.tilde, Z.tilde, XZcat, ZtZZt, Xcat.tilde, Xcat.bar, eta,
                                hat.beta, hat.gamma, hat.rr, tau0, m.expo, sigma, lam1, lam2, D, pos.Mat)
                hat.beta <- tmp$beta
                hat.gamma <- tmp$gamma
                hat.rr <- tmp$rr
                m.expo.pre <- m.expo
                m.expo <- rep(0, num.expo)
                for (jjj in 1 : num.expo) {
                    obj.mm <- rep(0, p.expo[jjj] + 1)
                    for (mm in 0 : p.expo[jjj]) {
                        obj.mm[mm + 1] <- eval_unimodal(hat.beta[pos.Mat[jjj, 1] : pos.Mat[jjj, 2]], mm)
                    }
                    m.expo[jjj] <- which.min(obj.mm) - 1
                }
                num.iter <- num.iter + 1
                obj.pre <- obj
                obj <- eval_obj_uni(y, Xcat, Z, hat.beta, hat.gamma, tau0, m.expo, lam1, lam2, D, pos.Mat) 
                if (verbose > 3) {
                    obj.Show <- sprintf("%.7f", obj)
                    cat("-----> Iteration:", num.iter, "m: ", m.expo, " obj :", obj.Show, "\n")
                }
                if (all(m.expo == m.expo.pre)) break
            }
            if (verbose > 2) {
                cat("num of iterations:", num.iter, "\n")
                cat("m:", m.expo, "\n")
                cat("Loss term: ", eval_loss(y, Xcat, Z, hat.beta, hat.gamma, tau0), "\n")
            }
            hat.beta.mat[, i + 1, iter2] <- hat.beta
            hat.gamma.mat[, i + 1, iter2] <- hat.gamma
            cvErr[i, iter2] <- eval_loss(yval, Xvalcat, Zval, hat.beta, hat.gamma, tau0)
            m.list[, i, iter2] <- m.expo
            
            if (verbose > 1) cat("Cross-validation error:", cvErr[i], "\n")
        }
    }

    cat("\n\nBy the End:", "\n")
    cat(cvErr, "\n\n")
    indices <- which(cvErr == min(cvErr), arr.ind = TRUE)
    if (nrow(indices) == 1) {
        inds <- indices
    } else { # TODO: CHECK THIS
        cat("Warning: Multiple minima detected!\n")
        sorted_indices <- indices[order(indices[,2], -indices[,1]), ]
        inds <- sorted_indices[1, ]
    }

    lam1.tuned <- lam1.vec[inds[1]]
    lam2.tuned <- lam2.vec[inds[2]]
    hat.beta <- hat.beta.mat[, inds[1] + 1, inds[2]]
    hat.gamma <- hat.gamma.mat[, inds[1] + 1, inds[2]]
    m.expo <- m.list[, inds[1], inds[2]]

    m.loc <- rep(0, num.expo)
    for (jjj in 1 : length(m.expo)) {
        if (m.expo[jjj] == p.expo[jjj]) {
            m.loc[jjj] <- m.expo[jjj]
        } else if (m.expo[jjj] == 0) {
            m.loc[jjj] <- 1
        } else if (hat.beta[pos.Mat[jjj, 1] + m.expo[jjj] - 1] > hat.beta[pos.Mat[jjj, 1] + m.expo[jjj]]) {
            m.loc[jjj] <- m.expo[jjj]
        } else {
            m.loc[jjj] <- m.expo[jjj] + 1
        }
    }


    cat("===========================================\n")
    cat("Selected Lambda1: ", lam1.tuned, "\n")
    cat("Selected Lambda2: ", lam2.tuned, "\n")
    cat("Estimated m: ", m.loc, "\n")
    cat("===========================================\n\n")

    return (list(
        "cvErr" = cvErr, "m.loc" = m.loc,
        "lam1.tuned" = lam1.tuned, "lam2.tuned" = lam2.tuned, 
        "hat.beta.mat" = hat.beta.mat, "hat.gamma.mat" = hat.gamma.mat,
        "hat.beta" = hat.beta, "hat.gamma" = hat.gamma,
        "hat.beta.qr" = hat.beta.qr, "hat.gamma.qr" = hat.gamma.qr))
}


# -----------------
# Known m
# -----------------
cv.q.unimodal0 <- function(y, X, Z, yval, Xval, Zval, lam1.vec, lam2.vec, tau0, m.expo, type = 'smooth', numFolds = 5, verbose = 0) {
    n <- length(y); nval <- length(yval)
    p2 <- ncol(Z)
    num.expo <- length(X)
    p.expo <- rep(0, num.expo)
    for (i in 1 : num.expo) {
        p.expo[i] <- ncol(X[[i]])
    }
    p.aggre <- sum(p.expo)

    nlam1 <- length(lam1.vec); nlam2 <- length(lam2.vec)
    hat.beta.mat <- array(0, dim = c(p.aggre, nlam1 + 1, nlam2))
    hat.gamma.mat <- array(0, dim = c(p2, nlam1 + 1, nlam2))
    cvErr <- matrix(0, nrow = nlam1, ncol = nlam2)

    # Position matrix
    pos.Mat <- matrix(nrow = num.expo, ncol = 2)
    start_pos <- 1
    for(i in 1:num.expo) {
        end_pos <- start_pos + p.expo[i] - 1
        pos.Mat[i, ] <- c(start_pos, end_pos)
        start_pos <- end_pos + 1
    }

    if (type == 'smooth') {
        D.list <- list()
        for (jjj in 1 : num.expo) {
            D.list[[jjj]] <- Dmat(p.expo[jjj])
        }
        D <- diagonal_matrix_bind(D.list)
    } else if (type == 'ridge') {
        D <- diag(p.aggre)
    } else {
        stop("'type' must be either 'smooth' or 'ridge'")
    }

    Xcat <- NULL
    Xvalcat <- NULL
    for (jjj in 1 : num.expo) {
        Xcat <- cbind(Xcat, X[[jjj]])
        Xvalcat <- cbind(Xvalcat, Xval[[jjj]])
    }

    xtest <- list(); xtrain <- list()
    ztest <- list(); ztrain <- list()
    ytest <- list(); ytrain <- list()

    rqfit0 <- rq(y ~ Xcat + Z - 1, tau = tau0)
    hat.beta.qr <- rqfit0$coef[1 : ncol(Xcat)]
    hat.gamma.qr <- rqfit0$coef[(ncol(Xcat) + 1) : (ncol(Xcat) + p2)]

    for (iter2 in 1 : nlam2) {
        if (verbose > 0) {
            cat("\n\n========================================\n")
            cat(iter2, "th lam2:", lam2.vec[iter2], "\n")
            cat("========================================\n")
        }
        lam2 <- lam2.vec[iter2]

        hat.beta.mat[, 1, iter2] <- hat.beta.qr
        hat.gamma.mat[, 1, iter2] <- hat.gamma.qr


        for (i in 1 : nlam1) {
            lam1 <- lam1.vec[i]
            if (verbose > 1) {
                cat("\n+++++++++++++++++++++++++++++\n")
                cat(i, "th lam1:", lam1, "\n")
                cat("+++++++++++++++++++++++++++++\n")
            }

            hat.beta <- hat.beta.mat[, i, iter2]
            hat.gamma <- hat.gamma.mat[, i, iter2]
            hat.rr <- y - Xcat %*% hat.beta - Z %*% hat.gamma
            if (verbose > 2) {
                cat("maximal loc: ", which.max(hat.beta), "\n")
                cat("m start:", m.expo, '\n')
            }
            obj.pre <- Inf
            obj <- eval_obj_uni(y, Xcat, Z, hat.beta, hat.gamma, tau0, m.expo, lam1, lam2, D, pos.Mat)
            num.iter <- 0
            # while (abs(obj.pre - obj) > 1e-6 && num.iter < 1e3) {
                tmp <- runADMM.uni(y, Xcat, Z, hat.beta, hat.gamma, hat.rr, tau0, m.expo, sigma = 1e-4, lam1, lam2, D, pos.Mat)
                hat.beta <- tmp$beta
                hat.gamma <- tmp$gamma
                hat.rr <- tmp$rr
                num.iter <- num.iter + 1

            if (verbose > 2) {
                cat("num of iterations:", num.iter, "\n")
                cat("m:", m.expo, "\n")
                cat("Loss term: ", eval_loss(y, Xcat, Z, hat.beta, hat.gamma, tau0), "\n")
            }
            hat.beta.mat[, i + 1, iter2] <- hat.beta
            hat.gamma.mat[, i + 1, iter2] <- hat.gamma
            cvErr[i, iter2] <- eval_loss(yval, Xvalcat, Zval, hat.beta, hat.gamma, tau0)
            
            if (verbose > 1) cat("Cross-validation error:", cvErr[i], "\n")
        }
    }

    cat("\n\nBy the End:", "\n")
    cat(cvErr, "\n\n")
    indices <- which(cvErr == min(cvErr), arr.ind = TRUE)
    if (nrow(indices) == 1) {
        inds <- indices
    } else {
        cat("Warning: Multiple minima detected!\n")
        sorted_indices <- indices[order(indices[,2], -indices[,1]), ]
        inds <- sorted_indices[1, ]
    }


    lam1.tuned <- lam1.vec[inds[1]]
    lam2.tuned <- lam2.vec[inds[2]]
    hat.beta <- hat.beta.mat[, inds[1] + 1, inds[2]]
    hat.gamma <- hat.gamma.mat[, inds[1] + 1, inds[2]]
    m.loc <- m.expo

    cat("===========================================\n")
    cat("Selected Lambda1: ", lam1.tuned, "\n")
    cat("Selected Lambda2: ", lam2.tuned, "\n")
    cat("Estimated m: ", m.loc, "\n")
    cat("===========================================\n\n")

    return (list(
        "cvErr" = cvErr, "m.loc" = m.loc,
        "lam1.tuned" = lam1.tuned, "lam2.tuned" = lam2.tuned, 
        "hat.beta.mat" = hat.beta.mat, "hat.gamma.mat" = hat.gamma.mat,
        "hat.beta" = hat.beta, "hat.gamma" = hat.gamma,
        "hat.beta.qr" = hat.beta.qr, "hat.gamma.qr" = hat.gamma.qr))
}


# ---------------------------------------------
# Fit quantile regression with concave penalty
# ---------------------------------------------
# y : n * 1 response 
# X : list of length num.expo. 
# Z : n * p2 covariates
# lam1.vec : tuning parameters for Unimodal
# lam2.vec : tuning parameters for Smoothness
# ---------------------------------------------
cv.q.concave <- function(y, X, Z, yval, Xval, Zval, lam1.vec, lam2.vec, tau0, sigma=1e-4, verbose = 0) {
    n <- length(y); nval <- length(yval)
    p2 <- ncol(Z)
    num.expo <- length(X)
    p.expo <- rep(0, num.expo)
    for (i in 1 : num.expo) {
        p.expo[i] <- ncol(X[[i]])
    }
    p.aggre <- sum(p.expo)

    nlam1 <- length(lam1.vec); nlam2 <- length(lam2.vec)
    hat.beta.mat <- array(0, dim = c(p.aggre, nlam1 + 1, nlam2))
    hat.gamma.mat <- array(0, dim = c(p2, nlam1 + 1, nlam2))
    cvErr <- matrix(0, nrow = nlam1, ncol = nlam2)

    # Position matrix
    pos.Mat <- matrix(nrow = num.expo, ncol = 2)
    start_pos <- 1
    for(i in 1:num.expo) {
        end_pos <- start_pos + p.expo[i] - 1
        pos.Mat[i, ] <- c(start_pos, end_pos)
        start_pos <- end_pos + 1
    }

    D.list <- list()
    for (jjj in 1 : num.expo) {
        D.list[[jjj]] <- Dmat(p.expo[jjj])
    }
    D <- diagonal_matrix_bind(D.list)

    Xcat <- NULL
    Xvalcat <- NULL
    for (jjj in 1 : num.expo) {
        Xcat <- cbind(Xcat, X[[jjj]])
        Xvalcat <- cbind(Xvalcat, Xval[[jjj]])
    }

    # ! ------------------------------------------
    # ! Store frequently used variables ahead
    # ! ------------------------------------------
    Z.tilde <- rbind(Z, matrix(0, nrow = nrow(D), ncol = p2))
    y.tilde <- c(y, rep(0, nrow(D)))
    XZcat <- cbind(Xcat, Z)
    ZtZZt <- solve(crossprod(Z.tilde), t(Z.tilde))
    # ! ------------------------------------------


    rqfit0 <- rq(y ~ Xcat + Z - 1, tau = tau0)
    hat.beta.qr <- rqfit0$coef[1 : ncol(Xcat)]
    hat.gamma.qr <- rqfit0$coef[(ncol(Xcat) + 1) : (ncol(Xcat) + p2)]

    for (iter2 in 1 : nlam2) {
        if (verbose > 0) {
            cat("\n\n========================================\n")
            cat(iter2, "th lam2:", lam2.vec[iter2], "\n")
            cat("========================================\n")
        }
        lam2 <- lam2.vec[iter2]

        hat.beta.mat[, 1, iter2] <- hat.beta.qr
        hat.gamma.mat[, 1, iter2] <- hat.gamma.qr

        # ! ------------------------------------------
        # ! Store frequently used variables ahead
        # ! ------------------------------------------
        Xcat.tilde <- rbind(Xcat, sqrt(2 * lam2 / sigma) * D)
        ZtZZtX <- ZtZZt %*% Xcat.tilde
        PzX.tilde <- Z.tilde %*% ZtZZtX
        Xcat.bar <- Xcat.tilde - PzX.tilde
        eta <- max(eigen(crossprod(Xcat.bar))$values)
        # ! ------------------------------------------

        for (i in 1 : nlam1) {
            lam1 <- lam1.vec[i]
            if (verbose > 1) {
                cat("\n+++++++++++++++++++++++++++++\n")
                cat(i, "th lam1:", lam1, "\n")
                cat("+++++++++++++++++++++++++++++\n")
            }

            hat.beta <- hat.beta.mat[, i, iter2]
            hat.gamma <- hat.gamma.mat[, i, iter2]
            hat.rr <- y - Xcat %*% hat.beta - Z %*% hat.gamma
            if (verbose > 2) {
                obj <- eval_obj_concave(y, Xcat, Z, hat.beta, hat.gamma, tau0,lam1, lam2, D, pos.Mat)
                cat("Initial obj: ", obj, "\n")
            }
            tmp <- runADMM.concave(y, Xcat, Z, 
                            y.tilde, Z.tilde, XZcat, ZtZZt, Xcat.tilde, Xcat.bar, eta,
                            hat.beta, hat.gamma, hat.rr, tau0, sigma, lam1, lam2, D, pos.Mat, verbose=verbose)
            hat.beta <- tmp$beta
            hat.gamma <- tmp$gamma
            hat.rr <- tmp$rr
            if (verbose > 2) {
                obj <- eval_obj_concave(y, Xcat, Z, hat.beta, hat.gamma, tau0,lam1, lam2, D, pos.Mat)
                cat("Final obj: ", obj, "\n")
            }
            hat.beta.mat[, i + 1, iter2] <- hat.beta
            hat.gamma.mat[, i + 1, iter2] <- hat.gamma
            cvErr[i, iter2] <- eval_loss(yval, Xvalcat, Zval, hat.beta, hat.gamma, tau0)

            if (verbose > 1) cat("Cross-validation error:", cvErr[i], "\n")
        }
    }

    cat("\n\nBy the End:", "\n")
    print(cvErr)
    indices <- which(cvErr == min(cvErr), arr.ind = TRUE)
    if (nrow(indices) == 1) {
        inds <- indices
    } else {
        cat("Warning: Multiple minima detected!\n")
        sorted_indices <- indices[order(indices[,2], -indices[,1]), ]
        inds <- sorted_indices[1, ]
    }

    lam1.tuned <- lam1.vec[inds[1]]
    lam2.tuned <- lam2.vec[inds[2]]
    hat.beta <- hat.beta.mat[, inds[1] + 1, inds[2]]
    hat.gamma <- hat.gamma.mat[, inds[1] + 1, inds[2]]


    cat("===========================================\n")
    cat("Selected Lambda1: ", lam1.tuned, "\n")
    cat("Selected Lambda2: ", lam2.tuned, "\n")
    cat("===========================================\n\n")

    return (list(
        "cvErr" = cvErr,
        "lam1.tuned" = lam1.tuned, "lam2.tuned" = lam2.tuned, 
        "hat.beta.mat" = hat.beta.mat, "hat.gamma.mat" = hat.gamma.mat,
        "hat.beta" = hat.beta, "hat.gamma" = hat.gamma,
        "hat.beta.qr" = hat.beta.qr, "hat.gamma.qr" = hat.gamma.qr))
}

