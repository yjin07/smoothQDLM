# ! -----------------------
# ! Multiple Exposure DLM
# ! -----------------------
library(neariso)

iso.update <- function(w, lambda, mm) {
    tt <- length(w)
    if (mm == 0) {
        res <- -neariso(-c(w), lambda = lambda)$beta
    } else if (mm == tt) {
        res <- neariso(c(w), lambda = lambda)$beta
    } else if (mm == 1) {
        res <- -neariso(-w[(mm + 1) : tt], lambda = lambda)$beta
        res <- c(w[1], res)
    } else if (mm == tt - 1) {
        res <- neariso(w[1 : mm], lambda = lambda)$beta
        res <- c(res, w[tt])
    } else {
        res1 <- neariso(w[1 : mm], lambda = lambda)$beta
        res2 <- -neariso(- w[(mm + 1) : tt], lambda = lambda)$beta
        res <- c(res1, res2)
    }
    return(res)
}

# ------------------------
# Update beta, gamma
# ------------------------
BG.update.uni <- function(y.tilde, X.bar, Z.tilde, ZtZZt, beta, rr, u, eta, sigma, lam1, lam2, m.vec, pos.Mat) {
    n <- length(y.tilde)
    nrr <- length(rr)
    p.aggre <- ncol(X.bar)
    rr.tilde <- c(rr, rep(0, n - nrr))
    u.tilde <- c(u, rep(0, n - nrr))
    yru <- y.tilde - rr.tilde + u.tilde / sigma
    ZtZZtyru <- ZtZZt %*% yru
    t.bar <- yru - Z.tilde %*% ZtZZtyru
    w <- beta + crossprod(X.bar, t.bar - X.bar %*% beta) / eta
    lam0 <- lam1 / (sigma * eta)
    res <- rep(0, p.aggre)
    for (jjj in 1 : nrow(pos.Mat)) {
        mm <- m.vec[jjj]
        res[pos.Mat[jjj, 1] : pos.Mat[jjj, 2]] <- iso.update(w[pos.Mat[jjj, 1] : pos.Mat[jjj, 2]], lambda = lam0, mm = mm)
    }
    return(list("beta" = res, "rr.tilde" = rr.tilde, "u.tilde" = u.tilde))
}

# ---------------------------
# ADMM for solving one ******
# iteration of \beta ********
# ---------------------------
runADMM.uni <- function(y, X, Z, 
                    y.tilde, Z.tilde, XZ, ZtZZt, X.tilde, X.bar, eta,
                    beta.pre, gamma.pre, rr.pre, tau0, m.vec, sigma, lam1, lam2, D, pos.Mat, verbose = 0, max.iter = 1e5, eps1 = 1e-4, eps2 = 1e-4) {
    n <- length(y)
    y.norm <- sqrt(sum(y ^ 2))
    p1 <- ncol(X); p2 <- ncol(Z)
    u.pre <- - sigma * (X %*% beta.pre + Z %*% gamma.pre + rr.pre - y)
    
    for (lll in 1 : max.iter) {
        tmp <- BG.update.uni(y.tilde, X.bar, Z.tilde, ZtZZt, beta.pre, rr.pre, u.pre, eta, sigma, lam1, lam2, m.vec, pos.Mat)
        beta <- tmp$beta
        gamma <- solve(crossprod(Z.tilde), t(Z.tilde)) %*% (y.tilde - tmp$rr.tilde + tmp$u.tilde / sigma - X.tilde %*% beta)
        rr <- rr_update(y, X, Z, beta, gamma, u.pre, sigma, tau0)
        u <- u.pre - sigma * (X %*% beta + Z %*% gamma + rr - y)

        rr.diff <- sqrt(sum((X %*% beta + Z %*% gamma + rr - y) ^ 2))
        ss.diff <- sigma * sqrt(sum(crossprod(XZ, rr - rr.pre) ^ 2))

        epri <- sqrt(n) * eps1 + eps2 * max(sqrt(sum((X %*% beta + Z %*% gamma) ^ 2)), sqrt(sum(rr ^ 2)), y.norm)
        edual <- sqrt(p1 + p2) * eps1 + eps2 * sqrt(sum(crossprod(XZ, u) ^ 2))

        if (verbose == 3) {
            ss.diff.show <- sprintf("%.16f", ss.diff)
            cat("-> Num of iter: ", lll, " sigma:", sigma, " rr: ", rr.diff, " ss: ", ss.diff.show, "\n")
            cat("-> Tolerence: ", epri, edual, "\n")
        }
        
        if (rr.diff < epri && ss.diff < edual || lll == max.iter) {
            # cat("** Num of ADMM inner iter: ", lll, " rr: ", rr.diff, " ss: ", ss.diff, "\n")
            break
        }

        beta.pre <- beta
        gamma.pre <- gamma
        rr.pre <- rr
        u.pre <- u

    }
    return(list(beta = beta, gamma = gamma, rr = rr))
}



# TODO: Wed May 29 14:17:50 EDT 2024: Compare with `ADMM-v5`

# ! -----------------------
# ! Multiple Exposure DLM
# ! -----------------------
# ------------------------
# Update beta, gamma
# ------------------------
BG.update.concave <- function(y.tilde, X.bar, Z.tilde, ZtZZt, beta, rr, u, eta, sigma, lam1, lam2, pos.Mat) {
    n <- length(y.tilde)
    nrr <- length(rr)
    p.aggre <- ncol(X.bar)
    rr.tilde <- c(rr, rep(0, n - nrr))
    u.tilde <- c(u, rep(0, n - nrr))
    yru <- y.tilde - rr.tilde + u.tilde / sigma
    ZtZZtyru <- ZtZZt %*% yru
    t.bar <- yru - Z.tilde %*% ZtZZtyru
    w <- beta + crossprod(X.bar, t.bar - X.bar %*% beta) / eta
    lam0 <- lam1 / (sigma * eta)
    res <- rep(0, p.aggre)
    for (jjj in 1 : nrow(pos.Mat)) {
        res[pos.Mat[jjj, 1] : pos.Mat[jjj, 2]] <- nearConcaveCpp(w[pos.Mat[jjj, 1] : pos.Mat[jjj, 2]], lambda = lam0, rho = 10)$beta
    }
    return(list("beta" = res, "rr.tilde" = rr.tilde, "u.tilde" = u.tilde))
}

# ---------------------------
# ADMM for solving one ******
# iteration of \beta ********
# ---------------------------
runADMM.concave<- function(y, X, Z, 
                    y.tilde, Z.tilde, XZ, ZtZZt, X.tilde, X.bar, eta,
                    beta.pre, gamma.pre, rr.pre, tau0, sigma, lam1, lam2, D, pos.Mat, verbose = 0, max.iter = 1e5, eps1 = 1e-4, eps2 = 1e-4) {
    n <- length(y)
    y.norm <- sqrt(sum(y ^ 2))
    p1 <- ncol(X); p2 <- ncol(Z)
    u.pre <- - sigma * (X %*% beta.pre + Z %*% gamma.pre + rr.pre - y)
    
    for (lll in 1 : max.iter) {
        tmp <- BG.update.concave(y.tilde, X.bar, Z.tilde, ZtZZt, beta.pre, rr.pre, u.pre, eta, sigma, lam1, lam2, pos.Mat)
        beta <- tmp$beta
        gamma <- solve(crossprod(Z.tilde), t(Z.tilde)) %*% (y.tilde - tmp$rr.tilde + tmp$u.tilde / sigma - X.tilde %*% beta)
        rr <- rr_update(y, X, Z, beta, gamma, u.pre, sigma, tau0)
        u <- u.pre - sigma * (X %*% beta + Z %*% gamma + rr - y)

        rr.diff <- sqrt(sum((X %*% beta + Z %*% gamma + rr - y) ^ 2))
        ss.diff <- sigma * sqrt(sum(crossprod(XZ, rr - rr.pre) ^ 2))

        epri <- sqrt(n) * eps1 + eps2 * max(sqrt(sum((X %*% beta + Z %*% gamma) ^ 2)), sqrt(sum(rr ^ 2)), y.norm)
        edual <- sqrt(p1 + p2) * eps1 + eps2 * sqrt(sum(crossprod(XZ, u) ^ 2))

        if (verbose == 3) {
            ss.diff.show <- sprintf("%.16f", ss.diff)
            cat("-> Num of iter: ", lll, " sigma:", sigma, " rr: ", rr.diff, " ss: ", ss.diff.show, "\n")
            cat("-> Tolerence: ", epri, edual, "\n")
        }
        
        if (rr.diff < epri && ss.diff < edual || lll == max.iter) {
            # cat("** Num of ADMM inner iter: ", lll, " rr: ", rr.diff, " ss: ", ss.diff, "\n")
            break
        }

        beta.pre <- beta
        gamma.pre <- gamma
        rr.pre <- rr
        u.pre <- u

    }
    return(list(beta = beta, gamma = gamma, rr = rr))
}