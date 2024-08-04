library("optparse")
option_list <- list(
    make_option(c("-m", "--model"), type = "character", default = "A", 
                help = "Model to run [A|B|C]", metavar = "character"),
    make_option(c("-e", "--error"), type = "character", default = "normal", 
                help = "Error distribution [normal|t]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!opt$model %in% c("A", "B", "C")) {
    stop("Invalid model. Choose from 'A', 'B', or 'C'.")
}

if (!opt$error %in% c("normal", "t")) {
    stop("Invalid error distribution. Choose from 'normal' or 't'.")
}

print(paste("Model:", opt$model))
print(paste("Error distribution:", opt$error))



library(logger)
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
t0 <- Sys.time(); t0
source("Functions/helper.R") 

combs <- expand.grid(tau = rep(c(0.10, 0.25, 0.3, 0.5, 1.0), each = 50), n = c(500, 750, 1000), quantile = c(0.05, 0.25, 0.5))
final <- list()

tau0 <- combs$quantile[uu]
n <- combs$n[uu]
ntest <- 1000
num.expo <- 6
p.all <- c(30, 30, 30, 30, 30, 30, 30)
m.expo <- c(12, 15, 18, 17, 15, 13)
p2 <- 5

SIGMA = 1e-4

gamma <- rep(1, p2)

# ! ------------------
# ! Generate beta
# ! ------------------
if (opt$model == "A") {
    for (i in 1:num.expo) {
        beta.all[[i]] <- generate_sequence3(5, p.all[i], m.expo[i])
    }
} else if (opt$model == "B") {
    for (i in 1:num.expo) {
        beta.all[[i]] <- generate_sequence3s(5, p.all[i], m.expo[i], 0.5)
} else if (opt$model == "C") {
    for (i in 1:num.expo) {
        beta.all[[i]] <- generate_parabolic_sequence(5, p.all[i], m.expo[i])
    }
}

beta.aggre <- NULL
for (i in 1 : num.expo) {
    beta.aggre <- c(beta.aggre, beta.all[[i]])
}


SigmaX <- list()
SigmaXsqrt <- list()
for (i in 1 : num.expo) {
    SigmaX[[i]] <- covAR(p.all[i], 0.8)
    eo <- eigen(SigmaX[[i]])
    SigmaXsqrt[[i]] <- eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)
}

SigmaX.cat <- diagonal_matrix_bind(SigmaX)

vars <- as.numeric(sqrt(t(beta.aggre) %*% SigmaX.cat %*% beta.aggre / combs$tau[uu]))


# ? ----------------------
# ? set random seed here
# ? ----------------------
set.seed(uu); cat("Seed: ", uu, "\n")
# ----------------------
# Training set
# ----------------------
X <- list()
Xp <- list()
BasisFunctions <- list()
tmp <- rep(0, n)
for (i in 1 : num.expo) {
    X[[i]] <- matrix(rnorm(n * p.all[i]), nrow = n) %*% SigmaXsqrt[[i]]
    tmp <- tmp + X[[i]] %*% beta.all[[i]]
    covx <- cov(X[[i]])
    pca <- refund::fpca.face(Y=covx, knots=10, pve=0.90)
    nBasis <- dim(pca$efunctions)[2]
    BasisFunctions[[i]] <- matrix(0, p.all[i], nBasis + 1)
    BasisFunctions[[i]][, 1] <- rep(1, p.all[i])
    BasisFunctions[[i]][, 2:(nBasis + 1)] <- pca$efunctions[,1:nBasis]
    designMat <- matrix(NA, dim(X[[1]])[1], nBasis + 1)
    designMat[, 1] <- apply(X[[i]], 1, sum)
    for (k in 2 : (nBasis + 1)) {
        designMat[, k] <- apply(t(t(X[[i]])*BasisFunctions[[i]][,k]), 1, sum)
    }
    Xp[[i]] <- designMat
}
Z <- matrix(rnorm(n * (p2 - 1)), nrow = n)
Z <- cbind(rep(1, n), Z)

if (opt$error == "t") {
    Y <- tmp + Z %*% gamma + vars / sqrt(2) * rt(n, 4)
} else {
    Y <- tmp + Z %*% gamma + vars * rnorm(n)
}

# ----------------------
# Validation set
# ----------------------
Xval <- list()
tmp <- rep(0, n)
for (i in 1 : num.expo) {
    Xval[[i]] <- matrix(rnorm(n * p.all[i]), nrow = n) %*% SigmaXsqrt[[i]]
    tmp <- tmp + Xval[[i]] %*% beta.all[[i]]
}
Zval <- matrix(rnorm(n * (p2 - 1)), nrow = n)
Zval <- cbind(rep(1, n), Zval)

if (opt$error == "t") {
    Yval <- tmp + Zval %*% gamma + vars / sqrt(2) * rt(n, 4)
} else {
    Yval <- tmp + Zval %*% gamma + vars * rnorm(n)
}

# ----------------------
# Testing set
# ----------------------
Xtest <- list()
tmp <- rep(0, ntest)
for (i in 1 : num.expo) {
    Xtest[[i]] <- matrix(rnorm(ntest * p.all[i]), nrow = ntest) %*% SigmaXsqrt[[i]]
    tmp <- tmp + Xtest[[i]] %*% beta.all[[i]]
}
Ztest <- matrix(rnorm(ntest * (p2 - 1)), nrow = ntest)
Ztest <- cbind(rep(1, ntest), Ztest)

if (opt$error == "t") {
    Ytest <- tmp + Ztest %*% gamma + vars / sqrt(2) * rt(ntest, 4)
} else {
    Ytest <- tmp + Ztest %*% gamma + vars * rnorm(ntest)
}

Xcat <- NULL
Xpcat <- NULL
Xvalcat <- NULL
Xtest.cat <- NULL
for (jjj in 1 : num.expo) {
    Xcat <- cbind(Xcat, X[[jjj]])
    Xpcat <- cbind(Xpcat, Xp[[jjj]])
    Xvalcat <- cbind(Xvalcat, Xval[[jjj]])
    Xtest.cat <- cbind(Xtest.cat, Xtest[[jjj]])
}


# ----------------------------
# Spline quantile regression
# ----------------------------
p.basis <- rep(0, num.expo)
for (jj in 1 : num.expo) {
    cat(dim(Xp[[jj]]), "\n")
    cat(dim(BasisFunctions[[jj]]), "\n")
    p.basis[jj] <- dim(BasisFunctions[[jj]])[2]
}
rqfit <- rq(Y ~ Z + Xpcat - 1, tau = tau0)
hat.gamma.spline <- rqfit$coef[1:ncol(Z)]
basisCoefs <- rqfit$coef[(ncol(Z) + 1):length(rqfit$coef)]
pos.Mat <- matrix(nrow = num.expo, ncol = 2)
start_pos <- 1
for(i in 1:num.expo) {
    end_pos <- start_pos + p.basis[i] - 1
    pos.Mat[i, ] <- c(start_pos, end_pos)
    start_pos <- end_pos + 1
}
hat.beta.spline.list <- list()
for (jj in 1 : num.expo) {
    hat.beta.spline.list[[jj]] <- BasisFunctions[[jj]] %*% basisCoefs[pos.Mat[jj, 1]:pos.Mat[jj, 2]]
}
hat.beta.spline <- NULL
for (i in 1:length(hat.beta.spline.list)) {
    hat.beta.spline <- c(hat.beta.spline, hat.beta.spline.list[[i]])
}
hat.Ytest.spline <- Xtest.cat %*% hat.beta.spline + Ztest %*% hat.gamma.spline


# ----------------------------
# Spline with Ridge penalty
# ----------------------------
ZXpcat <- cbind(Z, Xpcat)
penalty.factors <- c(rep(0, ncol(Z)), rep(1, ncol(Xpcat)))
nLambda <- 50
lambda.vec <- 10 ^ seq(from = 2, to = -4, length.out = 50)
cvErr.ridge <- rep(0, nLambda)
hat.beta.ridge <- list()
hat.gamma.ridge <- list()
hqfit <- hqreg_raw(ZXpcat, Y, method = "quantile", tau = tau0, alpha = 0, lambda = lambda.vec, intercept = FALSE, penalty.factor = penalty.factors)
for (jj in 1 : nLambda) {
    hat.gamma.ridge[[jj]] <- hqfit$beta[1 : ncol(Z), jj]
    basisCoefs <- hqfit$beta[(ncol(Z) + 1) : (dim(hqfit$beta)[1]), jj]
    hat.beta.ridge[[jj]] <- BasisFunctions[[1]] %*% basisCoefs[pos.Mat[1,1]:pos.Mat[1, 2]]
    for (pp in 2 : num.expo) {
        hat.beta.ridge[[jj]] <- c(hat.beta.ridge[[jj]], BasisFunctions[[pp]] %*% basisCoefs[pos.Mat[pp,1]:pos.Mat[pp, 2]])
    }
    cvErr.ridge[jj] <- eval_loss(Yval, Xvalcat, Zval, hat.beta.ridge[[jj]], hat.gamma.ridge[[jj]], tau0)
}
best.ind <- which.min(cvErr.ridge)
hat.beta.spline.ridge <- hat.beta.ridge[[best.ind]]
hat.gamma.spline.ridge <- hat.gamma.ridge[[best.ind]]
hat.Ytest.spline.ridge <- Xtest.cat %*% hat.beta.spline.ridge + Ztest %*% hat.gamma.spline.ridge

# ----------------------------
# Elastic Net quantile regression
# ----------------------------
ZXcat <- cbind(Z, Xcat)
penalty.factors <- c(rep(0, ncol(Z)), rep(1, ncol(Xcat)))
nLambda <- 50
lambda.vec <- 10 ^ seq(from = 2, to = -4, length.out = 50)
cvErr.elastic <- rep(0, nLambda)
hat.beta.elastic <- list()
hat.gamma.elastic <- list()
hqfit <- hqreg_raw(ZXcat, Y, method = "quantile", tau = tau0, alpha = 0.5, lambda = lambda.vec, intercept = FALSE, penalty.factor = penalty.factors)
for (jj in 1 : nLambda) {
    hat.gamma.elastic[[jj]] <- hqfit$beta[1 : ncol(Z), jj]
    hat.beta.elastic[[jj]] <- hqfit$beta[(ncol(Z) + 1) : (dim(hqfit$beta)[1]), jj]
    cvErr.elastic[jj] <- eval_loss(Yval, Xvalcat, Zval, hat.beta.elastic[[jj]], hat.gamma.elastic[[jj]], tau0)
}
best.ind <- which.min(cvErr.elastic)
hat.beta.elastic <- hat.beta.elastic[[best.ind]]
hat.gamma.elastic <- hat.gamma.elastic[[best.ind]]
hat.Ytest.elastic <- Xtest.cat %*% hat.beta.elastic + Ztest %*% hat.gamma.elastic

# ----------------------------
# Standard QR and our method
# ----------------------------
log_info("Start fitting nearly unimodal quantile regression...")
lambda.vec <- c(0, 10 ^ seq(from = -2,to = 1, length.out = 10)) # Unimodal penalty
lambda2.vec <- 10 ^ seq(from = -4, to = 4, length.out = 10) # Smoothness penalty
res <- cv.q.unimodal(Y, X, Z, Yval, Xval, Zval, lambda.vec, lambda2.vec, tau0, sigma=SIGMA, verbose = 1)
hat.beta.uni <- res$hat.beta
hat.gamma.uni <- res$hat.gamma
hat.Ytest.uni <- Xtest.cat %*% hat.beta.uni + Ztest %*% hat.gamma.uni
hat.beta.qr <- res$hat.beta.qr
hat.gamma.qr <- res$hat.gamma.qr
hat.Ytest.qr <- Xtest.cat %*% hat.beta.qr + Ztest %*% hat.gamma.qr


# ---------------------------
# Ridge quantile regression
# ---------------------------
ZXcat <- cbind(Z, Xcat)
penalty.factors <- c(rep(0, ncol(Z)), rep(1, ncol(Xcat)))
nLambda <- 50
lambda.vec <- 10 ^ seq(from = 2, to = -4, length.out = 50)
cvErr.ridge <- rep(0, nLambda)
hat.beta.ridge <- list()
hat.gamma.ridge <- list()
hqfit <- hqreg_raw(ZXcat, Y, method = "quantile", tau = tau0, alpha = 0, lambda = lambda.vec, intercept = FALSE, penalty.factor = penalty.factors)
for (jj in 1 : nLambda) {
    hat.gamma.ridge[[jj]] <- hqfit$beta[1 : ncol(Z), jj]
    hat.beta.ridge[[jj]] <- hqfit$beta[(ncol(Z) + 1) : (dim(hqfit$beta)[1]), jj]
    cvErr.ridge[jj] <- eval_loss(Yval, Xvalcat, Zval, hat.beta.ridge[[jj]], hat.gamma.ridge[[jj]], tau0)
}
best.ind <- which.min(cvErr.ridge)
hat.beta.ridge <- hat.beta.ridge[[best.ind]]
hat.gamma.ridge <- hat.gamma.ridge[[best.ind]]
hat.Ytest.ridge <- Xtest.cat %*% hat.beta.ridge + Ztest %*% hat.gamma.ridge

# ---------------------------
# Smooth quantile regression
# ! smoothness only
# ---------------------------
lambda.vec <- c(0)
lambda2.vec <- 10 ^ seq(from = -4, to = 4, length.out = 10) # Smoothness penalty
res2 <- cv.q.unimodal(Y, X, Z, Yval, Xval, Zval, lambda.vec, lambda2.vec, tau0, sigma=SIGMA, verbose = 1)
hat.beta.smooth <- res2$hat.beta
hat.gamma.smooth <- res2$hat.gamma
hat.Ytest.smooth <- Xtest.cat %*% hat.beta.smooth + Ztest %*% hat.gamma.smooth


final[['res']] <- res
final[['res2']] <- res2
final[['m.loc.uni']] <- res$m.loc




cat("\nEstimated mode location: ", res$m.loc, "\n")
cat("True parameters of beta: ", "\n")
cat("Quantile regression estimates of beta: ", "\n")
cat("Penalized estimates of beta: ", "\n")
cbind(beta.aggre, hat.beta.qr, hat.beta.uni, hat.beta.ridge, hat.beta.elastic, hat.beta.smooth,
hat.beta.spline, hat.beta.spline.ridge
)

final[['beta.all']] <- cbind(beta.aggre, hat.beta.qr, hat.beta.uni, hat.beta.ridge, hat.beta.elastic, hat.beta.smooth, 
hat.beta.spline, hat.beta.spline.ridge
)


resErr <- rep(0, 16)
resErr[1] <- sqrt(sum((hat.beta.qr - beta.aggre) ^ 2))
resErr[2] <- sqrt(sum((hat.beta.uni - beta.aggre) ^ 2))
resErr[3] <- sqrt(sum((hat.beta.ridge - beta.aggre) ^ 2))
resErr[4] <- sqrt(sum((hat.beta.elastic - beta.aggre) ^ 2))
resErr[5] <- sqrt(sum((hat.beta.smooth - beta.aggre) ^ 2))
# resErr[6] <- sqrt(sum((hat.beta.uni0 - beta.aggre) ^ 2))
resErr[7] <- sqrt(sum((hat.beta.spline - beta.aggre) ^ 2))
resErr[8] <- sqrt(sum((hat.beta.spline.ridge - beta.aggre) ^ 2))

resErr[9] <- pinball_loss(Ytest, hat.Ytest.qr, tau0)
resErr[10] <- pinball_loss(Ytest, hat.Ytest.uni, tau0)
resErr[11] <- pinball_loss(Ytest, hat.Ytest.ridge, tau0)
resErr[12] <- pinball_loss(Ytest, hat.Ytest.elastic, tau0)
resErr[13] <- pinball_loss(Ytest, hat.Ytest.smooth, tau0)
# resErr[14] <- pinball_loss(Ytest, hat.Ytest.uni0, tau0)
resErr[15] <- pinball_loss(Ytest, hat.Ytest.spline, tau0)
resErr[16] <- pinball_loss(Ytest, hat.Ytest.spline.ridge, tau0)


final[['err']] <- resErr

cat("Compare: \n")
resErr

Sys.time()
Sys.time() - t0

final[['time']] <- Sys.time() - t0

cat(paste("Results-v01s/01", opt$model, "-",  opt$error,"/Res_", uu, ".RDS", sep = ""))
saveRDS(final, file = paste("Results-v01s/01", opt$model, "-",  opt$error,"/Res_", uu, ".RDS", sep = ""))
