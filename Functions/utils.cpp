#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
// [[Rcpp::depends("RcppArmadillo")]]


double rhoTau(double x, double tau) {
    if (x > 0) {
        return tau * x;
    } else {
        return (tau - 1) * x;
    }
}

// [[Rcpp::export]]
double eval_loss(const arma::vec& y, 
                const arma::mat& X, 
                const arma::mat& Z, 
                const arma::vec& beta, 
                const arma::vec& gamma, 
                double tau0) {
    int n = y.n_elem;
    arma::vec residuals = y - X * beta - Z * gamma;
    
    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        res += rhoTau(residuals[i], tau0);
    }

    return res / n;
}

// eval.unimodal
// [[Rcpp::export]]
double eval_unimodal(const arma::vec& beta, int m) {
    double res = 0;
    int tt = beta.n_elem;

    if (m == 0) {
        res += arma::accu(arma::clamp(beta.subvec(1, tt-1) - beta.subvec(0, tt-2), 0.0, arma::datum::inf));
    } else if (m == tt) {
        res += arma::accu(arma::clamp(beta.subvec(0, tt-2) - beta.subvec(1, tt-1), 0.0, arma::datum::inf));
    } else if (m == 1) {
        res += arma::accu(arma::clamp(beta.subvec(2, tt-1) - beta.subvec(1, tt-2), 0.0, arma::datum::inf));
    } else if (m == tt - 1) {
        res += arma::accu(arma::clamp(beta.subvec(0, tt-3) - beta.subvec(1, tt-2), 0.0, arma::datum::inf));
    } else {
        res += arma::accu(arma::clamp(beta.subvec(0, m-2) - beta.subvec(1, m-1), 0.0, arma::datum::inf));
        res += arma::accu(arma::clamp(beta.subvec(m+1, tt-1) - beta.subvec(m, tt-2), 0.0, arma::datum::inf));
    }

    return res;
}




// eval.obj
// [[Rcpp::export]]
double eval_obj_uni(const arma::vec& y, 
                const arma::mat& X, 
                const arma::mat& Z, 
                const arma::vec& beta, 
                const arma::vec& gamma, 
                double tau0,
                const arma::ivec& m_vec,
                double lam1,
                double lam2,
                const arma::mat& Dmat,
                const arma::mat& pos_Mat) {
    int n = y.n_elem;
    int tt = beta.n_elem;
    double res = eval_loss(y, X, Z, beta, gamma, tau0);

    for (int jjj = 0; jjj < m_vec.n_elem; ++jjj) {
        int start = pos_Mat(jjj, 0) - 1;
        int end = pos_Mat(jjj, 1) - 1;
        res += lam1 * eval_unimodal(beta.subvec(start, end), m_vec[jjj]);
    }

    res += lam2 * arma::accu(arma::square(Dmat * beta));
    return res;
}


double prox(double x, double alpha, double t) {
    if (x > t / alpha) {
        return x - t / alpha;
    } else if (x < (t - 1) / alpha) {
        return x - (t - 1) / alpha;
    } else {
        return 0;
    }
}

// [[Rcpp::export]]
arma::vec rr_update(const arma::vec& y, 
                    const arma::mat& X, 
                    const arma::mat& Z, 
                    const arma::vec& beta, 
                    const arma::vec& gamma, 
                    const arma::vec& u, 
                    double sigma, 
                    double tau) {
    
    int n = y.n_elem;
    arma::vec residuals = y - X * beta - Z * gamma + u / sigma;

    arma::vec res = arma::zeros(n);
    for (int i = 0; i < n; ++i) {
        res[i] = prox(residuals[i], n * sigma, tau);
    }
    return res;
}



// [[Rcpp::export]]
arma::mat DDmat(int p) {
    arma::mat res = arma::zeros<arma::mat>(p - 2, p);

    for (int i = 0; i < p - 2; ++i) {
        res(i, i) = 1;
        res(i, i + 1) = -2;
        res(i, i + 2) = 1;
    }

    return res;
}


// eval.concave
// [[Rcpp::export]]
double eval_concave(const arma::vec& beta, const arma::mat& Dmat) {
    arma::vec theta = Dmat * beta;
    arma::vec pmax_theta = arma::clamp(theta, 0, arma::datum::inf);
    double result = arma::accu(pmax_theta);
    return result;
}


// eval.obj
// [[Rcpp::export]]
double eval_obj_concave(const arma::vec& y, 
                const arma::mat& X, 
                const arma::mat& Z, 
                const arma::vec& beta, 
                const arma::vec& gamma, 
                double tau0,
                double lam1,
                double lam2,
                const arma::mat& Dmat,
                const arma::mat& pos_Mat) {
    int n = y.n_elem;
    int tt = beta.n_elem;
    double res = eval_loss(y, X, Z, beta, gamma, tau0);

    for (int jjj = 0; jjj < pos_Mat.n_rows; ++jjj) {
        int start = pos_Mat(jjj, 0) - 1;
        int end = pos_Mat(jjj, 1) - 1;
        arma::mat dd = DDmat(end - start + 1);
        res += lam1 * eval_concave(beta.subvec(start, end), dd);
    }

    res += lam2 * arma::accu(arma::square(Dmat * beta));
    return res;
}


double prox(double x, double alpha, double t) {
    if (x > t / alpha) {
        return x - t / alpha;
    } else if (x < (t - 1) / alpha) {
        return x - (t - 1) / alpha;
    } else {
        return 0;
    }
}

// [[Rcpp::export]]
arma::vec rr_update(const arma::vec& y, 
                    const arma::mat& X, 
                    const arma::mat& Z, 
                    const arma::vec& beta, 
                    const arma::vec& gamma, 
                    const arma::vec& u, 
                    double sigma, 
                    double tau) {
    
    int n = y.n_elem;
    arma::vec residuals = y - X * beta - Z * gamma + u / sigma;

    arma::vec res = arma::zeros(n);
    for (int i = 0; i < n; ++i) {
        res[i] = prox(residuals[i], n * sigma, tau);
    }
    return res;
}


// [[Rcpp::export]]
Rcpp::List nearConcaveCpp(const arma::vec& y, double lambda, 
                        Rcpp::Nullable<arma::vec> beta_ = R_NilValue, 
                        Rcpp::Nullable<arma::vec> theta_ = R_NilValue, 
                        double rho = 1, double tol = 1e-6, int max_iter = 1e6, 
                        double mu = 10, double tau_inc = 2, double tau_dec = 2,
                        bool quiet = true) {
    
    int T = y.n_elem;
    arma::mat D = DDmat(T);
    arma::mat DtD = D.t() * D;
    arma::mat Q = arma::inv(arma::eye(T, T) + rho * DtD);
    
    arma::vec beta = beta_.isNull() ? y : Rcpp::as<arma::vec>(beta_);
    arma::vec theta = theta_.isNull() ? D * y : Rcpp::as<arma::vec>(theta_);
    arma::vec gamma(T - 2, arma::fill::zeros);
    
    for (int kk = 0; kk < max_iter; ++kk) {
        beta = Q * (y + rho * D.t() * (theta + gamma));
        arma::vec Dbeta = D * beta;
        arma::vec theta_prev = theta;
        theta = arma::clamp(Dbeta - gamma - lambda / rho, 0, arma::datum::inf);
        theta = arma::min(theta, Dbeta - gamma);
        
        gamma = gamma + theta - Dbeta;
        
        double rk = arma::norm(theta - Dbeta, 2);
        double sk = arma::norm(rho * D.t() * (theta - theta_prev), 2);
        
        // if (!quiet) {
        //     Rcpp::Rcout << sk << " " << rk << std::endl;
        // }
        
        if (rk < tol && sk < tol) {
            // Rcpp::Rcout << "Number of CPP iterations: " << kk << std::endl;
            break;
        }
    }
    
    return Rcpp::List::create(Rcpp::Named("beta") = beta,
                            Rcpp::Named("theta") = theta,
                            Rcpp::Named("gamma") = gamma);
}