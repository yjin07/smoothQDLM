generate_sequence <- function(n, p, m, s) {
    # n: maximal value
    # p: length of sequence
    # m: location of mode
    # s: difference
    before_n <- seq(from = n - (m - 1) * s, by = s, length.out = m - 1)
    after_n <- seq(from = n - s, by = -s, length.out = p - m)
    result <- c(before_n, n, after_n)

    # ---------
    #   Scale
    # ---------
    nmax <- max(result)
    nmin <- min(result)
    result <- result - (nmax + nmin) / 2
    result <- result / max(result) * n
    return(result)
}

generate_sequence2 <- function(n, p, m) {
    # n: maximal value
    # p: length of sequence
    # m: location of mode
    before_n <- seq(from = - n, to = n, length.out = m)
    after_n <- seq(from = n, to = -n, length.out = p - m + 1)
    result <- c(before_n, after_n[-1])
    return(result)
}

generate_sequence2s <- function(n, p, m, epsilon = 0.3) {
    # n: maximal value
    # p: length of sequence
    # m: location of mode
    before_n <- seq(from = - n, to = n, length.out = m)
    after_n <- seq(from = n, to = -n, length.out = p - m + 1)
    res <- c(before_n, after_n[-1])
    epsilon <- epsilon * n
    res[which(abs(res) < epsilon)] <- 0
    return(res)
}


generate_sequence3 <- function(n, p, m) {
    set.seed(m)
    increase_diffs <- runif(m - 1, min = 0, max = 1)
    increase_diffs <- increase_diffs / sum(increase_diffs) * n * 2
    increase_values <- c(-n, -n + cumsum(increase_diffs))

    decrease_diffs <- runif(p - m, min = 0, max = 1)
    decrease_diffs <- decrease_diffs / sum(decrease_diffs) * n * 2
    decrease_values <- n - cumsum(decrease_diffs)
    return(c(increase_values, decrease_values))
}


generate_sequence3s <- function(n, p, m, epsilon = 0.3) {
    set.seed(m)
    increase_diffs <- runif(m - 1, min = 0, max = 1)
    increase_diffs <- increase_diffs / sum(increase_diffs) * n * 2
    increase_values <- c(-n, -n + cumsum(increase_diffs))

    decrease_diffs <- runif(p - m, min = 0, max = 1)
    decrease_diffs <- decrease_diffs / sum(decrease_diffs) * n * 2
    decrease_values <- n - cumsum(decrease_diffs)
    res <- c(increase_values, decrease_values)
    epsilon <- epsilon * n
    res[which(abs(res) < epsilon)] <- 0
    return(res)
}

generate_sequence4 <- function(n, p, m, epsilon = 0.3) {
    set.seed(m)
    increase_diffs <- runif(m - 1, min = 0, max = 1)
    increase_diffs <- increase_diffs / sum(increase_diffs) * n * 2
    increase_values <- c(-n, -n + cumsum(increase_diffs))

    decrease_diffs <- runif(p - m, min = 0, max = 1)
    decrease_diffs <- decrease_diffs / sum(decrease_diffs) * n * 2
    decrease_values <- n - cumsum(decrease_diffs)
    res <- c(increase_values, decrease_values)
    epsilon <- epsilon * n
    res[which(abs(res) < epsilon)] <- 0
    return(res)
}

generate_sequence4s <- function(n, p, m, epsilon = 0.3) {
    set.seed(m)
    increase_diffs <- runif(m - 1, min = 0, max = 1)
    increase_diffs <- increase_diffs / sum(increase_diffs) * n * 2
    increase_values <- c(-n, -n + cumsum(increase_diffs))

    decrease_diffs <- runif(p - m, min = 0, max = 1)
    decrease_diffs <- decrease_diffs / sum(decrease_diffs) * n * 2
    decrease_values <- n - cumsum(decrease_diffs)
    res <- c(increase_values, decrease_values)
    epsilon <- epsilon * n
    res[which(abs(res) < epsilon)] <- 0
    res[which(res < 0)] <- 0
    return(res)
}

generate_sequence5 <- function(n, p, m) {
    res <- generate_sequence3(n, p, m)
    res[1 : round(0.6 * m)] <- res[round(0.6 * m)]
    return(res)
}

generate_sequence6 <- function(n, p, m) {
    res <- generate_sequence3(n, p, m)
    res[1 : round(0.4 * m)] <- res[round(0.4 * m)]
    res[(m + round(0.6 * (p - m))) : p] <- res[m + round(0.6 * (p - m))]
    return(res)
}


# * ----------------------------
# * Parabolic curve
# * ----------------------------
generate_parabolic_sequence <- function(n, p, m) {
    # n: maximal value
    # p: length of sequence
    # m: location of mode
    
    # Generate indices for the sequence
    x_before <- seq(from = -1, to = 0, length.out = m)
    x_after <- seq(from = 0, to = 1, length.out = p - m + 1)
    
    # Generate parabolic values for each part
    y_before <- n * (x_before^2)
    y_after <- n * (x_after^2)
    
    # Combine the two parts
    result <- -2 * (c(y_before, y_after[-1]) - n / 2)
    return(result)
}


Dmat <- function(p) {
    res <- matrix(0, nrow = p - 2, ncol = p)
    for (i in 1 : (p - 2)) {
        res[i, i] <- -1
        res[i, i + 1] <- 2
        res[i, i + 2] <- -1
    }
    return(res)
}

tranMat <- function(p) {
    res <- diag(p)
    for (j in 2 : p) {
        res[j, j - 1] <- -1
    }
    return(res)
}

diagonal_matrix_bind <- function(matrix_list) {
    total_rows <- sum(sapply(matrix_list, nrow))
    total_cols <- sum(sapply(matrix_list, ncol))

    big_matrix <- matrix(0, nrow = total_rows, ncol = total_cols)

    start_row <- 1
    start_col <- 1
    for (m in matrix_list) {
        rows <- nrow(m)
        cols <- ncol(m)

        big_matrix[start_row:(start_row + rows - 1), start_col:(start_col + cols - 1)] <- m

        start_row <- start_row + rows
        start_col <- start_col + cols
    }

    return(big_matrix)
}



covAR <- function(p, rho) {
    times <- 1:p
    sigma <- 1
    H <- abs(outer(times, times, "-"))
    V <- sigma * rho^H
    p <- nrow(V)
    V[cbind(1:p, 1:p)] <- V[cbind(1:p, 1:p)] * sigma
    return (V)
}


pinball_loss <- function(y, yhat, tau) {
    n <- length(y)
    diff <- y - yhat
    res <- sapply(diff, function(x) x * (tau - (x < 0)))
    return(sum(res) / n)
}


chkArray <- function(arr) {
    n <- length(arr)

    # 如果数组长度小于2，它就自然满足条件
    if (n < 2) {
        return(TRUE)
    }

    peak_found <- FALSE
    for (i in 2:n) {
        if (arr[i] < arr[i - 1]) {
            # 找到一个下降点
            peak_found <- TRUE
        } else if (peak_found && arr[i] > arr[i - 1]) {
            # 在下降点之后发现上升，返回FALSE
            return(FALSE)
        }
    }

    # 如果整个数组是非递减或非递增的，返回TRUE
    return(TRUE)
}

chkArrayAll <- function(arrs) {
    res <- apply(arrs, 2, chkArray)
    return(all(res))
}