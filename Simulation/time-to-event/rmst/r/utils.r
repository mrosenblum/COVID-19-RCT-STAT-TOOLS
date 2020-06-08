gendata <- function(data, n, eff) {
    df   <- data[sample(nrow(data), n, replace = TRUE), ]
    df$A <- rbinom(n, 1, 0.5)
    K <- max(df$T)
    df$T[df$D == 1 & df$A == 1] <- df$T[df$D == 1 & df$A == 1] +
        round(rchisq(length(df$T[df$D == 1 & df$A == 1]), df = eff), 0)
    df$T <- pmin(df$T, K)
    df$id <- 1:n
    censored <- runif(n) < 0.05
    C <- sample(1:max(df$T), n, replace = TRUE)
    df$T[C < df$T & censored == 1] <- C[C < df$T & censored == 1]
    df$D[C < df$T & censored == 1] <- 0
    return(df)
}

true <- function(data, eff, tau) {
    ## df <- gendata(data, n, eff)
    ## dlong <- transformData(df, 1)
    ## unad <- unadjusted(dlong, tau)
    ## return(diff(unad$km))
    mean(pmin(dat$T + round(rchisq(length(dat$T), df = eff), 0), tau) - pmin(dat$T, tau))
}
