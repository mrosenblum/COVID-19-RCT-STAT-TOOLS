rm(list = ls())
## library(devtools)
## install_github('idiazst/survtmlerct')
library(dplyr)
library(lubridate)
library(survtmlerct)
library(glm2)

load('dat.rds')
source('utils.r')

tau <- 7

set.seed(6235)
rep <- 10000
seeds <- sample(928397, rep)
tasks <- expand.grid(n = c(100, 200, 500, 1000), seed = seeds, eff = c(0, 2, 4))

sim <- function(i) {

    n <- tasks[i, 'n']
    eff <- tasks[i, 'eff']
    seed <- tasks[i, 'seed']

    set.seed(seed)

    data <- gendata(dat, n, eff)
    dlong <- transformData(data, 1)

    fitL <- glm(Lm ~ A * (m + sex + age + o2 + dyspnea + hyper + bilat),
                data = dlong, subset = Im == 1, family = binomial())
    fitR <- glm(Rm ~ A * (as.factor(m) + sex + age + o2 + dyspnea + hyper + bilat),
                data = dlong, subset = Jm == 1, family = binomial())
    fitA <- glm(A ~ sex + age + o2 + dyspnea + hyper + bilat,
                data = dlong, subset = m == 1, family = binomial())

    dlong <- mutate(
        dlong,
        gR1 = bound01(predict(fitR, newdata = mutate(dlong, A = 1), type = 'response')),
        gR0 = bound01(predict(fitR, newdata = mutate(dlong, A = 0), type = 'response')),
        h1  = bound01(predict(fitL, newdata = mutate(dlong, A = 1), type = 'response')),
        h0  = bound01(predict(fitL, newdata = mutate(dlong, A = 0), type = 'response')),
        gA1 = bound01(predict(fitA, newdata = mutate(dlong, A = 1), type = 'response')))

    tmle <- tmle_prob(dlong, tau)
    unad <- unadjusted_prob(dlong, tau)

    return(data.frame(estimator = c('tmle', 'km'),
                      estimate  = c(diff(tmle$prob), diff(unad$prob)),
                      se        = c(tmle$std.error.diff, unad$std.error.diff),
                      eff = eff,
                      seed = seed,
                      n = n))
}

funslave <- function(j){
    index <- (1:nrow(tasks))[(0:(nrow(tasks)-1)) %/% (nrow(tasks) / 500) + 1 == j]
    out <- lapply(index, function(k){
        cat('doing task ', k, ' of ', nrow(tasks), '\n', file = 'progress.txt')
        out <- try(sim(k))
        if(inherits(out, 'try-error'))
            cat('error in task ', k, '\n', file = 'errors.txt')
        return(out)
    })

    return(out)

}
