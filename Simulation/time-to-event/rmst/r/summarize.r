library(tidyverse)
library(data.table)
library(ggplot2)
library(ggthemes)
library(grid)
library(survtmlerct)

source('utils.r')
source('code.r')

res <- list()
for(i in 1:500){
    ss <- try(load(paste0('../out/out', i, '.rda')))
    if(inherits(ss, 'try-error')) {
        cat('error ', i)
    } else {
        rtemp <- get(paste0('r', i))
        res <- c(res, rtemp)
    }
}

table(sapply(res, length))

out <- rbindlist(res[sapply(res, length) != 1])

out <- out %>% filter(abs(estimate) < 2, se < 10)

load('dat.rds')
true0 <- true(dat, 0, tau)
true1 <- true(dat, 2, tau)
true2 <- true(dat, 4, tau)
true <- data.frame(eff = c(0, 2, 4), truth = c(true0, true1, true2))
save(true, file = 'true.rds')

alpha <- 0.05
dplot <- out %>% left_join(true, by = 'eff') %>%
    group_by(n, estimator, truth) %>%
    summarise(
        power = mean(abs(estimate / se) > qnorm(1 - alpha / 2)),
        mse   = mean((estimate - truth)^2),
        bias  = mean(estimate - truth),
        var   = var(estimate)) %>%
    arrange(n, truth, estimator)

ref <- dplot %>%
    filter(estimator == 'km') %>%
    ungroup() %>%
    select(n, truth, mse) %>%
    rename(ref = mse)

dplot <- dplot %>% left_join(ref, by = c('n', 'truth')) %>%
    mutate(releff = mse / ref) %>%
    select(-ref)

library(knitr)
kable(dplot, 'latex', digits = 3)
dplot
