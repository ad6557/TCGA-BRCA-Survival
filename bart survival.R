require(survbart)

## load survival package for the advanced lung cancer example
require(survival)

group <- -which(is.na(lung[ , 7])) ## remove missing row for ph.karno
times <- lung[group, 2]   ##lung$time
delta <- lung[group, 3]-1 ##lung$status: 1=censored, 2=dead
##delta: 0=censored, 1=dead

## this study reports time in days rather than months like other studies
## coarsening from days to months will reduce the computational burden
times <- ceiling(times/30)

summary(times)
table(delta)

x.train <- as.matrix(lung[group, c(4, 5, 7)]) ## matrix of observed covariates

## lung$age:        Age in years
## lung$sex:        Male=1 Female=2
## lung$ph.karno:   Karnofsky performance score (dead=0:normal=100:by=10)
##                  rated by physician

dimnames(x.train)[[2]] <- c('age(yr)', 'M(1):F(2)', 'ph.karno(0:100:10)')

summary(x.train[ , 1])
table(x.train[ , 2])
table(x.train[ , 3])

x.test <- matrix(nrow=84, ncol=3) ## matrix of covariate scenarios

dimnames(x.test)[[2]] <- dimnames(x.train)[[2]]

i <- 1

for(age in 5*(9:15)) for(sex in 1:2) for(ph.karno in 10*(5:10)) {
  x.test[i, ] <- c(age, sex, ph.karno)
  i <- i+1
}

## run one long MCMC chain in one process
set.seed(99)
post <- surv.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)

## run "mc.cores" number of shorter MCMC chains in parallel processes
## post <- mc.surv.bart(x.train=x.train, times=times, delta=delta, x.test=x.test,
##                      mc.cores=20, seed=99)

##saveRDS(object=post, file='post.rds')
## you can save time by reading in the posterior
## instead of re-generating it every time
## post <- readRDS(file='post.rds')

## let's look at some survival curves
## first, a younger group with a healthier KPS
## age 50 with KPS=90: males and females
## males: row 17, females: row 23
x.test[c(17, 23), ]

# length(post[["surv.test.mean"]]) = 2604 = 84 x 31
# K = length(unique(times)) = 31

# the 17th patient: 16x31+1 -- 17x31
low.risk.males <- 16*post$K+1:post$K ## K=unique times including censoring

low.risk.females <- 22*post$K+1:post$K

## second, an older group with a poor KPS
## age 70 with KPS=60: males and females
x.test[c(62, 68), ]

high.risk.males <- 61*post$K+1:post$K

high.risk.females <- 67*post$K+1:post$K

old.par <- par(mfrow=c(2, 2))

plot(post$times, post$surv.test.mean[low.risk.males], type='s', col='blue',
     main='Age 50 with KPS=90',
     xlab='t', ylab='S(t)', ylim=c(0, 1))
points(post$times, post$surv.test.mean[low.risk.females], type='s', col='red')

plot(post$times, post$surv.test.mean[high.risk.males], type='s', col='blue',
     main='Age 70 with KPS=60',
     xlab='t', ylab='S(t)', ylim=c(0, 1))
points(post$times, post$surv.test.mean[high.risk.females], type='s', col='red')

plot(post$times, post$surv.test.mean[low.risk.males], type='s', col='blue',
     main='Males',
     xlab='t', ylab='S(t)', ylim=c(0, 1))
points(post$times, post$surv.test.mean[high.risk.males], type='s')

plot(post$times, post$surv.test.mean[low.risk.females], type='s', col='red',
     main='Females',
     xlab='t', ylab='S(t)', ylim=c(0, 1))
points(post$times, post$surv.test.mean[high.risk.females], type='s')

par(old.par)

## calculate 95% credible interval
credible95 <- apply(post$surv.test[ , low.risk.females], 2, quantile,
                    probs=c(0.025, 0.975))

plot(post$times, post$surv.test.mean[low.risk.females], type='s', col='red',
     main='Females aged 50/KPS=90',
     xlab='t', ylab='S(t)', ylim=c(0, 1))
points(post$times, credible95[1, ], type='s')
points(post$times, credible95[2, ], type='s')
