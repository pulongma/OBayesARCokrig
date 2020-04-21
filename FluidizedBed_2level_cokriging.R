##############################################################################
##############################################################################

rm(list=ls())

############################################################################
## Install the ARCokrig Package
## install relevant packages if not installed 
install.packages(c("Rcpp", "RcppEigen", "RcppArmadillo", "RobustGaSP", "mvtnorm", "fields"))
## install ARCokrig package from local source
install.packages("ARCokrig_0.1.0.tar.gz", repos=NULL, type="source")
## one can also install the package from github
library(devtools)
install_github("pulongma/ARCokrig")
library(ARCokrig)


library(ggplot2)
library(reshape2)
library(plyr)


##############################################################################
##############################################################################
######### Test with Fluidized-bed processes
#rm(list=ls())
#setwd("~/Documents/manuscript/OBayesCokriging/code")

df = matrix(nrow=28, ncol=6)
colnames(df) <- c("Hr", "Tr", "Ta", "Rf", "Pa", "Vf")
df[ ,"Hr"] <- c(51, 46.4, 46.6, 53.1, 52, 45.6, 47.3, 53.3, 44, 52.3, 55, 54,
           50.8, 48, 42.8, 55.7, 55.2, 54.4, 55.4, 52.9, 28.5, 26.1, 24.2,
           25.4, 45.1, 43.1, 42.7, 38.7)
df[ ,"Tr"] <- c(20.7, 21.3, 19.2, 21.1, 20.4, 21.4, 19.5, 21.4, 20.1, 21.6,
                20.2, 20.6, 21.1, 21.2, 22.4, 20.8, 20.7, 20.7, 19.8, 20, 18.3, 
                19, 18.9, 18.5, 19.6, 20.3, 20.4, 21.6)
df[ ,"Ta"] <- c(50, 60, 70, 80, 90, 60, 70, 80, 70, 80,80,80,80,80,80,
                50,50,50,50,50, 80,80,80,80, 50,50,50,50)
df[ ,"Rf"] <- c(5.52, 5.53, 5.53, 5.51, 5.21, 7.25, 7.23, 7.23, 8.93, 8.91,
                7.57, 7.58, 7.40, 7.43, 7.51, 3.17, 3.18, 3.19, 3.20, 3.19,
                7.66, 7.69, 7.69, 7.70, 3.20, 3.23, 3.20, 3.22)
df[ ,"Pa"] <- c(rep(2.5,10), 1, 1.5, 2, 2.5, 3, 1, 1.5, 2, 2.5, 3, rep(2.5,8))
df[ ,"Vf"] <- c(rep(3,21), 4, 4.5, 5, 3, 4, 4.5, 5)
T2exp <- c(30.4, 37.6, 45.1, 50.2, 57.9, 32.9, 39.5, 45.6, 34.2, 41.1, 45.7,
           44.6, 44.7, 44, 43.3, 37, 37.2, 37.1, 36.9, 36.8, 46, 54.7, 57, 
           58.9, 35.9, 40.3, 41.9, 43.1)
T21 <- c(32.4, 39.5, 46.8, 53.8, 61.7, 35.2, 42.4, 49.5, 37.5, 45.5, 50.5,
         49.8, 49.8, 49.2, 48.6, 39.5,39.5,39.5,39.5, 37.7, 48.7, 57.7, 
         60.1, 62, 37.9, 41.7, 43, 43.9)
T22 <- c(31.5, 38.5, 45.5, 52.6, 59.9, 34.6, 41, 48.5, 36.6, 44.3, 49, 48.4,
         48.4, 48, 47.5, 38, 38.5, 37.5, 38.5, 37.2, 47.3, 56.2, 58.7, 60.5,  
         37.1, 40.8, 42.3, 42.3)
T23 <- c(30.2, 37, 43.7, 51, 58.2, 32.6, 39.1, 46.4, 34.8, 42, 47, 46.3, 46.3,
         45.7, 45.4, 37.7, 37.1, 36.7, 36.1, 36.2, 45.1, 54.2, 57, 58.7, 36.1,
         40.1, 41.4, 42.6)
#########################################################################
#########################################################################

# scale inputs
df.max = apply(df, 2, max)
df.min = apply(df, 2, min)
Len = df.max - df.min
df.scale = df
for(i in 1:6){
  df.scale[ ,i] = (df[ ,i] - df.min[i])/Len[i]
}
df = df.scale


#########################################################################
#########################################################################
##### 2-level cokriging using T22 and T2exp 
set.seed(1234)
ind.training = sample(1:length(T2exp), 20)
yf = T2exp[ind.training]
Df = df[ind.training, ]

yc = T22

Dc = df 

input.testing = df[-ind.training, ]
y.testing = T2exp[-ind.training]

################################################################################
####### Matern_5_2 Correlation
################################################################################
### Gratiet's method
library(MuFiCokriging)
NestedDesign = NestedDesignBuild(design=list(Dc, Df))

bh.model = MuFicokm(formula=list(~1,~1),
                    MuFidesign=NestedDesign,
                    response=list(yc, yf),
                    nlevel=2, 
                    covtype="matern5_2")
summary(bh.model)$Cov.Val
bh.model$cok[[2]]@covariance@shape.val 


bh.pred = MuFiCokriging::predict(object=bh.model,
                                 newdata=input.testing,
                                 type="UK")

z.new = y.testing

#### numerical measures
sqrt(mean((z.new-bh.pred$mean)^2))
# 2.218553

bh.pred$lower95 = bh.pred$mean - 2* bh.pred$sig2^0.5
bh.pred$upper95 = bh.pred$mean + 2* bh.pred$sig2^0.5

mean((bh.pred$lower95<z.new & bh.pred$upper95>z.new))
#0.875
mean(bh.pred$upper95 - bh.pred$lower95)
# 7.315

####### ARCokrig
library(ARCokrig)
hyper = list(a=0.5-6, b=1, nugget.UB=1)
prior = list(name="JR", hyperparam=list(hyper, hyper, hyper))

prior = list(name="Reference")
prior = list(name="Jeffreys")

opt = list(maxit=2000)
obj = cokm(formula=list(~1,~1), output=list(yc, yf),
           input=list(Dc, Df), prior=prior, opt=opt,
           param=list(rep(0.1,6), rep(0.1,6)), 
           cov.model="matern_5_2")


## update model parameters in the cokm object
obj = cokm.fit(obj)

model.param = cokm.param(obj)
model.param$var
model.param$corr
model.param$coef

out = cokm.predict(obj, input.testing)

pred.ND = out

z.new = y.testing
level = 2

### RMSPE 
sqrt(mean((pred.ND$mu[[level]]-z.new)^2))
# 0.5236 #JR
# 0.5139 #Reference
# 0.5633 #Jeffreys


### coverage probability
mean((pred.ND$lower95[[2]]<z.new & pred.ND$upper95[[2]]>z.new))
# 0.875 #JR
# 0.875 #Reference
# 0.875 #Jeffreys


### ALCI
mean(pred.ND$upper95[[2]] - pred.ND$lower95[[2]])
# 2.153 #JR
# 2.3579 #Reference
# 1.9413 #Jeffreys


################################################################################
####### Power-Exponential Correlation
################################################################################
### Gratiet's method
library(MuFiCokriging)
NestedDesign = NestedDesignBuild(design=list(Dc, Df))

bh.model = MuFicokm(formula=list(~1,~1),
                    MuFidesign=NestedDesign,
                    response=list(yc, yf),
                    nlevel=2, 
                    covtype="powexp")
summary(bh.model)$Cov.Val
bh.model$cok[[2]]@covariance@shape.val 


bh.pred = MuFiCokriging::predict(object=bh.model,
                                 newdata=input.testing,
                                 type="UK")

z.new = y.testing

#### numerical measures
sqrt(mean((z.new-bh.pred$mean)^2))
# 1.9684

bh.pred$lower95 = bh.pred$mean - 2* bh.pred$sig2^0.5
bh.pred$upper95 = bh.pred$mean + 2* bh.pred$sig2^0.5

mean((bh.pred$lower95<z.new & bh.pred$upper95>z.new))
# 0.875 # powexp
mean(bh.pred$upper95 - bh.pred$lower95)
# 6.0280 # powexp

####### ARCokrig
library(ARCokrig)
hyper = list(a=0.5-6, b=1, nugget.UB=1)
prior = list(name="JR", hyperparam=list(hyper, hyper, hyper))
prior = list(name="Reference")
prior = list(name="Jeffreys")

opt = list(maxit=2000)
obj = cokm(formula=list(~1,~1), output=list(yc, yf),
           input=list(Dc, Df), prior=prior, opt=opt,
           param=list(rep(0.1,6), rep(0.1,6)), 
           cov.model="powexp")


## update model parameters in the cokm object
obj = cokm.fit(obj)
obj@param 

model.param = cokm.param(obj)
model.param$var
model.param$corr
model.param$coef

out = cokm.predict(obj, input.testing)

pred.ND = out

z.new = y.testing
level = 2

### RMSPE 
sqrt(mean((pred.ND$mu[[level]]-z.new)^2))
# 0.52953 # JR powexp
# 0.551407 # Reference powexp
# 0.542258 # Jeffreys powexp

### coverage probability
mean((pred.ND$lower95[[2]]<z.new & pred.ND$upper95[[2]]>z.new))
# 0.875 # JR powexp
# 0.875 # Reference powexp
# 0.875 # Jeffreys powexp 

### ALCI
mean(pred.ND$upper95[[2]] - pred.ND$lower95[[2]])
# 2.3085 # JR powexp
# 2.239676 # Reference powexp
# 2.236596 # Jeffreys powexp
