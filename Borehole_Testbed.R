
rm(list=ls())

############################################################################
## Install the ARCokrig Package
## install relevant packages if not installed 
install.packages(c("Rcpp", "RcppEigen", "RcppArmadillo", "RobustGaSP", "mvtnorm", "fields"))
## install ARCokrig package from local source
install.packages("ARCokrig_0.1.0.tar.gz", repos=NULL, type="source")
## one can also install the package from github
library(devtools)
#install_github("pulongma/ARCokrig")
library(ARCokrig)


library(ggplot2)
library(reshape2)
library(plyr)






############################################################################
############################################################################
############# The Borehole Functions  

borehole <- function(xx)
{
  ##########################################################################
  #
  # BOREHOLE FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # OUTPUT AND INPUT:
  #
  # y  = water flow rate
  # xx = c(rw, r, Tu, Hu, Tl, Hl, L, Kw)
  #
  ##########################################################################
  
  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]
  
  frac1 <- 2 * pi * Tu * (Hu-Hl)
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}

boreholelc <- function(xx)
{ 
  ##########################################################################
  #
  # BOREHOLE FUNCTION, LOWER FIDELITY CODE
  # This function is used as the "low-accuracy code" version of the function
  # borehole.r.
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # OUTPUT AND INPUT:
  #
  # y = water flow rate
  # xx = c(rw, r, Tu, Hu, Tl, Hl, L, Kw)
  #
  #########################################################################
  
  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]
  
  frac1 <- 5 * Tu * (Hu-Hl)
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1.5+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}
############################################################################
############################################################################



############################################################################
############################################################################
########## The Borehole Example 

library(DiceDesign)
n=100
df.lf = lhsDesign(n=n,dimension=8,seed=(12345))
input = df.lf$design
## randomly held-out 20% for cross-validation
nc = n*0.8
set.seed(123)
indDc = sample(1:n, nc)
input.lf = input[indDc, ]

## convert toe original scale
XX.range = matrix(c(0.05, 0.15, 
                    100, 50000,
                    63070, 115600,
                    990, 1110,
                    63.1, 116,
                    700, 820,
                    1120, 1680,
                    9855, 12045), byrow=T, ncol=2)

Lens = (XX.range[ ,2]-XX.range[ ,1])

## create inputs for prediction
input.new = input[-indDc, ]
input.new = t(XX.range[ ,1] + t(input.new %*% diag(XX.range[ ,2]-XX.range[ ,1])))

# get true value from high-fidelity code
z.new = apply(input.new, 1, borehole)
z.new.hf = apply(input.new, 1, borehole)
z.new.lf = apply(input.new, 1, boreholelc)

## convert inputs to original scale
XX.lf = t(XX.range[ ,1] + t(input.lf %*% diag(XX.range[ ,2]-XX.range[ ,1])))

################################################################################
######## Matern-5/2 Correlation
################################################################################
###### Gratiet(2013)
library(MuFiCokriging)
nf = 30
set.seed(345)
DNest = NestedDesign(XX.lf, nlevel=2, n=nf)
XX.hf = ExtractNestDesign(DNest, 2)

## obtain data
zc = apply(XX.lf, 1, boreholelc)
zf = apply(XX.hf, 1, borehole)

bh.model = MuFicokm(formula=list(~1,~1),
                    MuFidesign=DNest,
                    response=list(zc,zf),
                    nlevel=2,
                    covtype="matern5_2")
summary(bh.model)$Cov.Val

bh.pred = MuFiCokriging::predict(object=bh.model,
                                 newdata=input.new,
                                 type="UK")

sqrt(mean((bh.pred$mean-z.new)^2))
# 1.936


bh.pred$lower95 = bh.pred$mean - 2* bh.pred$sig2^0.5
bh.pred$upper95 = bh.pred$mean + 2* bh.pred$sig2^0.5
mean((bh.pred$lower95<z.new & bh.pred$upper95>z.new))
# 1
mean(bh.pred$upper95 - bh.pred$lower95)
# 17.85


#########################################################################
#########################################################################
###### ARCokrig
## create the cokm object 
hyper = list(a=0.2, b=1, nugget.UB=1)
prior = list(name="JR", hyperparam=list(hyper, hyper, hyper))

#prior = list(name="Reference")
#prior = list(name="Jeffreys")
obj = cokm(formula=list(~1,~1), output=list(c(zc), c(zf)),
           input=list(XX.lf, XX.hf), prior=prior,
           param=list(Lens/2, Lens/2), cov.model="matern_5_2")

## update model parameters in the cokm object
obj = cokm.fit(obj)


## prediction 
out = cokm.predict(obj, input.new)

pred.ND = out

# RMSPE 
sqrt(mean((pred.ND$mu[[2]]-z.new)^2))
#0.464 #Reference
#0.466 #Jeffreys
#0.379 #JR

mean((pred.ND$lower95[[2]]<z.new & pred.ND$upper95[[2]]>z.new))
#0.85 #Reference
#0.85 #Jeffreys
#0.95 #JR

mean(pred.ND$upper95[[2]] - pred.ND$lower95[[2]])
#1.353 #Reference
#1.359 #Jeffreys
#1.436 #JR


#######################################################################
#######################################################################
######### figures


## compare ARCokrig versus Gratiet (2013)
df = data.frame(x=z.new.hf, y1=pred.ND$mu[[2]], y2=bh.pred$mux[[2]])
df.melt = melt(df, id="x")

df.melt$variable = revalue(df.melt$variable, c("y1"="New Formulas", "y2"="LG"))

g = ggplot(data=df.melt, aes(x=x, y=value, color=variable)) + 
  geom_point(shape=1, color="black") + 
  facet_wrap(variable ~ ., ncol=2) + 
  geom_abline(intercept=0, slope=1) + 
  theme(plot.title=element_text(size=14),
        axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=14),
        legend.text = element_text(size=12),
        #legend.direction = "horizontal",
        legend.position = c(0.6, 0.1)) + 
  xlab("model output") + ylab("predicted output")

print(g)

############################################################################
############################################################################


################################################################################
######## Power-Exponential Correlation
################################################################################
###### Gratiet(2013)
library(MuFiCokriging)
nf = 30
set.seed(345)
DNest = NestedDesign(XX.lf, nlevel=2, n=nf)
XX.hf = ExtractNestDesign(DNest, 2)

## obtain data
zc = apply(XX.lf, 1, boreholelc)
zf = apply(XX.hf, 1, borehole)

bh.model = MuFicokm(formula=list(~1,~1),
                    MuFidesign=DNest,
                    response=list(zc,zf),
                    nlevel=2,
                    covtype="powexp")
summary(bh.model)$Cov.Val

bh.pred = MuFiCokriging::predict(object=bh.model,
                                 newdata=input.new,
                                 type="UK")
bh.model$cok[[1]]@covariance@shape.val 
bh.model$cok[[2]]@covariance@shape.val

sqrt(mean((bh.pred$mean-z.new)^2))
#2.5873 # powexp

bh.pred$lower95 = bh.pred$mean - 2* bh.pred$sig2^0.5
bh.pred$upper95 = bh.pred$mean + 2* bh.pred$sig2^0.5
mean((bh.pred$lower95<z.new & bh.pred$upper95>z.new))
# 1 #powexp
mean(bh.pred$upper95 - bh.pred$lower95)
# 19.6826 # powexp

#########################################################################
#########################################################################
###### ARCokrig
## create the cokm object
hyper = list(a=0.2, b=1, nugget.UB=1)
prior = list(name="JR", hyperparam=list(hyper, hyper, hyper))
prior = list(name="Reference")
prior = list(name="Jeffreys")

obj = cokm(formula=list(~1,~1), output=list(c(zc), c(zf)),
           input=list(XX.lf, XX.hf), prior=prior,
           param=list(Lens/2, Lens/2), cov.model="powexp")

## update model parameters in the cokm object
obj = cokm.fit(obj)


## prediction 
out = cokm.predict(obj, input.new)

pred.ND = out

# RMSPE 
sqrt(mean((pred.ND$mu[[2]]-z.new)^2))
# 1  # JR with powexp
# 0.7994 # Reference with powexp
# 0.8378 # Jeffreys with powexp

mean((pred.ND$lower95[[2]]<z.new & pred.ND$upper95[[2]]>z.new))
# 1 # JR with powexp
# 1 # Reference with powexp
# 1 # Jeffreys with powexp

mean(pred.ND$upper95[[2]] - pred.ND$lower95[[2]])
# 11.61103 # JR with powexp
# 6.028 # Reference with powexp
# 6.1787 # Jeffreys with powexp


#######################################################################
#######################################################################
######### figures


## compare ARCokrig versus Gratiet (2013)
df = data.frame(x=z.new.hf, y1=pred.ND$mu[[2]], y2=bh.pred$mux[[2]])
df.melt = melt(df, id="x")

df.melt$variable = revalue(df.melt$variable, c("y1"="New Formulas", "y2"="LG"))

g = ggplot(data=df.melt, aes(x=x, y=value, color=variable)) + 
  geom_point(shape=1, color="black") + 
  facet_wrap(variable ~ ., ncol=2) + 
  geom_abline(intercept=0, slope=1) + 
  theme(plot.title=element_text(size=14),
        axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=14),
        legend.text = element_text(size=12),
        #legend.direction = "horizontal",
        legend.position = c(0.6, 0.1)) + 
  xlab("model output") + ylab("predicted output")

print(g)

