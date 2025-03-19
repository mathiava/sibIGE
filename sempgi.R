# SEM-PGI

library(lavaan)

# full sibling sem-pgi model

dat1 = dat[!is.na(dat$gm) & !is.na(dat$gp) & !is.na(dat$g1) & !is.na(dat$g2) & 
               !(is.na(dat$x1) & is.na(dat$x2)),]
   
model1 <-  "c =~ 1*x1 + 1*x2
   
   g1 + g2 ~ 0.5*gm + 0.5*gp
   
   x1 ~ g*g1 + ig*g2
   x2 ~ g*g2 + ig*g1
   
   x1 + x2 ~ fm*gm + fp*gp
   
   gm ~~ v_gm*gm
   gp ~~ v_gp*gp
   gp ~~ r*gm
   g1 ~~ 0.5*g1
   g2 ~~ 0.5*g2
   x1 ~~ e1*x1
   x2 ~~ e2*x2
   
   c ~~ vc*c
   
   gm + gp + g1 + g2 + x1 + x2 ~ 1
"

# run model and inspect 

result1 <- sem(model1, data=dat1, fixed.x=F, missing="FIML")
summary(result1, fit.measures=T, standardized=T)

# reduced sem-pgi model without sibling indirect genetic effects

model0 <- "
   
   c =~ 1*x1 + 1*x2
   
   g1 + g2 ~ 0.5*gm + 0.5*gp
   
   x1 ~ g*g1
   x2 ~ g*g2
   
   x1 + x2 ~ fm*gm + fp*gp
   
   gm ~~ v_gm*gm
   gp ~~ v_gp*gp
   gp ~~ r*gm
   g1 ~~ 0.5*g1
   g2 ~~ 0.5*g2
   x1 ~~ e1*x1
   x2 ~~ e2*x2
   
   c ~~ vc*c
   
   gm + gp + g1 + g2 + x1 + x2 ~ 1
"

result0 <- sem(model0, data=dats1, fixed.x=F, missing="FIML")
summary(result0, fit.measures=T, standardized=T)

# test model 0 and model 1 loglikelihood difference

-2*(as.numeric(result0@loglik$loglik)-as.numeric(result1@loglik$loglik))
pchisq(-2*(as.numeric(result0@loglik$loglik)-as.numeric(result1@loglik$loglik)), as.numeric(fitMeasures(result0, "df"))-as.numeric(fitMeasures(result1, "df")), lower.tail=FALSE)

# get variance due to different parameters in the full model

s = inspect(result1, what="est")$beta[3,4]
g = inspect(result1, what="est")$beta[2,4]
fm = inspect(result1, what="est")$beta[2,6]
fp = inspect(result1, what="est")$beta[2,7]
m = inspect(result1, what="est")$psi[6,7]
vm = inspect(result1, what="est")$psi[6,6]
vp = inspect(result1, what="est")$psi[7,7]
c = inspect(result1, what="est")$psi[1,1]
e1 = inspect(result1, what="est")$psi[2,2]

(g**2)*((2+2*m+vm+vp)/4)
(s**2)*((2+2*m+vm+vp)/4)
(fm**2)*vm
(fp**2)*vp
