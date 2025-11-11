library(data.table)
library(FactoMineR)
library(sleev)

#d.mi.A<-read.csv("d.mi.A.csv")
d.mi.A <- fread("d.mi.A.csv")


## sleev doesn't like factor variables, so I'm changing all factor variables to dummy 0-1 variables
## sleev needs to have the same number of x and x_unval variables. egaWk was not collected
##  during phase-1, so I artificially create egaWk.ph1 and assign everyone to the value of 40.
d.mi.A$egaWk.ph1<-40
d.mi.A.num<-d.mi.A
d.mi.A.num$mat3race.black.ph1<-ifelse(d.mi.A.num$mat3race.ph1=="Black",1,0)
d.mi.A.num$mat3race.other.ph1<-ifelse(d.mi.A.num$mat3race.ph1=="Other",1,0)
d.mi.A.num$mat3race.black<-ifelse(d.mi.A.num$mat3race=="Black",1,0)
d.mi.A.num$mat3race.other<-ifelse(d.mi.A.num$mat3race=="Other",1,0)
d.mi.A.num$insurancer<-ifelse(d.mi.A.num$insurancer=="Public/Other",1,0)
d.mi.A.num$insurancer.ph1<-ifelse(d.mi.A.num$insurancer.ph1=="Public/Other",1,0)
d.mi.A.num$tobacPregr<-ifelse(d.mi.A.num$tobacPregr==1,1,0)
d.mi.A.num$tobacPregr.ph1<-ifelse(d.mi.A.num$tobacPregr.ph1==1,1,0)
d.mi.A.num$asthma.ph1<-as.numeric(d.mi.A.num$asthma.ph1)
d.mi.A.num$mat3race.ph1<-NULL
d.mi.A.num$mat3race<-NULL

## I have a lot of X* and Z variables. So I am collapsing them down to 2 variables
##  using the FactoMineR package (factor analysis). This of course is a simplification,
##  but a necessary one. FactoMineR seems to like factor variables, so I'm using the
##  original dataset with factor variables.

phase1.vars.A<-d.mi.A[,c("estWtChangePerWk.ph1","est_bmi_preg_mother.ph1","mat_age_delivery","mat3race.ph1",
                         "maleA","insurancer.ph1","tobacPregr.ph1","mat_asthma.ph1","egaWk.ph1")]
famd.A<-FAMD(base=phase1.vars.A)
famd.A$eig
famd.A$var
famd2.A<-famd.A$ind$coord[,c(1:2)]

d.mi.smle1.A<-d.mi.A.num
d.mi.smle1.A$famd1<-famd2.A[,1]
d.mi.smle1.A$famd2<-famd2.A[,2]

data.famd.A <- spline2ph(x=c("famd1","famd2"), size=4, degree=3, data=d.mi.smle1.A)

library(peakRAM)
start.time <- Sys.time()
system.time(peakRAM(
  mod.logistic.smle1 <- logistic2ph(y = "asthma",
                                  y_unval = "asthma.ph1",
                                  x = c("estWtChangePerWk","est_bmi_preg_mother",
                                        "mat3race.black","mat3race.other",
                                        "insurancer","tobacPregr","mat_asthma",
                                        "egaWk"),
                                  x_unval = c("estWtChangePerWk.ph1","est_bmi_preg_mother.ph1",
                                              "mat3race.black.ph1","mat3race.other.ph1",
                                              "insurancer.ph1","tobacPregr.ph1","mat_asthma.ph1"),
                                              #"egaWk.ph1"),
                                  z = c("mat_age_delivery","maleA"),
                                  data = data.famd.A, hn_scale = 1,
                                  se = TRUE, tol = 1e-04, max_iter = 1000,
                                  verbose = TRUE)
))
paste0("Run time: ", round(difftime(Sys.time(), start.time,
                                    units = "secs"), 3), " sec")

mod.logistic.smle1
