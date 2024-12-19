'With the code in this  file the thresholds for both data analysis and simulation study can be generated.
For convenience the outcome of the below code is already included in the github directory and can directly be used 
in the simulation or data analysis'

# All calculations were done under R version 3.4.4.

setwd(".")  #current directory
source("functions.R")

# Thresholds for real data analysis -----------------

## Thresholds 256 x 256 Pixels ---------------------------------------------

# Number of simulations run and number of pixels in one row
sim<-10000

# Generation of Null images
null.images<-noise(sim=sim,pixels=256,seed=1002306) # Seed 256 x 256

# Diameter of circle and angles of inner strip
d<-20
alpha<-list()
# Angles for maximization over 3, 6 and 9 angles, respectively
alpha[[1]]<-seq(0,120,60)
alpha[[2]]<-seq(0,150,30)
alpha[[3]]<-seq(0,160,20)

# Arrays for maxima over various combinations of width of inner strip, sidedness of statistics and angles
# Naming convention: stat.[name1]pixel.[name2].[name3]
#
# name1: Number of Pixels in one row of an image
# 
# name2: Indicating the width of the inner strip. 
# h1: Inner strip has width 1, h3: Inner strip has width 3,...
#
# name3: Indication whether one-sided or two-sided statistics are used. 1down=One-sided, 2side=two-sided

stat.256pixel.h1.1down<-array(NA,dim=c(sim,length(alpha),3))
stat.256pixel.h5.1down<-array(NA,dim=c(sim,length(alpha),3))
stat.256pixel.h1.2side<-array(NA,dim=c(sim,length(alpha),3))
stat.256pixel.h5.2side<-array(NA,dim=c(sim,length(alpha),3))

for(j in 1:length(alpha)){
  for(i in 1:sim){
    
    ### One-sided statistics
    
    # Statistic for FnB1 with h=1, maximizing over angles alpha[[j]]
    stat.h1.1down<-FnB(null.images[[i]],d.x = d,h=1,var.est = "none",alpha=alpha[[j]],side=-1,sig = 1)
    # Maximum of F1-statistic over image
    stat.256pixel.h1.1down[i,j,1]<-max(stat.h1.1down[[1]])
    # Maximum of nB-statistic over image
    stat.256pixel.h1.1down[i,j,2]<-max(stat.h1.1down[[2]])
    # Maximum of FnB1-statistic over image
    stat.256pixel.h1.1down[i,j,3]<-max(stat.h1.1down[[3]])
    # print(c(i,3*j,"h=1, 1down"))
    
    # Statistic for FnB1 with h=5, maximizing over angles alpha[[j]]
    stat.h5.1down<-FnB(null.images[[i]],d.x = d,h=5,var.est = "none",alpha=alpha[[j]],side=-1,sig = 1)
    # Maximum of F1-statistic over image
    stat.256pixel.h5.1down[i,j,1]<-max(stat.h5.1down[[1]])
    # Maximum of nB-statistic over image
    stat.256pixel.h5.1down[i,j,2]<-max(stat.h5.1down[[2]])
    # Maximum of FnB1-statistic over image
    stat.256pixel.h5.1down[i,j,3]<-max(stat.h5.1down[[3]])
    # print(c(i,3*j,"h=5, 1down"))
    
    
    ### Two-sided statistics
    
    # Statistic for FnB2 with h=1, maximizing over angles alpha[[j]]
    stat.h1.2side<-FnB(null.images[[i]],d.x = d,h=1,var.est = "none",alpha=alpha[[j]],side=0,sig = 1)
    # Maximum of F2-statistic over image
    stat.256pixel.h1.2side[i,j,1]<-max(stat.h1.2side[[1]])
    # Maximum of nB-statistic over image
    stat.256pixel.h1.2side[i,j,2]<-max(stat.h1.2side[[2]])
    # Maximum of FnB2-statistic over image
    stat.256pixel.h1.2side[i,j,3]<-max(stat.h1.2side[[3]])
    # print(c(i,3*j,"h=1, 2side"))
    
    # Statistic for FnB2 with h=5, maximizing over angles alpha[[j]]
    stat.h5.2side<-FnB(null.images[[i]],d.x = d,h=5,var.est = "none",alpha=alpha[[j]],side=0,sig = 1)
    # Maximum of F2-statistic over image
    stat.256pixel.h5.2side[i,j,1]<-max(stat.h5.2side[[1]])
    # Maximum of nB-statistic over image
    stat.256pixel.h5.2side[i,j,2]<-max(stat.h5.2side[[2]])
    # Maximum of FnB2-statistic over image
    stat.256pixel.h5.2side[i,j,3]<-max(stat.h5.2side[[3]])
    # print(c(i,3*j,"h=5, 2side"))
  }
}


### Calculation of Thresholds
#
# Naming convention: thr.[name1].[name2].[name3]
# 
# name1: Name of statistics
# F1: One-sided F-statistic, F2: Two-sided F-statistic,...
#
# name2: Number of Pixels in one row of an image
# 
# name3: Indicating the width of the inner strip. 
# h1: Inner strip has width 1, h3: Inner strip has width 3,...


thr.F1.256.h1<-apply(stat.256pixel.h1.1down[,,1],2,quantile,prob=0.95)
thr.nB.256.h1<-apply(stat.256pixel.h1.1down[,,2],2,quantile,prob=0.95)
thr.FnB1.256.h1<-apply(stat.256pixel.h1.1down[,,3],2,quantile,prob=0.95)

thr.F1.256.h5<-apply(stat.256pixel.h5.1down[,,1],2,quantile,prob=0.95)
thr.nB.256.h5<-apply(stat.256pixel.h5.1down[,,2],2,quantile,prob=0.95)
thr.FnB1.256.h5<-apply(stat.256pixel.h5.1down[,,3],2,quantile,prob=0.95)

thr.F2.256.h1<-apply(stat.256pixel.h1.2side[,,1],2,quantile,prob=0.95)
thr.FnB2.256.h1<-apply(stat.256pixel.h1.2side[,,3],2,quantile,prob=0.95)

thr.F2.256.h5<-apply(stat.256pixel.h5.2side[,,1],2,quantile,prob=0.95)
thr.FnB2.256.h5<-apply(stat.256pixel.h5.2side[,,3],2,quantile,prob=0.95)

# Pooling of thresholds in one array for each width of the inner strip

thr.h1<-array(NA,dim=c(3,3,2))
thr.h5<-array(NA,dim=c(3,3,2))

thr.h1[1,,1]<-thr.F1.256.h1
thr.h1[2,,1]<-thr.nB.256.h1
thr.h1[3,,1]<-thr.FnB1.256.h1

thr.h1[1,,2]<-thr.F2.256.h1
thr.h1[2,,2]<-thr.nB.256.h1
thr.h1[3,,2]<-thr.FnB2.256.h1

thr.h5[1,,1]<-thr.F1.256.h5
thr.h5[2,,1]<-thr.nB.256.h5
thr.h5[3,,1]<-thr.FnB1.256.h5

thr.h5[1,,2]<-thr.F2.256.h5
thr.h5[2,,2]<-thr.nB.256.h5
thr.h5[3,,2]<-thr.FnB2.256.h5



## Thresholds 306 x 306 Pixels ---------------------------------------------

# Number of simulations run and number of pixels in one row
sim<-1

# Generation of Null images
null.images<-noise(sim=sim,pixels=256,seed=1002306) # Seed 306 x 306

# Arrays for maxima over various combinations of width of inner strip, sidedness of statistics and angles
# Naming convention: stat.[name1]pixel.[name2].[name3]
#
# name1: Number of Pixels in one row of an image
# 
# name2: Indicating the width of the inner strip. 
# h1: Inner strip has width 1, h3: Inner strip has width 3,...
#
# name3: Indication whether one-sided or two-sided statistics are used. 1down=One-sided, 2side=two-sided

stat.306pixel.h8.1down<-array(NA,dim=c(sim,length(alpha),3))
stat.306pixel.h8.2side<-array(NA,dim=c(sim,length(alpha),3))

for(j in 1:length(alpha)){
  for(i in 1:sim){
    ### One-sided statistics
    
    # Statistic for FnB1 with h=8, maximizing over angles alpha[[j]]
    stat.h8.1down<-FnB(null.images[[i]],d.x = d,h=8,var.est = "none",alpha=alpha[[j]],side=-1,sig = 1)
    # Maximum of F1-statistic over image
    stat.306pixel.h8.1down[i,j,1]<-max(stat.h8.1down[[1]])
    # Maximum of nB-statistic over image
    stat.306pixel.h8.1down[i,j,2]<-max(stat.h8.1down[[2]])
    # Maximum of FnB1-statistic over image
    stat.306pixel.h8.1down[i,j,3]<-max(stat.h8.1down[[3]])
    # print(c(i,3*j,"h=1, 1down"))
    
    ### Two-sided statistics
    
    # Statistic for FnB1 with h=8, maximizing over angles alpha[[j]]
    stat.h8.2side<-FnB(null.images[[i]],d.x = d,h=8,var.est = "none",alpha=alpha[[j]],side=0,sig = 1)
    # Maximum of F1-statistic over image
    stat.306pixel.h8.2side[i,j,1]<-max(stat.h8.2side[[1]])
    # Maximum of nB-statistic over image
    stat.306pixel.h8.2side[i,j,2]<-max(stat.h8.2side[[2]])
    # Maximum of FnB1-statistic over image
    stat.306pixel.h8.2side[i,j,3]<-max(stat.h8.2side[[3]])
    # print(c(i,3*j,"h=1, 1down"))
  }
}


### Calculation of Thresholds
#
# Naming convention: thr.[name1].[name2].[name3]
# 
# name1: Name of statistics
# F1: One-sided F-statistic, F2: Two-sided F-statistic,...
#
# name2: Number of Pixels in one row of an image
# 
# name3: Indicating the width of the inner strip. 
# h1: Inner strip has width 1, h3: Inner strip has width 3,...


thr.F1.306.h8<-apply(stat.306pixel.h8.1down[,,1],2,quantile,prob=0.95)
thr.nB.306.h8<-apply(stat.306pixel.h8.1down[,,2],2,quantile,prob=0.95)
thr.FnB1.306.h8<-apply(stat.306pixel.h8.1down[,,3],2,quantile,prob=0.95)

thr.F2.306.h8<-apply(stat.306pixel.h8.2side[,,1],2,quantile,prob=0.95)
thr.FnB2.306.h8<-apply(stat.306pixel.h8.2side[,,3],2,quantile,prob=0.95)

# Pooling of thresholds in one array

thr.h8<-array(NA,dim=c(3,3,2))

thr.h8[1,,1]<-thr.F1.306.h8
thr.h8[2,,1]<-thr.nB.306.h8
thr.h8[3,,1]<-thr.FnB1.306.h8

thr.h8[1,,2]<-thr.F2.306.h8
thr.h8[2,,2]<-thr.nB.306.h8
thr.h8[3,,2]<-thr.FnB2.306.h8


save(thr.h1,thr.h5,thr.h8,file="Thresholds_real_data.RData")


# Thresholds for simulation study -------------------------


# Number of simulations run and number of pixels in one row
sim<-10000
pixels<-100

# Generate images with only noise
null.images<-noise(sim=sim,pixels=pixels,48299281)

# Diameter of circle, width and angle of inner strip
d<-10
h.in<-2
alpha<-0

max.null.images<-rep(NA,sim)

# Calculate maximum of one-sided FnB-statistic for each image with true variance
for(i in 1:sim){
  max.null.images[i]<-max(FnB(null.images[[i]],d.x = d,h=h.in,var.est = "none",alpha=alpha,side=-1,sig = 1)[[3]])
  print(i)
}

# Calculate Threshold
thr<-quantile(max.null.images,prob=0.95)
save(thr,file="Threshold_simulation_study.RData")

