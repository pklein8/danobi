###
###
### The raw images being analyzed are "1Fissure_width1.png", "1Fissure_width5.png" and "Steel_fiber.png". 
###
### The images "1Fissure_width1.png" and "1Fissure_width5.png" are 2D slices of semi-synthetic 3D images published in
###
### Barisin, T., Jung, C., Muesebeck, F., Redenbach, C., and Schladitz, K. (2022a). 
### "Methods for segmenting cracks in 3d images of concrete: A comparison based on semi-synthetic images." In: Pattern Recognition 129
###
### where artificial cracks were added to 3D images of uncracked concrete by Franziska Muesebeck (TU Kaiserslautern). 
### This is open data under the CC BY license http://creativecommons.org/licenses/by/4.0/.
###
### The image "Steel_fiber.png" is one of the concrete samples investigated in 
### 
### Schuler, F. (2006). "Richtungsanalyse von Fasern in Beton und Charakterisierung von rissquerenden Fasern mittels Computer-Tomografie." Dissertation. TU Kaiserslautern.
###
### imaged by Franz Schreiber at ITWM.
###
###

# All calculations were done under R version 4.1.2.

setwd(".")              # set working directory to current path
source("functions.R")
 


 load("Thresholds_real_data.RData")                # Load pre-calculated thresholds

 'If preferred, the thresholds can be generated using the file Generate_Thresholds.R'

# Graphics generation -----------------------------------------------------


# install.packages("png")
# install.packages("latex2exp")
# install.packages("fields")
# install.packages("colorspace")

# Package to load png-files
library(png)
# Package to use function image.plot()
library(fields)
# Package to define customized color palettes
library(colorspace)


###
### Define color palettes for images of statistic and significant pixels
###

load("Palettes.RData")          # load pre-defined palettes

### Manual definition of palette for statistics
#
# Type of palette: Basic: Sequential (single-hue) -> leftmost palette
#
# H1=0
# C1=0
# L1=0
# L2=97
# P1=1.3
# Correct all colors to valid RGB color model values: Yes
# n=50
# Reverse colors

# pal.basis1<-choose_palette()
# pal.stat<-pal.basis1(100)


### Manual definition of palette for significant pixels
#
# Type of palette: Advanced: Sequential (single-hue) -> middle palette
#
# H1=10
# C1=115
# CMAX=150
# L1=50
# L2=90
# P1=1.1
# P2=1.3
# Correct all colors to valid RGB color model values: Yes
# n=5
# Reverse colors

# pal.basis2<-choose_palette()
# pal.signif<-pal.basis2(2)

###
### Disclosure: All graphic titles were added via PDF-Annotator, Version 9.0.0.916
###

# Graphics for image containing crack of width 5 --------------------------

# Loading the image

data<-readPNG("../Graphics/raw_images/1Fissure_width5.png")

# Diameter of circle, width and angles of inner strip
d<-20
h.in<-5
alpha<-list()
alpha[[1]]<-seq(0,120,60)
alpha[[2]]<-seq(0,150,30)
alpha[[3]]<-seq(0,160,20)
side<-c(-1,0)

# Variables for naming of graphics

side.name<-c("1down","2side")
alpha.name<-c("3Angles","6Angles","9Angles")
stat.name<-matrix(c("F1","nB","FnB1","F2","nB","FnB2"),2,3,byrow = T)

# Naming convention
#
# File name: 1F_w5_h5_[name1]_[name2]Angles_[statname]_[side].pdf
#
# 1F_w5:    As graphic contains one artificial fissure of width 5
# h5:       As inner strip used for scan has width 5
# name1:    From {stat,signif}.
#           name1=stat:   (Rescaled) statistic is plotted
#           name2=signif: Significant pixels are displayed
# name2:    Number of angles over which the statistic is maximized
# statname: From {F1,F2,nB,FnB1,FnB2}: Name of the statistic displayed
# side:     From {1down,2side,1up}
#           side=1down: One-sided F- and FnB-statistics are calculated scanning for anomalies with lower expected gray values
#           side=1up: One-sided F- and FnB-statistics are calculated scanning for anomalies with higher expected gray values
#           side=2side: Two-sided F- and FnB-statistics are calculated


windows(28,32)

for(k in 1:length(side)){
  for(j in 1:length(alpha)){
    # Calculate F-, nB- and FnB-statistic
    stat<-FnB(data,d.x = d,h=h.in,var.est = "robust_global",alpha=alpha[[j]],side=side[k])
    for(i in 1:3){
      # Normalize statistics by threshold, cap values at 2*threshold
      stat.resc<-matrix(pmin(2,stat[[i]]/thr.h5[i,j,k]),256,256)
      # Determine significant pixels
      stat.signif<-stat.resc>1
      
      # Plot normalized statistic and save image
      image(stat.resc)
      image.plot(stat.resc,zlim=c(0,2),add=T,smallplot = c(.90, .92, .2, .8),col=pal.stat,axis.args=list( at=c(0,2), labels=c("low","high")))
      dev.print(device = pdf,file=paste("../Graphics/output_images/1F_w5_h",h.in,"_stat_",alpha.name[j],"_",stat.name[k,i],"_",side.name[k],".pdf",sep=""))
      # Plot significant pixels and save image
      image(stat.signif,col=pal.signif)
      dev.print(device = pdf,file=paste("../Graphics/output_images/1F_w5_h",h.in,"_signif_",alpha.name[j],"_",stat.name[k,i],"_",side.name[k],".pdf",sep=""))
    }
  }
}

### For the sake of better visibility, the scales of the following two graphics have been adjusted

# Calculate F-, nB- and FnB-statistic
stat<-FnB(data,d.x = d,h=h.in,var.est = "robust_global",alpha=seq(0,160,20),side=0)

# Save F2- and nB-statistic
stat.F2<-stat[[1]]
stat.nB<-stat[[2]]

# Plot statistics and save images
image(stat.F2)
image.plot(stat.F2,zlim=c(0,max(stat.F2)),add=T,smallplot = c(.90, .92, .2, .8),col=pal.stat,axis.args=list( at=c(0,max(stat.F2)), labels=c("low","high")))
dev.print(device=pdf,file="../Graphics/output_images/1F_w5_h5_stat_9Angles_F2_2side.pdf")

image(stat.nB)
image.plot(stat.nB,zlim=c(0,max(stat.nB)),add=T,smallplot = c(.90, .92, .2, .8),col=pal.stat,axis.args=list( at=c(0,max(stat.nB)), labels=c("low","high")))
dev.print(device=pdf,file="../Graphics/output_images/1F_w5_h5_stat_9Angles_nB_2side.pdf")

# Graphics for image containing crack of width 1 --------------------------

# setwd(".")

# Loading the image
data<-readPNG("../Graphics/raw_images/1Fissure_width1.png")

# Diameter of circle, width and angles of inner strip
d<-20
h.in<-1
alpha<-list()
alpha[[1]]<-seq(0,120,60)
alpha[[2]]<-seq(0,150,30)
alpha[[3]]<-seq(0,160,20)
side<-c(-1,0)

# Variables for naming of graphics

side.name<-c("1down","2side")
alpha.name<-c("3Angles","6Angles","9Angles")
stat.name<-matrix(c("F1","nB","FnB1","F2","nB","FnB2"),2,3,byrow = T)

# Naming convention
#
# File name: 1F_w1_h1_[name1]_[name2]Angles_[statname]_[side].pdf
#
# 1F_w1:    As graphic contains one artificial fissure of width1
# h1:       As inner strip used for scan has width 1
# name1:    From {stat,signif}.
#           name1=stat:   (Rescaled) statistic is plotted
#           name2=signif: Significant pixels are displayed
# name2:    Number of angles over which the statistic is maximized
# statname: From {F1,F2,nB,FnB1,FnB2}: Name of the statistic displayed
# side:     From {1down,2side,1up}
#           side=1down: One-sided F- and FnB-statistics are calculated scanning for anomalies with lower expected gray values
#           side=1up: One-sided F- and FnB-statistics are calculated scanning for anomalies with higher expected gray values
#           side=2side: Two-sided F- and FnB-statistics are calculated


###
### Disclosure: As the significant pixels are barely visible, we have enlarged them in the paper by adding larger squares
### centered around the spots of the significant pixels with PDF-Annotator, Version 9.0.0.916.
###



windows(28,32)

for(k in 1:length(side)){
  for(j in 1:length(alpha)){
    # Calculate F-, nB- and FnB-statistic
    stat<-FnB(data,d.x = d,h=h.in,var.est = "robust_global",alpha=alpha[[j]],side=side[k])
    for(i in 1:3){
      # Normalize statistics by threshold, cap values at 2*threshold
      stat.resc<-matrix(pmin(2,stat[[i]]/thr.h1[i,j,k]),256,256)
      # Determine significant pixels
      stat.signif<-stat.resc>1
      
      # Plot normalized statistic and save image
      image(stat.resc)
      image.plot(stat.resc,zlim=c(0,2),add=T,smallplot = c(.90, .92, .2, .8),col=pal.stat,axis.args=list( at=c(0,2), labels=c("low","high")))
      dev.print(device = pdf,file=paste("../Graphics/output_images/1F_w1_h",h.in,"_stat_",alpha.name[j],"_",stat.name[k,i],"_",side.name[k],".pdf",sep=""))
      # Plot significant pixels and save image
      image(stat.signif,col=pal.signif)
      dev.print(device = pdf,file=paste("../Graphics/output_images/1F_w1_h",h.in,"_signif_",alpha.name[j],"_",stat.name[k,i],"_",side.name[k],".pdf",sep=""))
    }
  }
}





# Graphics for image containing steelfibers --------------------------

# setwd(".") #current directory

# Loading the image
data<-t(readPNG("../Graphics/raw_images/Steel_fiber.png")[,,1])[,306:1]

# Diameter of circle, width and angles of inner strip
d<-20
h.in<-8
alpha<-list()
alpha[[1]]<-seq(0,120,60)
alpha[[2]]<-seq(0,150,30)
alpha[[3]]<-seq(0,160,20)
side<-c(-1,0)

# Variables for naming of graphics

side.name<-c("1down","2side")
alpha.name<-c("3Angles","6Angles","9Angles")
stat.name<-matrix(c("F1","nB","FnB1","F2","nB","FnB2"),2,3,byrow = T)

# Naming convention
#
# File name: Steelfiber_h8_[name1]_[name2]Angles_[statname]_[side].pdf
#
# Steelfiber:    As graphic contains steelfibers
# h8:            As inner strip used for scan has width 8
# name1:         From {stat,signif}.
#                name1=stat:   (Rescaled) statistic is plotted
#                name2=signif: Significant pixels are displayed
# name2:         Number of angles over which the statistic is maximized
# statname:      From {F1,F2,nB,FnB1,FnB2}: Name of the statistic displayed
# side:          From {1down,2side,1up}
#                side=1down: One-sided F- and FnB-statistics are calculated scanning for anomalies with lower expected gray values
#                side=1up: One-sided F- and FnB-statistics are calculated scanning for anomalies with higher expected gray values
#                side=2side: Two-sided F- and FnB-statistics are calculated


windows(28,32)

for(k in 1:length(side)){
  for(j in 1:length(alpha)){
    # Calculate F-, nB- and FnB-statistic
    stat<-FnB(data,d.x = d,h=h.in,var.est = "robust_global",alpha=alpha[[j]],side=side[k])
    for(i in 1:3){
      # Normalize statistics by threshold, cap values at 2*threshold
      stat.resc<-matrix(pmin(2,stat[[i]]/thr.h8[i,j,k]),306,306)
      # Determine significant pixels
      stat.signif<-stat.resc>1
      
      # Plot normalized statistic and save image
      image(stat.resc)
      image.plot(stat.resc,zlim=c(0,2),add=T,smallplot = c(.90, .92, .2, .8),col=pal.stat,axis.args=list( at=c(0,2), labels=c("low","high")))
      dev.print(device = pdf,file=paste("../Graphics/output_images/Steelfiber_h",h.in,"_stat_",alpha.name[j],"_",stat.name[k,i],"_",side.name[k],".pdf",sep=""))
      # Plot significant pixels and save image
      image(stat.signif,col=pal.signif)
      dev.print(device = pdf,file=paste("../Graphics/output_images/Steelfiber_h",h.in,"_signif_",alpha.name[j],"_",stat.name[k,i],"_",side.name[k],".pdf",sep=""))
    }
  }
}

