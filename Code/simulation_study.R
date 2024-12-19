# All calculations were done under R version 3.4.4.

setwd(".") # set working directory as current directory

source("functions.R")

load("Threshold_simulation_study.RData")  

'If preferred, the thresholds can be generated using the file Generate_Thresholds.R'

# Empirical Power  -------------------------------------------------------------------


# Number of simulations run and number of pixels in one row
sim<-1000
pixels<-100

# Noise of the images
noise.image<-noise(sim=sim,pixels=pixels,seed=74928121)

# Length of rectangles
len<-30
# Widths of rectangles
width<-1:3

# x- and y-Coordinates of rectangle anchor point
wherex.rect<-35
wherey.rect<-60
# Angle in degrees by which the rectangles are turned
alpha.fissure<-40

# Diameter of circle, width and angle misspecification of inner strips
d<-10
h<-2
delta<-seq(0,25,5)
# Angles of inner strips
alpha.strip<-40-delta
# Signal of gray values of rectangles
signal<--seq(1,3,0.5)

# Array for maximum of statistic of pixels adjacent to fissure
max.stat.rectangle<-array(NA,dim=c(sim,length(delta),length(width),length(signal)))


canvas<-matrix(0,pixels,pixels)

for(w in width){
  # Define fissure
  type<-rectangle(canvas,1,wherex.rect,wherey.rect,len,wid=w,alpha = alpha.fissure)
  # Adjacent pixels to fissure
  adjacent<-adjacent.pixels(pixels=pixels,type=type,d=d)
  wch.adj<-which(adjacent==1,arr.ind = T)
  
  # Noisy images with fissures
  for(s in 1:length(signal)){
    fields<-list()
    for(i in 1:sim){
      fields[[i]]<-noise.image[[i]]+type*signal[s]
    }
    
    for(j in 1:length(alpha.strip)){
      for(i in 1:sim){
        # Calculate FnB1-statistic and Maximum of all pixels adjacent to fissure
        stat<-FnB(fields[[i]],d.x = d,h=h,var.est = "robust_global",alpha=alpha.strip[j],side=-1)[[3]]
        max.stat.rectangle[i,j,w,s]<-max(stat[wch.adj])
        print(c(i,delta[j],w,signal[s]))
      }
    }
    
  }
}
# save(max.stat.rectangle,file="Maxima_Simulation_Study_Adjacent_Pixels.RData")

# Array for detection rates
detection.rates<-array(NA,dim=c(length(delta),length(width),length(signal)))

for(w in width){
  for(s in 1:length(signal)){
    # Determine whether a maximum is above the threshold (hence the fissure detected)
    detected<-max.stat.rectangle[,,w,s]>thr
    # Calculate detection rates for fixed combination of width w and signal s
    detection.rates[,w,s]<-apply(detected,2,mean)
  }
}

# Generate graphics to display detection rates
# Colors
farbe<-c("blue","green3","orange","red","firebrick")
# point characters used
pch<-c(16,17,18,15,8)
# Sizes of graphic, labels for axes
cx<-2
cx.leg<-2
cx.axis<-3
xlb<-""
ylb<-""


###
### Disclosure: The labeling in the paper ("Detection rates", "w=0.01", "w=0.02", "w=0.03", "\Delta") has been added via PDF-Annotator, Version 9.0.0.916
###

# Plotting of detection rates
windows(32,12)
par(mfrow=c(1,3))
for(w in 1:3){
  plot(0,type="n",xlim=c(0,25),ylim=c(0,1),xlab=xlb,ylab=ylb,main="",xaxt="n",yaxt="n")
  axis(1,delta,delta,cex.axis=cx.axis,mgp=c(3,3,0))
  axis(2,c(0,0.5,1),c(0,0.5,1),cex.axis=cx.axis,mgp=c(3,1,0))
  for(s in 1:length(signal)){
    points(delta,detection.rates[,w,s],col=farbe[s],cex=cx,pch=pch[s])
    
    for(j in 1:5){
      lines(delta[c(j,j+1)],detection.rates[,w,s][c(j,j+1)],col=farbe[s])
    }
  }
}

# Save image
dev.print(device=pdf,file="../Graphics/output_images/Detection_rates_1angle.pdf")
