
# Required functions for calculation of FnB statistic ---------------------



###
### FnB: Function to calculate the FnB-statistic
###
### Input:
###
### data     - square matrix 
### ball     - Indicator whether a circle is used as scan window 
###            If ball=F, a rectangle is used as scan window
### d.x, d.y - If ball=T: Radius of circle (must coincide)
###            If ball=F: length and width of rectangle, respectively
### h        - Width of "Inner strip"
### var.est  - Method for variance estimation
###            Possible input: c("none","global","robust_global","minimum","robust_local")
### alpha    - Angles of inner strip given in degree
### rescale  - Indicator whether scan statistics with sample means or rescaled sample means according to asymptotic is used
###            If rescale=F: Sample means are used
###            If rescale=T: Sample means are multiplied by t (number of pixels in one dimension)
### sig      - Standard deviation if var.est="none"
### side     - Which version of FnB-statistic is used
###            side=1: One-sided FnB-statistic searching for fissures with higher expected gray values
###            side=0: Two-sided FnB-statistic
###            side=-1: One-sided FnB-statistic searching for fissures with lower expected gray values
###
### Output: List of 3
###         1st entry: Either F1- or F2-statistic (F2 when side=0, else F1)
###         2nd entry: nB-statistic
###         3rd entry: FnB-statistc (FnB2 when side=0, else FnB1)


FnB<-function(data,d.x,d.y=d.x,h,var.est="robust_global",alpha,ball=T,rescale=F,sig=1,side=0){
  
  # Abortion criteria for which the function shall stop
  if(!is.matrix(data) | nrow(data)!=ncol(data)){stop("data must be given as square matrix.")}
  if(!is.numeric(d.x) | !is.numeric(d.y) | !is.numeric(h) |d.x<=0 | d.y<=0 | h<=0){stop("d.x, d.y and h must be given as positive numeric values")}
  if(any(c(d.x,d.y,h)<1) & any(c(d.x,d.y,h)>=1)){stop("Bandwidths must either all be given as proportion of T or as positive values >= 1")}
  if(ball & d.x!=d.y){stop("If ball=T, d.x must be equal to d.y")}
  if(all(var.est!=c("none","global","robust_global","minimum","robust_local"))){stop("var.est must be from c('none','global','robust_global','minimum','robust_local')")}
  if(!is.numeric(alpha)){stop("alpha must be numeric")}
  if(all(ball!=c(T,F))){stop("ball must be either TRUE or FALSE.")}
  if(all(rescale!=c(T,F))){stop("rescale must be either TRUE or FALSE.")}
  if(var.est=="none" & (!is.numeric(sig) | sig<=0)){stop("For var.est='none', sig must be positive and numeric.")}
  if(all(side!=c(-1,0,1))){stop("side must be either -1, 0 or 1.")}
  
  # Number of pixels in one row
  t<-nrow(data)
  
  # Potential rescaling of d.x, d.y and h
  if(all(c(d.x,d.y,h)<1)){
    d.x<-t*d.x
    d.y<-t*d.y
    h<-t*h
  }
  
  # Conversion of alpha from degree to radians
  alpha<- alpha %% 180
  alpha<-sort(alpha-180*(alpha>90))
  alpha<-alpha*pi/180
  
  # Creation of list with pixels belonging to scan sets A_1,...,A_5 for all angles
  pixel.scansets.list<-list()
  
  for(i in 1:length(alpha)){
    pixel.scansets.list[[i]]<-pixel.scansets(d.x,d.y,h,alpha[i],ball)
  }
  
  fissure.stat<-matrix(0,t,t) # F-statistic
  bubble.stat<-matrix(0,t,t)  # nB-statistic
  combi<-matrix(0,t,t)        # FnB-statistic
  

    if(var.est=="none" |var.est=="global"|var.est=="robust_global"){
      
      # Calculation of Standard deviation
      if(var.est=="none"){sigma<-sig}
      if(var.est=="global"){sigma<-sd(data)}
      if(var.est=="robust_global"){sigma<-(quantile(data,3/4)-quantile(data,1/4))/(2*qnorm(3/4))}
      
      # Calculation of F- and nB-statistic for each pixel
      for(i in (floor(d.x/2)+1):ceiling(t-d.x/2)){
        for(j in (floor(d.y/2)+1):ceiling(t-d.y/2)){
          
          # Choose matrix for scan window with anchoring point (i,j)
          data.local<-data[ceiling(i-d.x/2):floor(i+d.x/2),ceiling(j-d.y/2):floor(j+d.y/2)]
          
          # Calculation of FnB-statistic for all angles
          for(k in 1:length(alpha)){
            
            # Pixels beloning to scan windows A_1(i,j),...,A_5(i,j)
            bound.A1<-pixel.scansets.list[[k]][[1]]
            bound.A2<-pixel.scansets.list[[k]][[2]]
            bound.A3<-pixel.scansets.list[[k]][[3]]
            bound.A4<-pixel.scansets.list[[k]][[4]]
            bound.A5<-pixel.scansets.list[[k]][[5]]
            
            # Calculation of sample means in scan windows
            mean.A1<-mean(data.local[bound.A1])
            mean.A2<-mean(data.local[bound.A2])
            mean.A3<-mean(data.local[bound.A3])
            mean.A4<-mean(data.local[bound.A4])
            mean.A5<-mean(data.local[bound.A5])
            
            # Calculation of contrast between A_1 and A_2, as well as between A_1 and A_3
            contrast12<-(mean.A1-mean.A2)*side+abs(mean.A1-mean.A2)*(1-abs(side))
            contrast13<-(mean.A1-mean.A3)*side+abs(mean.A1-mean.A3)*(1-abs(side))
            
            m1<-min(contrast12,contrast13)
            m2<-abs(mean.A4-mean.A5)
            
            # Maximization over angles
            fissure.stat[i,j]<-max(fissure.stat[i,j],m1)
            bubble.stat[i,j]<-max(bubble.stat[i,j],m2)
          }
        }
      }
      
      # Potential rescaling of statistics
      fissure.stat<-fissure.stat/sigma*t^(rescale)
      bubble.stat<-bubble.stat/sigma*t^(rescale)
      combi<-pmax(fissure.stat-bubble.stat,0)
    }
    
    if(var.est=="minimum"){
      
      # Calculation of F- and nB-statistic for each pixel
      for(i in (floor(d.x/2)+1):ceiling(t-d.x/2)){
        for(j in (floor(d.y/2)+1):ceiling(t-d.y/2)){
          
          # Choose matrix for scan window with anchoring point (i,j)
          data.local<-data[ceiling(i-d.x/2):floor(i+d.x/2),ceiling(j-d.y/2):floor(j+d.y/2)]
          
          # Calculation of FnB-statistic for all angles
          for(k in 1:length(alpha)){
            
            # Pixels beloning to scan windows A_1(i,j),...,A_5(i,j)
            bound.A1<-pixel.scansets.list[[k]][[1]]
            bound.A2<-pixel.scansets.list[[k]][[2]]
            bound.A3<-pixel.scansets.list[[k]][[3]]
            bound.A4<-pixel.scansets.list[[k]][[4]]
            bound.A5<-pixel.scansets.list[[k]][[5]]
            
            # Calculation of sample means in scan windows
            mean.A1<-mean(data.local[bound.A1])
            mean.A2<-mean(data.local[bound.A2])
            mean.A3<-mean(data.local[bound.A3])
            mean.A4<-mean(data.local[bound.A4])
            mean.A5<-mean(data.local[bound.A5])
            
            # Calculation of contrast between A_1 and A_2, as well as between A_1 and A_3
            contrast12<-(mean.A1-mean.A2)*side+abs(mean.A1-mean.A2)*(1-abs(side))
            contrast13<-(mean.A1-mean.A3)*side+abs(mean.A1-mean.A3)*(1-abs(side))
            
            m1<-min(contrast12,contrast13)
            m2<-abs(mean.A4-mean.A5)
            
            # Minima of empirical standard deviations in A_2 and A_3, as well as A_4 and A_5 
            s1<-min(sd(data.local[bound.A2]),sd(data.local[bound.A3]))
            s2<-min(sd(data.local[bound.A4]),sd(data.local[bound.A5]))
            
            # Maximization over angles
            fissure.stat[i,j]<-max(fissure.stat[i,j],m1/s1)
            bubble.stat[i,j]<-max(bubble.stat[i,j],m2/s2)
          }
        }
      }
      # Potential rescaling of statistics
      fissure.stat<-fissure.stat*t^(rescale)
      bubble.stat<-bubble.stat*t^(rescale)
      combi<-pmax(fissure.stat-bubble.stat,0)
    }
    
    if(var.est=="robust_local"){
      
      # Pixels belonging to whole circle
      bound.outer<-pixel.scansets.list[[1]][[7]]
      
      # Calculation of F- and nB-statistic for each pixel
      for(i in (floor(d.x/2)+1):ceiling(t-d.x/2)){
        for(j in (floor(d.y/2)+1):ceiling(t-d.y/2)){
          data.local<-data[ceiling(i-d.x/2):floor(i+d.x/2),ceiling(j-d.y/2):floor(j+d.y/2)]
          s2<-(quantile(data.local[bound.outer],3/4)-quantile(data.local[bound.outer],1/4))/(2*qnorm(3/4))
          for(k in 1:length(alpha)){
            # Pixels beloning to scan windows A_1(i,j),...,A_5(i,j) and A_2(i,j)\cup A_3(i,j)
            bound.A1<-pixel.scansets.list[[k]][[1]]
            bound.A2<-pixel.scansets.list[[k]][[2]]
            bound.A3<-pixel.scansets.list[[k]][[3]]
            bound.A4<-pixel.scansets.list[[k]][[4]]
            bound.A5<-pixel.scansets.list[[k]][[5]]
            bound.A6<-pixel.scansets.list[[k]][[6]]
            
            # Calculation of sample means in scan windows
            mean.A1<-mean(data.local[bound.A1])
            mean.A2<-mean(data.local[bound.A2])
            mean.A3<-mean(data.local[bound.A3])
            mean.A4<-mean(data.local[bound.A4])
            mean.A5<-mean(data.local[bound.A5])
            
            # Calculation of contrast between A_1 and A_2, as well as between A_1 and A_3
            contrast12<-(mean.A1-mean.A2)*side+abs(mean.A1-mean.A2)*(1-abs(side))
            contrast13<-(mean.A1-mean.A3)*side+abs(mean.A1-mean.A3)*(1-abs(side))
            
            m1<-min(contrast12,contrast13)
            m2<-abs(mean.A4-mean.A5)
            
            # Local robust estimator of standard deviation based on empirical quantiles in A_2(i,j)\cup A_3(i,j)
            s1<-(quantile(data.local[bound.A6],3/4)-quantile(data.local[bound.A6],1/4))/(2*qnorm(3/4))
            
            # Maximization over angles
            fissure.stat[i,j]<-max(fissure.stat[i,j],m1/s1)
            bubble.stat[i,j]<-max(bubble.stat[i,j],m2/s2)
          }
        }
      }
      
      # Potential rescaling of statistics
      fissure.stat<-fissure.stat*t^(rescale)
      bubble.stat<-bubble.stat*t^(rescale)
      combi<-pmax(fissure.stat-bubble.stat,0)
    }

  
  return(list(fissure.stat,bubble.stat,combi))
}


###
### pixel.scansets: Function to calculate pixels belonging to scan sets
###
### Input:
###
### ball     - Indicator whether a circle is used as scan window 
###            If ball=F, a rectangle is used as scan window
### d.x, d.y - If ball=T: Radius of circle (must coincide)
###            If ball=F: length and width of rectangle, respectively
### h        - Width of "Inner strip"
### alpha    - Angle of inner strip given in radians
###
### Output: List of 7
###         1st-5th entry: Pixels belonging to A_1,...,A_5
###         6th entry: Pixels belonging to union of A_2 and A_3
###         7th entry: Pixels belonging to whole circle (if ball=T) or rectangle (if ball=F)


pixel.scansets<-function(d.x,d.y,h,alpha,ball){

  # Anchor point of scan sets and upper bound of rectangle
  p.x<-floor(d.x/2)+1
  oben.x<-floor(d.x/2+floor(d.x/2))+1
  
  p.y<-floor(d.y/2)+1
  oben.y<-floor(d.y/2+floor(d.y/2))+1

  ind<-matrix(rep(1:oben.y,oben.x),oben.x,oben.y,byrow=T)
  
  x<-1:oben.x
  
  # Calculate boundary of circle
  if(ball==T){
    grenze.unten<-p.x-sqrt(d.x^2/4-(p.x-x)^2)
    grenze.oben<-p.x+sqrt(d.x^2/4-(p.x-x)^2)
  }else{
    grenze.unten<-rep(-Inf,length(x))
    grenze.oben<-rep(Inf,length(x))
  }

  
  # Lines to separate A_1 and A_2, A_1 and A_3 as well as A_4 and A_5
  y1<-tan(alpha)*x+p.y-p.x*tan(alpha)-h/2*(cos(alpha)+sin(alpha)*tan(alpha))  
  y2<-tan(alpha)*x+p.y-p.x*tan(alpha)+h/2*(cos(alpha)+sin(alpha)*tan(alpha))
  y3<-tan(alpha)*x+p.y-p.x*tan(alpha)
  
  # Pixels belonging to A_1,...,A_5
  # A6: Union of A_2 and A_3
  # A7: Either whole circle or whole rectangle
  bound.A1<-which((ind<=grenze.oben & ind>=grenze.unten & ind>=y1 & ind<=y2),arr.ind=T)
  bound.A2<-which((ind<=grenze.oben & ind>=grenze.unten & ind<y1),arr.ind=T)
  bound.A3<-which((ind<=grenze.oben & ind>=grenze.unten & ind>y2),arr.ind=T)
  bound.A6<-which((ind<=grenze.oben & ind>=grenze.unten & (ind>y2| ind<y1)),arr.ind=T)
  
  bound.A4<-which((ind<=grenze.oben & ind>=grenze.unten & ind>y3),arr.ind=T)
  bound.A5<-which((ind<=grenze.oben & ind>=grenze.unten & ind<=y3),arr.ind=T)
  bound.A7<-which((ind<=grenze.oben & ind>=grenze.unten),arr.ind=T)

  return(list(bound.A1,bound.A2,bound.A3,bound.A4,bound.A5,bound.A6,bound.A7))
}


# Artificial generation of 'fissures' and 'bubbles' -----------------------


###
### bubble, rectangle: Functions to add a 'bubble' (=circle) and 'fissure' (=rectangle) to an image
###
### Input:
###
### data           - square matrix
### value          - signal strength of circle/rectangle
### wherex, wherey - x- and y-components of anchoring point of circle/rectangle
### radius         - radius of circle
### alpha          - Angle of rectangle turned against x-axis in degree
### len, wid       - length and width of rectangle
###
### Output: image with added bubble/fissure
###

bubble<-function(data,value,wherex,wherey,radius){
  if(!is.matrix(data)){stop("data must be given as matrix")}
 
   # Lower and upper boundaries of y-components of circle
  grenze.unten<-wherey-floor(sqrt(radius^2-(seq(-radius,radius,1))^2))
  grenze.oben<-wherey+floor(sqrt(radius^2-(seq(-radius,radius,1))^2))
  
  for(i in 1:length(seq(-radius,radius,1))){
    # Add signal
    data[floor(wherex-radius+i-1),(grenze.unten[i]):(grenze.oben[i])]<-value+data[floor(wherex-radius+i-1),(grenze.unten[i]):(grenze.oben[i])]
  }
  
  return(data)
}

rectangle<-function(data,value,alpha,wherex,wherey,len,wid){
  l.data<-nrow(data)
  
  for(i in 1:l.data){
    for(j in 1:l.data){
      # Add signal
      data[i,j]<-data[i,j]+is.rectangle(c(i,j),alpha,wherex,wherey,wid,len)*value
    }
  }
  return(data)
}

###
### is.rectangle: Function to decide whether a point p belongs to a rectangle
###
### Input: p -  point to be tested to lie in a given rectangle
###        all other inputs: as in rectangle()
###
### Output: TRUE/FALSE whether p belongs to rectangle
###

is.rectangle<-function(p,alpha,wherex,wherey,wid,len){
  alpha<-(alpha*pi/180)%%(2*pi)
  
  if(alpha==0){
    rect<-((p[1]-wherex>=0) & (p[1]-wherex<=len) & (p[2]-wherey>=0) & (p[2]-wherey<=wid))
  }
  
  if(alpha==pi/2){
    rect<-((p[1]-wherex<=0) & (p[1]-wherex>=-wid) & (p[2]-wherey>=0) & (p[2]-wherey<=len))
  }
  
  if(alpha==pi){
    rect<-((p[1]-wherex<=0) & (p[1]-wherex>=-len) & (p[2]-wherey<=0) & (p[2]-wherey>=-wid))
    
  }
  
  if(alpha==3*pi/2){
    rect<-((p[1]-wherex>=0) & (p[1]-wherex<=wid) & (p[2]-wherey<=0) & (p[2]-wherey>=-len))
  }
  
  
  if(alpha>0 & alpha<pi/2){
    rect<-((sign(g1(p,alpha,wherex,wherey,wid,len))<=0) & (sign(g2(p,alpha,wherex,wherey,wid,len))<=0) & (sign(g3(p,alpha,wherex,wherey,wid,len))>=0) & (sign(g4(p,alpha,wherex,wherey,wid,len))>=0))
  }
  
  if(alpha>pi/2 & alpha<pi){
    rect<-((sign(g1(p,alpha,wherex,wherey,wid,len))>=0) & (sign(g2(p,alpha,wherex,wherey,wid,len))<=0) & (sign(g3(p,alpha,wherex,wherey,wid,len))<=0) & (sign(g4(p,alpha,wherex,wherey,wid,len))>=0))
  }
  
  if(alpha>3/2*pi & alpha<pi){
    rect<-((sign(g1(p,alpha,wherex,wherey,wid,len))<=0) & (sign(g2(p,alpha,wherex,wherey,wid,len))>=0) & (sign(g3(p,alpha,wherex,wherey,wid,len))>=0) & (sign(g4(p,alpha,wherex,wherey,wid,len))<=0))
  }
  
  if(alpha>pi & alpha<3/2*pi){
    rect<-((sign(g1(p,alpha,wherex,wherey,wid,len))>=0) & (sign(g2(p,alpha,wherex,wherey,wid,len))>=0) & (sign(g3(p,alpha,wherex,wherey,wid,len))<=0) & (sign(g4(p,alpha,wherex,wherey,wid,len))<=0))
  }
  
  return(rect)
  
}

###
### g1-g4 Parametrization of lines defining the boundary of a rectangle
###

g1<-function(p,alpha,wherex,wherey,wid,len){
  tan(alpha)*p[1]+wherey-wherex*tan(alpha)-p[2]
}

g2<-function(p,alpha,wherex,wherey,wid,len){
  -1/tan(alpha)*p[1]+wherey+wherex/tan(alpha)-p[2]
}

g3<-function(p,alpha,wherex,wherey,wid,len){
  tan(alpha)*p[1]+wherey-wherex*tan(alpha)+wid*(cos(alpha)+tan(alpha)*sin(alpha))-p[2]
}

g4<-function(p,alpha,wherex,wherey,wid,len){
  -1/tan(alpha)*p[1]++wherey+wherex/tan(alpha)+len*(sin(alpha)+cos(alpha)/tan(alpha))-p[2]
}

###
### ball.location: Function to determine all pixels being 'adjacent' to an anomaly, i.e.
###                intersection between scan window and anomaly is non-empty
###
### Input: 
###
### t    - Number of pixels in one row
### type - square matrix with entries of 0 and 1 designating whether pixels belong to anomaly
### h    - radius of scan window
###
### Output: Matrix with entries 0 and 1 displaying whether pixels are adjacent to anomaly
###


adjacent.pixels<-function(pixels,type,d){
  
  if(!(is.matrix(type) & nrow(type)==ncol(type))){stop("type must be given as square matrix")}
  if(!all(type==1 | type==0)){stop("type must contain only entries of 0 or 1")}
  
  part<-matrix(0,pixels,pixels)
  
  # Pixels belonging to anomaly
  wch.fiss<-which(type==1,arr.ind=T)
  
  # Determine for all pixels (k,m) whether they are adjacent to anomaly
  for(m in (ceiling(d/2)):(floor(pixels-d/2))){
    for(k in (ceiling(d/2)):(ceiling(pixels-d/2))){
      
      
      # Create circle of diameter d around (k,m)
      l<-floor(k+d/2)-floor(k-d/2)
      ms<-m-floor(m-d/2)
      
      ind<-matrix(rep(1:l,l),l,l,byrow=T)
      grenze.unten<-floor(ms-sqrt(ms^2-(ms-1:l)^2))+1
      grenze.oben<-floor(ms+sqrt(ms^2-(ms-1:l)^2))
      wch<-which((ind>=grenze.unten & ind<=grenze.oben),arr.ind=T)
      add<-matrix(NA,length(wch[,1]),2)
      add[,1]<-floor(k-d/2)
      add[,2]<-floor(m-d/2)
      wch<-wch+add
      
      # Comparison of pixels in circle around (k,m) with pixels in anomaly
      vgl<-0
      for(i in 1:length(wch[,1])){
        for(j in 1:length(wch.fiss[,1])){
          vgl<-vgl+all(wch[i,]==wch.fiss[j,])
        }
      }
      part[k,m]<-(vgl>=1)
    }
  }
  
  return(part)
}



# Generation of artificial images -----------------------------------------


noise<-function(sim,pixels,seed){
  # Set seed
  set.seed(seed)
  
  noise<-list()
  
  for(i in 1:sim){
    # Create [pixels] x [pixels]-matrix with i.i.d. N(0,1)-distributed entries
    noise[[i]]<-matrix(rnorm(pixels^2),pixels,pixels)
  }
  
  return(noise)
}

