rm(list=ls())
setwd("~/Calcium-transient-analysis/")
inputfile="wt-2.xlsx"

# setting cutoff for detect peaks
# According to your data and requirement to set the following values:

MinPeakHeight =25 
MinPeakWidth=10
MinPeakDistance=50
#======================================================
#install.packages("TrenchR")
library("TrenchR")
#install.packages("gsignal")
library("gsignal")
#Functions==============================================
cal_ang=function(point_1, point_2, point_3){
  #  calculate angle
  #  :param point_1: 
  #  :param point_2: 
  #  :param point_3: 
  #  :return: angle
  a=sqrt((point_2[1]-point_3[1])*(point_2[1]-point_3[1])+(point_2[2]-point_3[2])*(point_2[2] - point_3[2]))
  b=sqrt((point_1[1]-point_3[1])*(point_1[1]-point_3[1])+(point_1[2]-point_3[2])*(point_1[2] - point_3[2]))
  c=sqrt((point_1[1]-point_2[1])*(point_1[1]-point_2[1])+(point_1[2]-point_2[2])*(point_1[2]-point_2[2]))
  B=radians_to_degrees(acos((b*b-a*a-c*c)/(-2*a*c))) 
  return(B)
}

#Read Data and smooth Data==============================
{
  dd=xlsx::read.xlsx(inputfile,1)
  x=dd$x.s.
  y=dd$Y.intensity.
  plot(x,y,type="p",cex=.3) 
  dat=data.frame(x,y)
  data_df=data.frame(id=1:nrow(dat),x=dat$x,y=dat$y)
  library(stats)
  # smooth data
  smoothed_y <- stats::filter(data_df$y, rep(1/10, 10), sides = 2)  
  data_df$yy=as.numeric(smoothed_y)
  data_df=na.omit(data_df)
  plot(data_df$x,data_df$y,type="p",cex=.3,xlab = "Time", ylab = "Intensity")
  lines(data_df$x,data_df$yy,col="red",type="l") 
  dat=data_df
}

#Find peaks==========================
{
  peaks=gsignal::findpeaks(data_df$yy,MinPeakHeight = MinPeakHeight,MinPeakDistance=MinPeakDistance,MinPeakWidth=MinPeakWidth)
  #points(data_df$x[peaks$loc], peaks$pks, col = "red", pch = 1)
  data_df$point_type=""
  data_df$point_type[peaks$loc]="summit_point"
  data_df$peakID=""
  data_df[peaks$loc,"peakID"]=paste("peak",seq(1,length(peaks$loc)),sep="")
  data_df$peakID[is.na(data_df$peakID)]=""
  
  data_df1=data_df[!is.na(data_df$yy),]
  data_df1$point_type[1]="start_point"
  data_df1$point_type[nrow(data_df1)]="end_point"
  peak_points=data_df1[data_df1$point_type!="",]
}

#Find left and right botton points======
{
  for(i in 2:(nrow(peak_points)-1)){
    point1=peak_points[i-1,]
    point2=peak_points[i,]
    point3=peak_points[i+1,]
    mypoints=rbind(point1,point2,point3)
    points(mypoints[,c(2,3)],col="red",cex=1)
    #to find left start point and median of peak
    {
      temp_points=data_df[point1$id<data_df$id & data_df$id<point2$id,]
      middle_point1_point2=temp_points[temp_points$yy==min(temp_points$yy),][1,]
      if(middle_point1_point2$id>10){
        middle_point1_point2=temp_points[temp_points$id==(middle_point1_point2$id-10),]
      }
      angle1=360
      k=0
      temp_points$angle=NA
      for(ii in (middle_point1_point2$id+1):max(temp_points$id)){
        moving_point=temp_points[temp_points$id==ii,]
        angle = cal_ang(middle_point1_point2[,c(1,4)],moving_point[,c(1,4)],point2[,c(1,4)])
        angle=angle$id
        temp_points$angle[temp_points$id==ii]=angle
        if(!is.nan(angle1) && angle<angle1){
          angle1=angle
          k=ii
        }
      }
      data_df$point_type[data_df$id==k]="left_start_point"
      data_df$peakID[data_df$id==k]=point2$peakID
      points(data_df[data_df$id==k,c(2,3)],col="blue")
      
      temp_points=data_df[k<data_df$id & data_df$id<point2$id,]
      temp_points$dis2mean=abs(temp_points$y-mean(c(max(temp_points$y),min(temp_points$y))))
      median_point=temp_points[temp_points$dis2mean==min(temp_points$dis2mean),][1,]
      data_df$point_type[data_df$id==median_point$id]="left_middle_point"
      data_df$peakID[data_df$id==median_point$id]=point2$peakID
      median_point_height=data_df$y[data_df$id==median_point$id]
    }    
    
    #to find right end point and median of peak
    {
      temp_points=data_df[point2$id<data_df$id & data_df$id<point3$id,]
      if(nrow(temp_points)==0){next}
      middle_point2_point3=temp_points[temp_points$y==min(temp_points$y),][1,]
      if(middle_point2_point3$y>median_point_height){next}
      temp_points=data_df[point2$id<data_df$id & data_df$id<middle_point2_point3$id,]
      temp_points$dis2mean=abs(temp_points$y-median_point_height)
      median_point=temp_points[temp_points$dis2mean==min(temp_points$dis2mean),][1,]
      data_df$point_type[data_df$id==median_point$id]="right_middle_point"
      data_df$peakID[data_df$id==median_point$id]=point2$peakID
    } 
  }
}

key_points=data_df[data_df$point_type!="",]
points(key_points[grep("left_start_point",key_points$point_type),c(2,3)],col="blue",cex=1)
points(key_points[grep("middle_point",key_points$point_type),c(2,3)],col="green",cex=1)
#peak2peak distance==============
{
  peaks=key_points[key_points$point_type=="summit_point",]
  peaks$Peak_peak=NA
  peaks$Peak2Peak_distance=NA
  i=2
  while(i<=nrow(peaks)){
    peaks$Peak_peak[i]=paste(as.character(peaks$peakID[i-1]),"-",as.character(peaks$peakID[i]),sep="")
    peaks$Peak2Peak_distance[i]=peaks$x[i]-peaks$x[i-1]
    i=i+1
  }
  peaks=na.omit(peaks[,c(7,8)])
  write.table(peaks,file=paste(inputfile,"peak2peak_distance.txt",sep="_"),quote=F,row.names = F,sep="\t")
}

#Duration_time======================
{
  peak_middle_width=key_points[key_points$point_type %in% c("left_middle_point","right_middle_point"),]
  peak_middle_width$Peak_middle_width=NA
  peak_middle_width=peak_middle_width[order(peak_middle_width$x),]
  i=2
  while(i<=nrow(peak_middle_width)){
    peak_middle_width$Peak_middle_width[i]=peak_middle_width$x[i]-peak_middle_width$x[i-1]
    i=i+2
  }
  peak_middle_width=na.omit(peak_middle_width[,c(6,7)])
  write.table(peak_middle_width,file=paste(inputfile,"Duration_time.txt",sep="_"),quote=F,row.names = F,sep="\t")
}


#Rise_Time======================
{
  rise_time=key_points[key_points$point_type %in% c("left_start_point","summit_point"),]
  rise_time$rise_time=NA
  rise_time=rise_time[order(rise_time$x),]
  i=2
  while(i<=nrow(rise_time)){
    rise_time$rise_time[i]=rise_time$x[i]-rise_time$x[i-1]
    i=i+2
  }
  rise_time=rise_time[,c(6,7)]
  rise_time=na.omit(rise_time)
  write.table(rise_time,file=paste(inputfile,"rise_time.txt",sep="_"),quote=F,row.names = F,sep="\t")
}

#decay time======================
{
  decay_time=key_points[key_points$point_type %in% c("summit_point","left_start_point"),]
  decay_time$decay_time=NA
  decay_time=decay_time[order(decay_time$x),]
  decay_time=decay_time[2:nrow(decay_time),]
  decay_time$type=NA
  i=2
  while(i<=nrow(decay_time)){
    decay_time$decay_time[i]=decay_time$x[i]-decay_time$x[i-1]
    decay_time$type[i]=paste(decay_time$peakID[i-1],"_",decay_time$point_type[i-1],"--to--",decay_time$peakID[i],"_",decay_time$point_type[i],sep="")
    i=i+2
  }
  decay_time=decay_time[,c(8,7)]
  decay_time=na.omit(decay_time)
  write.table(decay_time,file=paste(inputfile,"decay_time.txt",sep="_"),quote=F,row.names = F,sep="\t")
}

#Delta_F======================
{
  deltaF=key_points[key_points$point_type %in% c("left_start_point","summit_point"),]
  deltaF$deltaF=NA
  deltaF=deltaF[order(deltaF$x),]
  i=2
  while(i<=nrow(deltaF)){
    deltaF$deltaF[i]=(deltaF$y[i]-deltaF$y[i-1])/deltaF$y[i-1]
    i=i+2
  }
  deltaF=na.omit(deltaF[,c(2,3,6,7)])
  write.table(decay_time,file=paste(inputfile,"Delta_F.txt",sep="_"),quote=F,row.names = F,sep="\t")
}

