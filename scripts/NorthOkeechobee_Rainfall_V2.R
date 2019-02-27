## Kissimmee River Rainfall
## 
## Code was compiled by Paul Julian
## contact infor: paul.julian@dep.state.fl.us 

#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

#Libraries
library(plyr);
library(reshape);
library(RColorBrewer)
library(RcppRoll)
library(mblm)

#GIS libraries
library(maptools)
library(classInt)
library(GISTools)
library(rgdal)
library(sp)
library(tmap)
library(tmaptools)
library(raster)
library(spatstat)
library(sf)
library(HatchedPolygons)
library(spatialEco)

#install.packages("tlocoh", repos="http://R-Forge.R-project.org")
library(tlocoh)

#Custom Functions
#source("Y:/CommonlyUsedFunctions.r")
source("D:/CommonlyUsedFunctions.r")

thessian_create.v2=function(data,clip,plot=T,full.extent=T){
  require(dismo)
  require(sp)
  require(raster)
  require(rgeos)
  if(full.extent==T){bbox.da=c(bbox(clip)[1,1],bbox(clip)[1,2],bbox(clip)[2,1],bbox(clip)[2,2])}else{bbox.da=NULL}
  th=voronoi(data,ext=bbox.da)
  th.z=sp::over(th,data)
  th.z.spdf=sp::SpatialPolygonsDataFrame(th,th.z)
  th.clp=raster::intersect(clip,th.z.spdf)
  th.clp$area=rgeos::gArea(th.clp,byid=T)
  if(plot==T){plot(th.clp)}
  return(th.clp)
}

WWI.index.fun=function(RF){
  require(zoo)
  dif.val=c(NA,diff(RF))
  tot.val=c(NA,rollapply(RF,width=2,FUN=sum,na.rm=T))
  WWI=dif.val/tot.val
  return(WWI)
}

display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE,colorblindFriendly=T)
#Paths
#setwd("//Fldep1/OWPER/EVG/Programmatic/WaterQuality/2_WQ Project Analyses/Basin_RainfallLoad")
setwd("D:/UF/WeatherWhip")

paths=paste(getwd(),c("/Data/","/export/","/plots/","/GIS/"),sep="")
#Folder.Maker(paths);#One and done. Creates folders in working directory.
data.path=paths[1]
export.path=paths[2]
plot.path=paths[3]
gis.path=paths[4]

##
nad83.pro=CRS("+init=epsg:4269")
state.plane=CRS("+proj=tmerc +lat_0=24.33333333333333 +lon_0=-81 +k=0.999941177 +x_0=200000.0001016002 +y_0=0 +ellps=GRS80 +to_meter=0.3048006096012192 +no_defs")
SFWMD.proj=CRS("+proj=tmerc +lat_0=24.33333333333333 +lon_0=-81 +k=0.9999411764705882 +x_0=199999.9999999999 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=us-ft +no_defs")

#GIS Data
basins.all=readOGR(paste0(gis.path,"SHP"),"SFWMD_WATERSHED") 
n.ok.basins=data.frame(NAME=c(c("ALLIGATOR LAKE", "BOGGY CREEK", "CATFISH CREEK", 
                               "EAST LAKE TOHOPEKALIGA", "ECONLOCKHATCHEE RIVER SWAMP", "HORSE CREEK (CLOSED BASIN)", 
                               "LAKE CONLIN", "LAKE CYPRESS", "LAKE GENTRY", "LAKE HART", "LAKE HATCHINEHA", 
                               "LAKE JACKSON", "LAKE KISSIMMEE", "LAKE MARIAN", "LAKE MARION", 
                               "LAKE MYRTLE", "LAKE PIERCE", "LAKE ROSALIE", "LAKE TOHOPEKALIGA", 
                               "LAKE WEOHYAKAPKA", "LOWER REEDY CREEK", "MARION CREEK", "S63A", 
                               "SHINGLE CREEK", "TIGER LAKE", "UPPER REEDY CREEK"),c("S-65BC","S-65D","S-65E","S-65A"),
                             c("L-59E","L-60W","L-60E","C-41AN","L-59W","C-41AS","C-40","C-41S","L-49","L-48","C-41N","S-131"),
                             c("FISHEATING CREEK/L-61"),c("S-154C","S-154","S-133","S191","TAYLOR CREEK STA","S-135","LAKESIDE RANCH STA","NUBBIN SLOUGH STA")),
                      Basin_gen=c(rep("Upper Kissimmee",26),rep("Lower Kissimmee",4),rep("HarneyPond_IndianPrairie",12),c("Fisheating Creek"),rep("NubbinSlough_TaylorCreek",8)))
n.ok=merge(basins.all,n.ok.basins,"NAME")
n.ok=subset(n.ok,is.na(Basin_gen)==F)

n.ok.dis=unionSpatialPolygons(n.ok,n.ok@data$Basin_gen)
n.ok.dis=SpatialPolygonsDataFrame(spTransform(n.ok.dis,SFWMD.proj),data.frame(Basin_gen=names(n.ok.dis),row.names=names(n.ok.dis),stringsAsFactors = F))
plot(subset(n.ok.dis,Basin_gen=="Upper Kissimmee"))
#plot(n.ok.dis)
#n.ok.dis=clean_slivers(n.ok.dis,min_nodes=9)$sp
#plot(n.ok.dis)
#plot(subset(n.ok.dis,Basin_gen=="Lower Kissimmee"))

active.hydromet=read.csv(paste0(gis.path,"Active Hydrometeorologic Stations.csv"))
active.wx=subset(active.hydromet,ACTIVITY_SUBTYPE%in%c("Rain"))
active.wx=SpatialPointsDataFrame(coords=active.wx[,c("LON","LAT")],data=active.wx,proj4string =CRS("+init=epsg:4269")) 
active.wx=spTransform(active.wx,SFWMD.proj)

rf.sites=active.wx[n.ok.dis,]
rf.sites.over=point.in.poly(rf.sites,n.ok.dis)

tmap_mode("view")
tm_shape(basins.all)+tm_polygons()+  
  tm_shape(n.ok.dis)+tm_polygons(alpha=0.5,col="grey",border.col = "red")+
    tm_shape(active.wx)+tm_dots()

tmap_mode("plot")
wx_maps=tm_shape(n.ok.dis)+tm_polygons(alpha=0.5)+
  tm_shape(active.wx)+tm_symbols(shape=21,size=0.1)+
  tm_shape(rf.sites)+tm_symbols(shape=21,size=0.1,col="red")
wx_maps

#tmap_save(wx_maps,paste0(plot.path,"WX_basin_maps.png"),height=4,units="in")

### Rainfall Data
sdate=as.Date("1997-05-01")
edate=as.Date("2018-04-30")

rf.sites
#paste(rf.sites@data$STATION,collapse="/")
dbkey.ts.list=read.csv(paste0(data.path,"timeseries_list.csv"));#from DBHYDRO
dbkey.ts.list$Start.Date2=date.fun(as.POSIXct(as.character(dbkey.ts.list$Start.Date),format="%m/%d/%Y"))
dbkey.ts.list$End.Date2=date.fun(as.POSIXct(as.character(dbkey.ts.list$End.Date),format="%m/%d/%Y"))
dbkey.ts.list=subset(dbkey.ts.list,Data.Type=="RAIN"&Frequency=="DA")
dbkey.ts.list=subset(dbkey.ts.list,Start.Date2>=date.fun(sdate))

dbkey.ts.list=merge(dbkey.ts.list,rf.sites.over@data[,c("STATION","Basin_gen")],by.x="Station",by.y="STATION")
head(dbkey.ts.list)

subset(dbkey.ts.list,Basin_gen=="Fisheating Creek")

rf.dat.all=data.frame()
for(i in 1:nrow(dbkey.ts.list)){
  tmp=SFWMD.DBHYDRO.Data.daily(sdate,edate,dbkey.ts.list$Dbkey[i])
  tmp$DBKEY=as.character(dbkey.ts.list$Dbkey[i])
  rf.dat.all=rbind(rf.dat.all,tmp)
  print(paste0(i,": ",dbkey.ts.list$Dbkey[i]))
}
head(rf.dat.all)
rf.dat.all=merge(rf.dat.all,rf.sites.over@data[,c("STATION","Basin_gen")],by.x="Station",by.y="STATION")     
head(rf.dat.all)
rf.dat.all$DATE.EST=date.fun(rf.dat.all$Date)
rf.dat.basin=ddply(rf.dat.all,c("Basin_gen","DATE.EST"),summarise,RF.in=mean(Data.Value,na.rm=T))
rf.dat.basin$month=as.numeric(format(rf.dat.basin$DATE.EST,"%m"))
rf.dat.basin$CY=as.numeric(format(rf.dat.basin$DATE.EST,"%Y"))
rf.dat.basin$WY=WY(rf.dat.basin$DATE.EST)
rf.dat.basin$Hydro=FL.Hydroseason(rf.dat.basin$DATE.EST)
#write.csv(rf.dat.basin,paste0(export.path,"/basin_daily_RF.csv"),row.names=F)
#attributes(rf.dat.basin$DATE.EST)

ylim.val=c(0,8);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(RF.in*-1~DATE.EST,subset(rf.dat.basin,Basin_gen=="Lower Kissimmee"),type="h",yaxt="n")
with(subset(rf.dat.basin,Basin_gen=="Lower Kissimmee"),points(DATE.EST,RF.in*-1,pch=19,cex=0.8))
axis_fun(2,ymaj*-1,ymin*-1,ymaj,1)

dates.val=seq(date.fun(sdate),date.fun(edate),"1 months")
basin.val=sort(rep(as.character(unique(rf.dat.basin$Basin_gen)),length(dates.val)))
fill.dat=data.frame(date.moncy=rep(dates.val,length(unique(rf.dat.basin$Basin_gen))),Basin_gen=basin.val,fill=1)

rf.dat.mon=ddply(rf.dat.basin,c("CY","month","WY","Basin_gen"),summarise,TRF.in=sum(RF.in,na.rm=T))
rf.dat.mon$date.moncy=with(rf.dat.mon,date.fun(paste(CY,month,"1",sep="-")))
rf.dat.mon=merge(rf.dat.mon,fill.dat,c("date.moncy","Basin_gen"),all.y=T)
rf.dat.mon=rf.dat.mon[order(rf.dat.mon$Basin_gen,rf.dat.mon$date.moncy),]

plot(RF.cm~date.moncy,subset(rf.dat.mon,Basin_gen=="NubbinSlough_TaylorCreek"),type="l")

rf.dat.mon$WWI=with(rf.dat.mon,ave(RF.cm,Basin_gen,FUN=function(x)WWI.index.fun(x)))
rf.dat.mon$WWI_6mon=with(rf.dat.mon,ave(WWI,Basin_gen,FUN=function(x)c(rep(NA,5),roll_mean(x,6))))

plot(WWI~date.moncy,subset(rf.dat.mon,Basin_gen=="Lower Kissimmee"&WY==2015),type="l")
hist(subset(rf.dat.mon,Basin_gen=="Lower Kissimmee"&WY==2015)$WWI)
hist(subset(rf.dat.mon,Basin_gen=="Lower Kissimmee"&WY==2018)$WWI)

###
##
#dbkeys=data.frame(SITE=c("S65E","S65","S191"),DBKEY=c("15631","H0289","15639"),Priority=c("P1","P1","P1"),WQ.Site=c("S65E","S65","S191"))
#dbkeys=subset(dbkeys,SITE=="S65E")
flow.dbkeys=data.frame(SITE=c(rep("S65",2),rep("S65A",2),rep("S65B",2),rep("S65C",2),"S65CX",rep("S65D",2),"S65DX1",rep("S65DX2",2),rep("S65E",2),"S65EX1"),
                       DBKEY=c("H0289","WN230","J9202","91646","HG238","91649","04458","91651","15332","WN344","91655","88300","AI354","91653","P1020","KO585","AL760"),
                       Priority=c(rep(c("P1","P2"),4),"P1",c("P1","P2"),"P1",rep(c("P1","P2"),2),"P1"),
                       WQ=c(rep("S65",2),rep("S65A",2),rep("S65B",2),rep("S65C",2),"KEA101",rep("S65D",2),"S65DX1",rep("S65E",2),rep("S65E",3)),
                       Struct=c(rep("S65",2),rep("S65A",2),rep("S65B",2),rep("S65C",3),rep("S65D",5),rep("S65E",3)))
flow.dbkeys=subset(flow.dbkeys,SITE!="S65B")
flow.dbkeys=subset(flow.dbkeys,DBKEY!="04458")
flow.dbkeys$Priority=with(flow.dbkeys,ifelse(DBKEY=="15338","P1",as.character(Priority)))

flow.dbkeys=subset(flow.dbkeys,Struct=="S65E")

#Flow Data
flow=data.frame()
for(i in 1:nrow(flow.dbkeys)){
  tmp=SFWMD.DBHYDRO.Data.daily(sdate,edate,flow.dbkeys$DBKEY[i])
  tmp$DBKEY=as.character(flow.dbkeys$DBKEY[i])
  flow=rbind(tmp,flow)
  print(paste(i,":",flow.dbkeys$SITE[i],"(",flow.dbkeys$DBKEY[i],")"))
}
flow.data=merge(flow,flow.dbkeys,"DBKEY")
flow.data$WY=WY(flow.data$Date)

flow.xtab=data.frame(cast(flow.data,Date+WY+SITE+Struct+WQ~Priority,value="Data.Value",fun.aggregate=function(x) ifelse(sum(x,na.rm=T)==0,NA,sum(x,na.rm=T))))
flow.xtab$P3=as.numeric(NA)
flow.xtab$Date.EST=date.fun(flow.xtab$Date)
flow.xtab$fflow.cfs=with(flow.xtab,ifelse(is.na(P1)==T&is.na(P2)==T,P3,ifelse(is.na(P2)==T&is.na(P3)==T,P1,ifelse(is.na(P1)==T&is.na(P3)==T,P2,P1))));#final flow value for analysis
range(flow.xtab$fflow.cfs,na.rm=T)
flow.xtab$flow.kacft=cfs.to.acftd(flow.xtab$fflow.cfs)/1000

##
## WQ Data Download

wq.sites=ddply(flow.dbkeys,c("WQ","Struct"),summarise,N.val=N(SITE))[,c("WQ","Struct")]
wq.param=data.frame(Test.Number=c(16,18,20,21,23,25,80,89,100),
                    Param=c("TSS","NOx","NH4","TKN","OP","TP","TN","DOC","TOC"))

wq.dat=SFWMD.DBHYDRO.Data.WQ(sdate,edate,wq.sites$WQ,subset(wq.param,Param=="TP")$Test.Number)
#ddply(wq.dat,c("Test.Number","Test.Name"),summarise,N.val=N(HalfMDL))
unique(wq.dat$Collection.Method)
wq.dat$sample.date=with(wq.dat,ifelse(Collection.Method%in%c("ACT","ACF")&is.na(First.Trigger.Date)==F,as.character(First.Trigger.Date),as.character(Collection_Date)))
wq.dat$sample.date=date.fun(wq.dat$sample.date,tz="America/New_York")
wq.dat=merge(wq.dat,wq.sites,by.x="Station.ID",by.y="WQ")
wq.dat=merge(wq.dat,wq.param,"Test.Number")
wq.dat$HalfMDL.ugL=wq.dat$HalfMDL*1000

wq.dat.xtab=data.frame(cast(wq.dat,Struct+Station.ID+Param+sample.date~Collection.Method,value="HalfMDL",mean))
wq.dat.xtab$ACF.flag=ifelse(is.nan(wq.dat.xtab$ACF),0,1)
wq.dat.xtab$ACT.flag=ifelse(is.nan(wq.dat.xtab$ACT),0,1)
wq.dat.xtab$ADT.flag=ifelse(is.nan(wq.dat.xtab$ADT),0,1)
wq.dat.xtab$auto.flag=rowSums(wq.dat.xtab[,c("ACF.flag","ACT.flag","ADT.flag")])
wq.dat.xtab$G.flag=ifelse(is.nan(wq.dat.xtab$G),0,1)
subset(wq.dat.xtab,ADT.flag==1&ACF.flag==1&ACT.flag==1)
subset(wq.dat.xtab,auto.flag>1)
wq.dat.xtab$sumflag=rowSums(wq.dat.xtab[,c("auto.flag","G.flag")])
wq.dat.xtab$auto=with(wq.dat.xtab,ifelse(is.nan(ACF),ACT,ACF))
#wq.dat.xtab$auto=with(wq.dat.xtab,ifelse(is.nan(ACF),ADT,ACF))
wq.dat.xtab$finalwq=with(wq.dat.xtab,ifelse(is.nan(auto),G,auto))
wq.dat.xtab$Date.EST=date.fun(wq.dat.xtab$sample.date)
wq.dat.xtab$WY=WY(wq.dat.xtab$Date.EST)

flow.vars=c("Date.EST","WY","SITE","Struct","WQ","fflow.cfs")
wq.vars=c("Date.EST","Struct","Station.ID","Param","finalwq","G","ACF","ACT","ADT")

flow.wq.TP=merge(flow.xtab[,flow.vars],subset(wq.dat.xtab,Param=="TP")[,wq.vars],by.x=c("Date.EST","Struct","WQ"),by.y=c("Date.EST","Struct","Station.ID"),all.x=T)
flow.wq.TP$finalwq.int=with(flow.wq.TP,ave(finalwq,WQ,FUN=function(x)dat.interp(x)))
flow.wq.TP$Load.kg=with(flow.wq.TP,Load.Calc.kg(fflow.cfs,finalwq.int))
flow.wq.struct.sum.TP=ddply(flow.wq.TP,c("Struct","Date.EST","WY"),summarise,Q.acft=sum(cfs.to.acftd(fflow.cfs)),TLoad.kg=sum(Load.kg,na.rm=T),TLoad.mt=sum(kg.to.mt(Load.kg),na.rm=T))
flow.wq.struct.sum.TP$hydro.day=hydro.day(flow.wq.struct.sum.TP$Date.EST,"FL")
flow.wq.struct.sum.TP$DOY=as.numeric(format(flow.wq.struct.sum.TP$Date.EST,"%j"))
flow.wq.struct.sum.TP$Q.kacft.plot=with(flow.wq.struct.sum.TP,ifelse(is.na(Q.acft)==T,0,Q.acft/1000))

###
###
rainfall.lk=subset(rf.dat.basin,Basin_gen=="Lower Kissimmee")
WYs=seq(2003,2018,1)
season.start=date.fun(c(paste0(WYs,"-05-01"),paste0(WYs,"-11-01")))
#tiff(filename=paste0(plot.path,"/S65E_RF_TP_FLOW_ALL.tiff"),width=10,height=4.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",cex.axis=1,mar=c(1,2,0.75,2),oma=c(4,2,0.5,2))
layout(matrix(1:2,2,1,byrow=T),heights=c(0.3,0.7))

ylim.val.rf=c(0,9);by.y=2;ymaj.rf=seq(ylim.val.rf[1],ylim.val.rf[2],by.y);ymin.rf=seq(ylim.val.rf[1],ylim.val.rf[2],by.y/2)
ylim.val.q=c(0,45);by.y=10;ymaj.q=seq(ylim.val.q[1],ylim.val.q[2],by.y);ymin.q=seq(ylim.val.q[1],ylim.val.q[2],by.y/2)
ylim.val.wq=c(0,450);by.y=100;ymaj.wq=seq(ylim.val.wq[1],ylim.val.wq[2],by.y);ymin.wq=seq(ylim.val.wq[1],ylim.val.wq[2],by.y/2)
xlim.val=date.fun(c("2003-05-01","2018-05-01"));xmaj=seq(xlim.val[1],xlim.val[2],"12 months");xmin=seq(xlim.val[1],xlim.val[2],"6 months")

plot(RF.in*-1~DATE.EST,rainfall.lk,type="n",yaxt="n",lwd=1,ylim=rev(ylim.val.rf*-1),xlim=xlim.val,ylab=NA,xlab=NA,xaxt="n",xaxs="i",yaxs="i",col="dodgerblue")
abline(h=ymaj.rf*-1,v=xmaj,lty=3,col="grey")
for(i in 1:nrow(rainfall.lk)){with(rainfall.lk[i,],lines(rep(DATE.EST,2),c(RF.in*-1,0),col="dodgerblue"))}
abline(v=season.start,lty=2,col="red")
axis_fun(2,ymaj.rf*-1,ymin.rf*-1,ymaj.rf,0.8)
axis_fun(1,xmaj,xmin,NA,1);box(lwd=1)
mtext(side=2,line=2,"Rainfall (inches)")

plot(Q.acft~Date.EST,flow.wq.struct.sum.TP,type="n",xaxt="n",yaxt="n",ylim=ylim.val.q,xlim=xlim.val,ylab=NA,xlab=NA,yaxs="i",xaxs="i")
abline(h=ymaj.q,v=xmaj,lty=3,col="grey")
with(subset(flow.wq.struct.sum.TP,Struct=="S65E"),shaded.range(Date.EST,rep(0,length(Date.EST)),Q.kacft.plot,bg="grey",col="dodgerblue",lty=1))
axis_fun(4,ymaj.q,ymin.q,ymaj.q,0.9)
par(new=T)
plot(HalfMDL.ugL~Date.EST,wq.dat,type="n",xaxt="n",yaxt="n",ylim=ylim.val.wq,xlim=xlim.val,ylab=NA,xlab=NA,xaxs="i",yaxs="i")
with(subset(wq.dat,Collection.Method=="ACF"),points(Date.EST,HalfMDL.ugL,pch=21,bg="grey",cex=0.8))
with(subset(wq.dat,Collection.Method=="G"),points(Date.EST,HalfMDL.ugL,pch=21,bg="black",cex=0.8))
abline(v=season.start,lty=2,col="red")
axis_fun(2,ymaj.wq,ymin.wq,ymaj.wq,0.9)
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m/%Y"),0.9);box(lwd=1)
mtext(side=2,line=2.25,"TP Concentration (\u03BCg L\u207B\u00B9)")
mtext(side=4,line=2,"Daily Flow (kac-ft day\u207B\u00B9)")
mtext(side=1,line=2,"Date (Month-Year)")
leg.x=xlim.val[1]+diff(xlim.val)/2
leg.y=ylim.val.wq[1]-125
legend.text=c("TP Grab","TP Auto (ACF)  ","Flow")
pt.cols=c("black","white")
legend(leg.x,leg.y,legend=legend.text,pch=c(21,21,NA),pt.bg=pt.cols,col=c("black","black","dodgerblue1"),lty=c(0,0,1),lwd=1,pt.cex=1.5,ncol=3,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend(leg.x,leg.y,legend=legend.text,pch=c(21,21,NA),pt.bg=pt.cols,col=c("black","black","dodgerblue1"),lty=c(0,0,0),lwd=1,pt.cex=1.5,ncol=3,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()

##
##
##
##
##
##
kr.basin.rf=subset(rf.dat.basin,Basin_gen=="Lower Kissimmee")
kr.basin.rf$RF.cm=in.to.cm(kr.basin.rf$RF.in)

wq.dat.da=ddply(subset(wq.dat,Collection.Method=="G"),c("Date.EST","Station.ID","Test.Name"),summarise,mean.TP=mean(HalfMDL*1000,na.rm=T))
wq.dat.da$date.diff=c(NA,diff(wq.dat.da$Date.EST))

rftot.wq.samp=data.frame()
for(i in 2:nrow(wq.dat)){
  date.seq=seq(wq.dat.da$Date.EST[i-1],wq.dat.da$Date.EST[i],"1 days")
  TRF=sum(subset(kr.basin.rf,DATE.EST%in%date.seq)$RF.cm,na.rm=T)
  Tflow=sum(subset(flow.wq.struct.sum.TP,Date.EST%in%date.seq)$Q.acft,na.rm=T)
  tmp=data.frame(RF.cm=TRF,Q.acft=Tflow,Date.EST=wq.dat.da$Date.EST[i])
  rftot.wq.samp=rbind(rftot.wq.samp,tmp)
}

wq.da.rf=merge(wq.dat.da,rftot.wq.samp,"Date.EST")
wq.da.rf$Hydro=FL.Hydroseason(wq.da.rf$Date.EST)
wq.da.rf$WWI=WWI.index.fun(wq.da.rf$RF.cm)
plot(mean.TP~WWI,wq.da.rf)
with(wq.da.rf,cor.test(mean.TP,WWI,method="spearman"))


#Bin the data
#test$TP.bin=findInterval(test$mean.TP,seq(0,600,50))
#test$RF.bin=findInterval(test$RF.cm,seq(0,30,5))

plot(mean.TP~Q.acft,wq.da.rf)
plot(mean.TP~Q.acft,subset(wq.da.rf,Hydro=="A_Wet"))
plot(mean.TP~Q.acft,subset(wq.da.rf,Hydro=="B_Dry"))
with(subset(wq.da.rf,Hydro=="B_Dry"),cor.test(mean.TP,Q.acft,method="spearman"))
with(subset(wq.da.rf,Hydro=="A_Wet"),cor.test(mean.TP,Q.acft,method="spearman"))
with(wq.da.rf,cor.test(mean.TP,Q.acft,method="spearman"))
mod.q=mblm(mean.TP~Q.acft,wq.da.rf)
fit.mod.q=data.frame(x.val=seq(0,4.25e5,100),predict(mod.q,data.frame(Q.acft=seq(0,4.25e5,100)),interval="confidence"))


plot(mean.TP~RF.cm,wq.da.rf)
plot(mean.TP~RF.cm,subset(wq.da.rf,Hydro=="A_Wet"))
plot(mean.TP~RF.cm,subset(wq.da.rf,Hydro=="B_Dry"))
with(subset(wq.da.rf,Hydro=="B_Dry"),cor.test(mean.TP,RF.cm,method="spearman"))
with(subset(wq.da.rf,Hydro=="A_Wet"),cor.test(mean.TP,RF.cm,method="spearman"))
with(wq.da.rf,cor.test(mean.TP,RF.cm,method="spearman"))
mod=mblm(mean.TP~RF.cm,subset(wq.da.rf,is.na(RF.cm)==F))
fit.mod=data.frame(x.val=seq(0,30,0.5),predict(mod,data.frame(RF.cm=seq(0,30,0.5)),interval="confidence"))



xlim.val=c(0,4.25e5);by.x=5e4;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,600);by.y=100;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2);
#tiff(filename=paste0(plot.path,"/antecedent_Q_wq.tiff"),width=5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.5,0.75,0.5),oma=c(4,2,0.75,0.25),mgp=c(3,1,0));

plot(mean.TP~Q.acft,wq.da.rf,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(subset(wq.da.rf,Hydro=="B_Dry"),points(Q.acft,mean.TP,pch=22,bg=adjustcolor("khaki",0.75),lwd=0.1,cex=0.9))
with(subset(wq.da.rf,Hydro=="A_Wet"),points(Q.acft,mean.TP,pch=21,bg=adjustcolor("dodgerblue",0.5),lwd=0.1,cex=0.9))
#with(fit.mod.q,shaded.range(x.val,lwr,upr,"grey30"))
#with(fit.mod.q,lines(x.val,fit,lwd=2.5,lty=2))
axis_fun(1,xmaj,xmin,xmaj/1000,1,line=-0.5);axis_fun(2,ymaj,ymin,ymaj,1);box(lwd=1)
mtext(side=1,line=1.5,"Antecedent Total Flow Volume (k-AcFt)")
mtext(side=2,line=2.25,"Total Phosphorus Concentration (\u03BCg L\u207B\u00B9)")
leg.x=xlim.val[1]+diff(xlim.val)/2
leg.y=ylim.val[1]-125
pt.cols=c(adjustcolor(c("khaki","dodgerblue"),0.5),adjustcolor("grey30",0.25))
legend.text=c("Dry Season","Wet Season")#, "Thiel-Sen 95% CI", "Thiel-Sen")
#legend(leg.x,leg.y,legend=legend.text,pch=c(22,21,22,NA),pt.bg=pt.cols,col="black",lty=c(0,0,0,2),lwd=1,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
#legend(leg.x,leg.y,legend=legend.text,pch=c(22,21,22,NA),pt.bg=pt.cols,col="black",lty=c(0,0,0,0),lwd=1,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
legend(leg.x,leg.y,legend=legend.text,pch=c(22,21),pt.bg=pt.cols,col="black",lty=c(0,0),lwd=1,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend(leg.x,leg.y,legend=legend.text,pch=c(22,21),pt.bg=pt.cols,col="black",lty=c(0,0),lwd=1,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)

dev.off()


xlim.val=c(0,30);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,600);by.y=100;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2);
#tiff(filename=paste0(plot.path,"/antecedent_RF_wq.tiff"),width=5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.5,0.75,0.5),oma=c(4,2,0.75,0.25),mgp=c(3,1,0));

plot(mean.TP~RF.cm,wq.da.rf,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(subset(wq.da.rf,Hydro=="B_Dry"),points(RF.cm,mean.TP,pch=22,bg=adjustcolor("khaki",0.75),lwd=0.1,cex=0.9))
with(subset(wq.da.rf,Hydro=="A_Wet"),points(RF.cm,mean.TP,pch=21,bg=adjustcolor("dodgerblue",0.5),lwd=0.1,cex=0.9))
with(fit.mod,shaded.range(x.val,lwr,upr,"grey30"))
with(fit.mod,lines(x.val,fit,lwd=2.5,lty=2))
axis_fun(1,xmaj,xmin,xmaj,1,line=-0.5);axis_fun(2,ymaj,ymin,ymaj,1);box(lwd=1)
mtext(side=1,line=1.5,"Antecedent Total Rain (cm)")
mtext(side=2,line=2.25,"Total Phosphorus Concentration (\u03BCg L\u207B\u00B9)")
leg.x=xlim.val[1]+diff(xlim.val)/2
leg.y=ylim.val[1]-125
pt.cols=c(adjustcolor(c("khaki","dodgerblue"),0.5),adjustcolor("grey30",0.25))
legend.text=c("Dry Season","Wet Season", "Thiel-Sen 95% CI", "Thiel-Sen")
legend(leg.x,leg.y,legend=legend.text,pch=c(22,21,22,NA),pt.bg=pt.cols,col="black",lty=c(0,0,0,2),lwd=1,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend(leg.x,leg.y,legend=legend.text,pch=c(22,21,22,NA),pt.bg=pt.cols,col="black",lty=c(0,0,0,0),lwd=1,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()

####
####
####
####
## Seasonal Rainfall - FWM Comparison
flow.wq.TP$Hydro=FL.Hydroseason(flow.wq.TP$Date.EST)

rf.vars=c("DATE.EST","WY","Hydro","RF.cm")
flow.wq.rf=merge(flow.wq.TP,kr.basin.rf[,rf.vars],by.x=c("Date.EST","WY","Hydro"),by.y=c("DATE.EST","WY","Hydro"))
flow.wq.rf.WY=ddply(flow.wq.rf,c("WY","Hydro"),summarise,TLoad.kg=sum(Load.kg,na.rm=T),TFlow.cmd=sum(cfs.to.m3d(fflow.cfs),na.rm=T),TRF.cm=sum(RF.cm,na.rm=T))
flow.wq.rf.WY$FWM.ugL=with(flow.wq.rf.WY,(TLoad.kg*1e9)/(TFlow.cmd*1000))
flow.wq.rf.WY$WWI=WWI.index.fun(flow.wq.rf.WY$TRF.cm)

plot(FWM.ugL~WWI,flow.wq.rf.WY)
with(flow.wq.rf.WY,cor.test(FWM.ugL,WWI,method="spearman"))
with(subset(flow.wq.rf.WY,Hydro=="B_Dry"),cor.test(FWM.ugL,WWI,method="spearman"))
with(subset(flow.wq.rf.WY,Hydro=="A_Wet"),cor.test(FWM.ugL,WWI,method="spearman"))
mod.wwi=mblm(FWM.ugL~WWI,flow.wq.rf.WY)
fit.mod.wwi=data.frame(x.val=seq(-0.75,0.75,0.1),predict(mod.wwi,data.frame(WWI=seq(-0.75,0.75,0.1)),interval="confidence"))



plot(FWM.ugL~TRF.cm,flow.wq.rf.WY)
with(flow.wq.rf.WY,cor.test(FWM.ugL,TRF.cm,method="spearman"))
with(subset(flow.wq.rf.WY,Hydro=="B_Dry"),cor.test(FWM.ugL,TRF.cm,method="spearman"))
with(subset(flow.wq.rf.WY,Hydro=="A_Wet"),cor.test(FWM.ugL,TRF.cm,method="spearman"))
mod.fwm=mblm(FWM.ugL~TRF.cm,flow.wq.rf.WY)
fit.mod.fwm=data.frame(x.val=seq(10,140,10),predict(mod.fwm,data.frame(TRF.cm=seq(10,140,10)),interval="confidence"))


flow.wq.rf.WY$WY.plot=with(flow.wq.rf.WY,ifelse(Hydro=="A_Wet",WY-0.15,WY+0.15))
WYs=seq(1998,2018,1)
xlim.val=c(1998,2018);by.x=3;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(50,250);by.y=50;ymaj=seq(min(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(min(c(0,ylim.val[1])),ylim.val[2],by.y/2);
#tiff(filename=paste0(plot.path,"/FWM_hydro.tiff"),width=5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.5,0.5,0.5),oma=c(3.5,2,0.5,0.25),mgp=c(3,1,0));

plot(FWM.ugL~WY,flow.wq.rf.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
#with(flow.wq.rf.WY,lines(WY.plot,FWM.ugL,lty=1,lwd=1.25))
#for(i in 1:length(WYs)){with(subset(flow.wq.rf.WY,WY==WYs[i]),lines(WY.plot,FWM.ugL,lty=2,lwd=1))}
for(i in 1:nrow(flow.wq.rf.WY)){with(flow.wq.rf.WY[i,],lines(rep(WY.plot,2),c(0,FWM.ugL),lty=1,lwd=2.5,col=adjustcolor(ifelse(Hydro=="A_Wet","dodgerblue1","khaki"),0.75)))}
with(subset(flow.wq.rf.WY,Hydro=="A_Wet"),points(WY.plot,FWM.ugL,pch=21,bg="dodgerblue",lwd=0.25,cex=1))
with(subset(flow.wq.rf.WY,Hydro=="B_Dry"),points(WY.plot,FWM.ugL,pch=21,bg="khaki",lwd=0.25,cex=1))
axis_fun(1,xmaj,xmin,xmaj,1,line=-0.5);axis_fun(2,ymaj,ymin,ymaj,1);box(lwd=1)
mtext(side=1,line=1.5,"Water Year")
mtext(side=2,line=2.3,"TP FWM Concentration (\u03BCg L\u207B\u00B9)")
leg.x=xlim.val[1]+diff(xlim.val)/2
leg.y=ylim.val[1]-40
pt.cols=c(adjustcolor(c("dodgerblue","khaki"),0.5))
legend.text=c("Wet Season","Dry Season")
legend(leg.x,leg.y,legend=legend.text,pch=c(21,21),pt.bg=pt.cols,col="black",lty=c(0,0),lwd=1,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend(leg.x,leg.y,legend=legend.text,pch=c(21,21),pt.bg=pt.cols,col="black",lty=c(0,0),lwd=1,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()


xlim.val=c(10,135);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(40,250);by.y=50;ymaj=seq(min(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(min(c(0,ylim.val[1])),ylim.val[2],by.y/2);
#tiff(filename=paste0(plot.path,"/RF_FWM_hydro.tiff"),width=5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.5,0.5,0.5),oma=c(3.5,2,0.5,0.25),mgp=c(3,1,0));

plot(FWM.ugL~TRF.cm,flow.wq.rf.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(fit.mod.fwm,shaded.range(x.val,lwr,upr,"grey30"))
with(fit.mod.fwm,lines(x.val,fit,lwd=2.5,lty=2))
with(subset(flow.wq.rf.WY,Hydro=="B_Dry"),points(TRF.cm,FWM.ugL,pch=22,bg=adjustcolor("khaki",0.75),lwd=0.1,cex=1))
with(subset(flow.wq.rf.WY,Hydro=="A_Wet"),points(TRF.cm,FWM.ugL,pch=21,bg=adjustcolor("dodgerblue",0.5),lwd=0.1,cex=1))
axis_fun(1,xmaj,xmin,xmaj,1,line=-0.5);axis_fun(2,ymaj,ymin,ymaj,1);box(lwd=1)
mtext(side=1,line=1.5,"Total Rain (cm)")
mtext(side=2,line=2.25,"TP FWM Concentration (\u03BCg L\u207B\u00B9)")
leg.x=xlim.val[1]+diff(xlim.val)/2
leg.y=ylim.val[1]-40
pt.cols=c(adjustcolor(c("khaki","dodgerblue"),0.5),adjustcolor("grey30",0.25))
legend.text=c("Dry Season","Wet Season", "Thiel-Sen 95% CI", "Thiel-Sen")
legend(leg.x,leg.y,legend=legend.text,pch=c(22,21,22,NA),pt.bg=pt.cols,col="black",lty=c(0,0,0,2),lwd=1,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend(leg.x,leg.y,legend=legend.text,pch=c(22,21,22,NA),pt.bg=pt.cols,col="black",lty=c(0,0,0,0),lwd=1,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()

xlim.val=c(-0.75,0.75);by.x=0.25;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(40,250);by.y=50;ymaj=seq(min(c(0,ylim.val[1])),ylim.val[2],by.y);ymin=seq(min(c(0,ylim.val[1])),ylim.val[2],by.y/2);
#tiff(filename=paste0(plot.path,"/WWI_FWM_hydro.tiff"),width=5,height=4,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.5,0.5,0.5),oma=c(3.5,2,0.5,0.25),mgp=c(3,1,0));

plot(FWM.ugL~TRF.cm,flow.wq.rf.WY,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(fit.mod.wwi,shaded.range(x.val,lwr,upr,"grey30"))
with(fit.mod.wwi,lines(x.val,fit,lwd=2.5,lty=2))
with(subset(flow.wq.rf.WY,Hydro=="B_Dry"),points(WWI,FWM.ugL,pch=22,bg=adjustcolor("khaki",0.75),lwd=0.1,cex=1))
with(subset(flow.wq.rf.WY,Hydro=="A_Wet"),points(WWI,FWM.ugL,pch=21,bg=adjustcolor("dodgerblue",0.5),lwd=0.1,cex=1))
axis_fun(1,xmaj,xmin,xmaj,1,line=-0.5);axis_fun(2,ymaj,ymin,ymaj,1);box(lwd=1)
mtext(side=1,line=1.5,"Weather Whiplash Index")
mtext(side=2,line=2.25,"TP FWM Concentration (\u03BCg L\u207B\u00B9)")
leg.x=xlim.val[1]+diff(xlim.val)/2
leg.y=ylim.val[1]-40
pt.cols=c(adjustcolor(c("khaki","dodgerblue"),0.5))
legend.text=c("Dry Season","Wet Season", "Thiel-Sen")
legend(leg.x,leg.y,legend=legend.text,pch=c(22,21,NA),pt.bg=pt.cols,col="black",lty=c(0,0,2),lwd=1,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend(leg.x,leg.y,legend=legend.text,pch=c(22,21,NA),pt.bg=pt.cols,col="black",lty=c(0,0,0),lwd=1,pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()


####
####
####
####
####
####
####
# OLD CODE

wq.mon=ddply(wq.dat,c("CY","month","WY","Station.ID"),summarise,mean.TP.ugL=mean(HalfMDL*1000,na.rm=T))
wq.mon$date.moncy=with(wq.mon,date.fun(paste(CY,month,"1",sep="-")))
plot(mean.TP.ugL~date.moncy,wq.mon)

vars=c("date.moncy","WY","TRF.in","WWI")
vars2=c("date.moncy","mean.TP.ugL")
WWI.TP=merge(subset(rf.dat.mon,Basin_gen=="Lower Kissimmee")[,vars],wq.mon[,vars2],"date.moncy")
WWI.TP$trans.cat=with(WWI.TP,ifelse(WWI<0,"dry-wet","wet-dry"))


plot(mean.TP.ugL~TRF.in,WWI.TP)
with(WWI.TP,cor.test(mean.TP.ugL,RF.cm,method="spearman"))
mod=mblm(mean.TP.ugL~TRF.in,WWI.TP)
mod.pred=data.frame(x.val=seq(0,20,0.5),predict(mod,data.frame(TRF.in=seq(0,20,0.5)),interval="prediction"))
with(mod.pred,lines(x.val,fit,lty=2,lwd=2))
library(psych)

with(WWI.TP,ellipses(RF.cm,mean.TP.ugL,smooth=F,lm=F,n=1,add=T))

plot(mean.TP.ugL~WWI,WWI.TP,log="y")
with(WWI.TP,cor.test(mean.TP.ugL,WWI,method="spearman"))
mod=mblm(mean.TP.ugL~WWI,WWI.TP)
mod.pred=data.frame(x.val=seq(-1,1,0.1),predict(mod,data.frame(WWI=seq(-1,1,0.1)),interval="prediction"))
with(mod.pred,lines(x.val,fit,lty=2,lwd=2))


boxplot(mean.TP.ugL~trans.cat,WWI.TP,outline=F)
plot(RF.cm~date.moncy,subset(rf.dat.mon,Basin_gen=="Lower Kissimmee"),type="l")
par(new=T);
plot(HalfMDL*1000~Date.EST,wq.dat,pch=19)

plot(WWI~date.moncy,subset(rf.dat.mon,Basin_gen=="Lower Kissimmee"),type="l")
with(subset(rf.dat.mon,Basin_gen=="Lower Kissimmee"),points(date.moncy,WWI,pch=21,bg=ifelse(WWI<0,"dodgerblue1","khaki")))
par(new=T);
plot(HalfMDL*1000~Date.EST,wq.dat,pch=19)

low.kr.rf=subset(rf.dat.mon,Basin_gen=="Lower Kissimmee")

xlim.val=date.fun(c(sdate,edate));xmaj=seq(xlim.val[1],xlim.val[2],"3 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
#tiff(filename=paste0(plot.path,"/lower_KR_WWI.tiff"),width=6,height=5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",mar=c(1,1.5,0.75,0.5),oma=c(3,2,0.75,0.25),mgp=c(3,1,0));
layout(matrix(1:2,2,1,byrow=F))

ylim.val=c(0,20);by.y=5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(RF.cm~date.moncy,low.kr.rf,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA,yaxs="i")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
for(i in 1:nrow(low.kr.rf)){
  with(subset(low.kr.rf,date.moncy==low.kr.rf$date.moncy[i]),lines(rep(date.moncy,2),c(0,RF.cm),lty=1,lwd=1,col="grey50"))
}
with(low.kr.rf,points(date.moncy,RF.cm,pch=19,cex=0.5))
axis_fun(1,xmaj,xmin,NA,1)
axis_fun(2,ymaj,ymin,format(ymaj),1);box(lwd=1)
mtext(side=2,"Monthly Rainfall Total (inches)",line=2.25)
mtext(side=3, "Lower Kissimmee")

ylim.val=c(-1,1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(WWI~date.moncy,low.kr.rf,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",type="n",ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(low.kr.rf,lines(date.moncy,WWI,lwd=1,lty=1,col="black"))
with(low.kr.rf,points(date.moncy,WWI,pch=21,bg=ifelse(WWI<0,"dodgerblue1","khaki"),cex=0.8))
#with(low.kr.rf,lines(date.moncy,WWI_6mon,col="red",lty=1))
axis_fun(1,xmaj,xmin,format(xmaj,"%m-%Y"),line=-0.50,0.8)
axis_fun(2,ymaj,ymin,format(ymaj),1);box(lwd=1)
mtext(side=2,"Weather Wiplash Index",line=2.25)
mtext(side=1,"Date (Month-Year)",line=2)
dev.off()

