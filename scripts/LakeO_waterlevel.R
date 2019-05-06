## 
## Weather/management whiplash
## 
##
## Code was compiled by Paul Julian
## contact info: paul.julian@floridadep.gov

#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

#Libraries
library(AnalystHelper);#devtools::install_github("SwampThingPaul/AnalystHelper")
library(plyr)
library(reshape)
library(zoo);

#Paths
setwd("D:/UF/WeatherWhip")

paths=paste(getwd(),c("/Data/","/export/","/plots/","/GIS/"),sep="")
#Folder.Maker(paths);#One and done. Creates folders in working directory.
data.path=paths[1]
export.path=paths[2]
plot.path=paths[3]
gis.path=paths[4]

shell.exec(plot.path)

# Lake Okeechobee Water Level ---------------------------------------------

da.dbks=data.frame(SITE=c("L001","L005","L006","LZ40","S133TW","S352HW","S4TW"),DBKEY=c("16022","12509","12519","16265","15826","FF579","15732"),type="Daily")

dates=date.fun(c("1997-05-01","2017-05-01"))

WL.dat=data.frame()
for(i in 1:nrow(da.dbks)){
tmp=DBHYDRO_daily(dates[1],dates[2],da.dbks$DBKEY[i])
WL.dat=rbind(WL.dat,tmp)
print(i)
}
attributes(WL.dat$Date)
WL.dat$Date.EST=date.fun(WL.dat$Date)
WL.dat$WY=WY(WL.dat$Date.EST)
WL.dat=merge(WL.dat,da.dbks,"DBKEY",all.x=T)
#sum(is.na(WL.dat$SITE));#look for orphans
subset(WL.dat,Data.Value==0)
WL.dat$Data.Value[WL.dat$Data.Value==0]=NA

LakeO.xtab=cast(WL.dat,Date.EST+WY~SITE,value="Data.Value",mean)
LakeO.xtab$Mean=rowMeans(LakeO.xtab[,c("L001","L005","L006","LZ40","S133TW","S352HW","S4TW")],na.rm=T)
LakeO.xtab$Mean.m=ft.to.m(LakeO.xtab$Mean)
range(LakeO.xtab$Mean.m)

ylim.val=c(2.5,6.0);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
ylim.val2=c(8,19.5);by.y2=2;ymaj2=seq(ylim.val2[1],ylim.val2[2],by.y2);ymin2=seq(ylim.val2[1],ylim.val2[2],by.y2/2)
xlim.val=date.fun(dates);xmaj=seq(xlim.val[1],xlim.val[2],"4 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
#tiff(filename=paste0(plot.path,"LakeO_WaterLevel.tiff"),width=6,height=4.5,units="in",res=200,type="windows",compression=c("lzw"),bg="white")
par(family="serif",oma=c(2.5,1.75,1,0.25),mar=c(1.5,2,0.5,3))
plot(Mean.m~Date.EST,LakeO.xtab,ylim=ylim.val,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n",xaxs="i")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(LakeO.xtab,shaded.range(Date.EST,rep(0,N(Mean.m)),Mean.m,"dodgerblue1",lty=1))
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"))
axis_fun(2,ymaj,ymin,ymaj)
par(new=T);plot(Mean.m~Date.EST,LakeO.xtab,ylim=ylim.val2,xlim=xlim.val,yaxt="n",xaxt="n",ylab=NA,xlab=NA,type="n",xaxs="i")
axis_fun(4,ymaj2,ymin2,ymaj2)
mtext(side=2,line=2.25,"Water Level (m, NGVD29)")
mtext(side=4,line=1.75,"Water Level (ft, NGVD29)")
mtext(side=1,line=2,"Date (Month-Year)")
dev.off()