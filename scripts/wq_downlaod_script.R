## St Lucie Estuary WQ
## Code was compiled by Paul Julian
## contact infor: pjulian@ufl.edu

#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

#Libraries
library(plyr)
library(reshape)

##Custom Functions
setwd("D:/UF/EstuaryWQ") # Different on your PC

#Function to download WQ data from DBHYDRO
SFWMD.DBHYDRO.Data.WQ=function(SDATE,EDATE,WQ.Station,tests){
  SDATE=paste(format(SDATE,"%d"),toupper(format(SDATE,"%b")),format(SDATE,"%Y"),sep="-");#"%d-%b-%Y
  EDATE=paste(format(EDATE,"%d"),toupper(format(EDATE,"%b")),format(EDATE,"%Y"),sep="-");#"%d-%b-%Y
  WQ.Station=paste("'",WQ.Station,"'",collapse=",",sep="")
  tests=paste("(",paste(tests,collapse=",",sep=""),")",sep="")
  WQ.link=paste0("http://my.sfwmd.gov/dbhydroplsql/water_quality_data.report_full?v_where_clause=where+station_id+in+(",WQ.Station,")+and+test_number+in+",tests,"+and+date_collected+>=+'",SDATE,"'+and+date_collected+<+'",EDATE,"'+and+sample_type_new+=+'SAMP'&v_exc_qc=Y&v_exc_flagged=Y&v_target_code=file_csv")
  REPORT=read.csv(WQ.link);
  REPORT=subset(REPORT,is.na(Test.Number)==F)
  REPORT$Collection_Date=as.POSIXct(strptime(REPORT$Collection_Date,"%d-%b-%Y %R"),tz="America/New_York")
  REPORT$First.Trigger.Date=as.POSIXct(strptime(REPORT$First.Trigger.Date,"%d-%b-%Y %R"),tz="America/New_York")
  REPORT$Date=as.POSIXct(strptime(REPORT$Collection_Date,"%F"),tz="America/New_York")
  REPORT$DateTime.EST=REPORT$Collection_Date
  attr(REPORT$DateTime.EST,"tzone")<-"EST"
  REPORT$Date.EST=as.POSIXct(strptime(REPORT$DateTime.EST,"%F"),tz="EST")
  #REPORT$Date.EST=REPORT$Date
  #attributes(REPORT$Date.EST)$tzone="EST"
  REPORT$HalfMDL=with(REPORT,ifelse(Value<0,abs(Value)/2,Value))
  return(REPORT)
}

## Molar ratios 
N.mw=14.0067
P.mw=30.973762
C.mw=12.0107

uf.cols.muted=c(rgb(108/255,154/255,195/255,1,"gator.blue"),rgb(226/255,143/255,65/255,1,"gator.orange"))
uf.cols=c(rgb(0/255,33/255,165/255,1,"gator.blue"),rgb(250/255,70/255,22/255,1,"gator.orange"))
tz.val="EST"

#Paths
paths=paste(getwd(),c("/export/","/plots/","/data/"),sep="")
export.path=paths[1]
plot.path=paths[2]
data.path=paths[3]

dates=as.Date(c("1978-05-01","2016-04-30"))
##
params=data.frame(parameter=c("NH4","TOC","Chla","Chla_cor","ChlB","ChlC","DO","TKN","NOx","pH","Pheo","DP","OrthoP","TP","Sal","SiO4","SPC","Temp","TSS","Turb","VSS","DOC","TN"),
                  Test.Number=c(20,100,61,112,62,113,8,21,18,10,64,26,23,25,98,27,9,7,6,12,77,89,80))
wq.sites=data.frame(StationID=c("SE 13","SE 12","SE 06","HR1","SE 08B","SE 03","SE 02","SE 01","SE 11","IRL12","IRL12B","IRL18B","S308C","C44S80","GORDYRD","C24S49","C23S48"),
                    ALIAS=c("SE13","SE12","sE06","HR1","SE08B","SE03","SE02","SE01","SE11","IRL12","IRL12","ITL18B","S308C","S80","GORDYRD","S49","S48"),
                    Region=c(rep("N_Fork",4),"S_Fork",rep("Estuary",3),"Marine",rep("IRL",3),rep("FW",2),rep("FW",3)))

pb=txtProgressBar(min=1,max=nrow(wq.sites),style=3)
wq.dat=NA
for( i in 1:nrow(wq.sites)){
  tmp=SFWMD.DBHYDRO.Data.WQ(dates[1],dates[2],wq.sites$StationID[i],params$Test.Number)
  wq.dat=rbind(wq.dat,tmp)
  setTxtProgressBar(pb,i)
}
wq.dat=merge(wq.dat,wq.sites,by.x="Station.ID",by.y="StationID")
#write.csv(wq.dat,paste(export.path,"DBHYDRO_wq_all.csv",sep=""),row.names = F)
