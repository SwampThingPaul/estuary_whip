## Project specific tools
# Weather Whiplash Index from
# Loecke TD, et al (2017)
## Weather whiplash in agricultural regions drives deterioration of water quality.
## Biogeochemistry 133:7–15.

WWI.index.fun=function(RF){
  require(zoo)
  dif.val=c(NA,diff(RF))
  tot.val=c(NA,rollapply(RF,width=2,FUN=sum,na.rm=T))
  WWI=dif.val/tot.val
  return(WWI)
}
