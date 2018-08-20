# Estuary whiplash: The role of climate and water management on estuary water quality and conditions.

This repository will be used to organize data, ideas and related information to the project.

## General Introduction
Loceke et al (2017) comprehensively reviewed how climate change (CC) can potentially drive and significanly deteriorate water quality in the intensive agricultural region of the North American Midwest. It is expected that CC will lead to increased frequency and severity of extreme climatic conditions such as drought and floods. This potential for extremes in weather conditions can lead to large swings in rainfall and subsequent surface water discharges. A transition from drought to flood can lead to the mobilization of nutrients and material to downstream systems and depedning on upstream landuses can potentially causing significant impacts to water quality conditions.  This drought-to-flood transition can happen rather rapid given the expectation of climate intensification resulting in a "weather whiplash" effect (Loceke et al 2017). 

Loceke et al (2017) using actual precipiation data and downscaled climate projection data developed a Weather Whiplash Index (WWI) and describes it as: 
>> _The weather whiplash index was calculated as the total precipitation from January to June of each year (1951 - 2099) minus the total precipiation from July to December of the previous year (1950 - 2098), divided by the total precipitation over that entire period._

Weather Whiplash Index provides a relative value to the shift in precipitation regimes between two different periods thereby providing a quantitative metric to CC severity. A positive WWI indicates shifts from dry to wet conditions, while negative WWI values indicate shifts from wet to dry. The relative value provides an estimate of shift severity.

Not stated by Loceke et al (2017), it is assumed that the January to June time period is considered a relatively dry period and the July to December is considered a relatively wet period for the American Midwest. Further inspection of [regional precipitation data](https://w2.weather.gov/climate/xmacis.php?wfo=eax) indicates that long-term mean rainfall between these periods are relatively equal. Below is a summary of regional precipitation values expressed as inches of precipitation.

| Calendar Year | Jan to Jun | Jul to Dec | WWI   | 
|:---------------:|:------------:|:------------:|:-------:| 
| 2000          | 18.36      | 16.60      |       | 
| 2001          | 29.91      | 23.59      | 0.29  | 
| 2002          | 16.36      | 8.41       | -0.18 | 
| 2003          | 15.94      | 12.01      | 0.31  | 
| 2004          | 19.57      | 18.02      | 0.24  | 
| 2005          | 23.97      | 20.17      | 0.14  | 
| 2006          | 10.02      | 20.85      | -0.34 | 
| 2007          | 18.28      | 14.74      | -0.07 | 
| 2008          | 19.59      | 25.07      | 0.14  | 
| 2009          | 22.31      | 22.64      | -0.06 | 
| 2010          | 20.29      | 21.62      | -0.05 | 
| 2011          | 18.56      | 18.36      | -0.08 | 
| 2012          | 13.14      | 9.14       | -0.17 | 
| 2013          | 18.05      | 16.43      | 0.33  | 
| 2014          | 15.53      | 24.51      | -0.03 | 
| 2015          | 24.18      | 22.41      | -0.01 | 
| 2016          | 21.51      | 27.14      | -0.02 | 
| 2017          | 22.70      | 23.32      | -0.09 | 
|               |            |            |       | 
| Mean          | 19.35      | 19.17      | 0.02  | 
| SD            | 4.52       | 5.42       | 0.19  | 



Below is the mathematical equivalent to the text explaining WWI by Loceke et al (2017) where P is total precipitation during each respective period. 

>><a href="https://www.codecogs.com/eqnedit.php?latex=WWI&space;=&space;\frac{P_{Dry,i}-P_{Wet,i-1}}{P_{Dry,i}&plus;P_{Wet,i-1}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?WWI&space;=&space;\frac{P_{Dry,i}-P_{Wet,i-1}}{P_{Dry,i}&plus;P_{Wet,i-1}}" title="WWI = \frac{P_{Dry,i}-P_{Wet,i-1}}{P_{Dry,i}+P_{Wet,i-1}}" /></a>

This approach can be adapted to other regions at any time scale. 

## Analysis in R
This analysis was performed in R (Ver 3.4.1). To faciliate analysis of rainfall data with the Weather Whiplash Index  a custom R-function was developed to have reproduceable results. Below is the WWI funtion using a `data.frame` with monthly total rainfall. 

```
WWI.index.fun=function(RF){
  require(zoo)
  dif.val=c(NA,diff(RF))
  tot.val=c(NA,rollsum(RF,2))
  WWI=dif.val/tot.val
  return(WWI)
}
```

## References
  + Loecke TD, Burgin AJ, Riveros-Iregui DA, et al (2017) Weather whiplash in agricultural regions drives deterioration of water quality. Biogeochemistry 133:7–15. doi: 10.1007/s10533-017-0315-z
