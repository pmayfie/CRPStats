library(dplyr)
library(sf)
library(terra)
library(sp)
library(ggplot2)

##
## Parse Data
##

setwd("~/Desktop/Applied Bayseian Modeling and Prediction")
CRP_HG_dataset <- read.csv("CRP_HG_dataset_Final.csv")
counties <- unique(CRP_HG_dataset$county)

county_HG21 <- lapply(counties, function(cty) {CRP_HG_dataset$HG_2021[CRP_HG_dataset$county == cty]})
names(county_HG21) <- counties

county_HG22 <- lapply(counties, function(cty) {CRP_HG_dataset$HG_2022[CRP_HG_dataset$county == cty]})
names(county_HG22) <- counties

n <- 30000

##
##  MCMC function for sampling from county probability distribution
##

MCMC <- function(lower, upper, data, n, init){
  p <- matrix(,n,1)
    k <- 1
    p[1] <- init
    for(k in 2:n){
      p.star <- runif(1, lower, upper)
      mh1 <- log(prod(dbinom(data, 1, p.star))*dunif(p.star, lower, upper))
      mh2 <- log(prod(dbinom(data, 1, p[k-1]))*dunif(p[k-1], lower, upper))
      mh <- exp(mh1 - mh2)
      keep <- ifelse(mh > 1, 1, rbinom(1, 1, mh))
      p[k] <- ifelse(keep, p.star, p[k-1])
    }
  return(p);
}

##
##  2021 data analysis
##

post21 <- matrix(,length(county_HG21), n)
burn.in <- 1000;
for(i in 1:length(county_HG21)){
  post21[i,] <- MCMC(0, 1, county_HG21[[i]], n, 0.1)
  hist(post21[1,burn.in:n], freq = FALSE, main = paste(counties[i], "21"))
}

summarypost21 <- matrix(,4,length(counties))
summarypost21[1,] <- counties
summarypost21[2,] <- rowMeans(post21)
for(i in 1:length(counties)) {
summarypost21[c(3, 4),i] <- quantile(post21[i,], probs=c(0.025, 0.975))
}

##
##  Raster stuff below
##

unzip("Bayes_Final_Shapefiles.zip")
sf.ks <- vect("KS_counties.shp")
sf.ks <- project(sf.ks, "EPSG:32614")
plot(sf.ks)
rast(sf.ks)

sf.lpc <- vect("LPC_range_ks.shp")
plot(sf.lpc)
rast(sf.lpc)

sf.kslpc <- rbind(sf.lpc, sf.ks)
plot(sf.kslpc)

r_data21 <- rep(NA, ncell(sf.ks$name))
index <- matrix(,length(sf.ks$name),1)
index <- match(sf.ks$name, summarypost21[1,])
r_data21 <- data.frame(sf.ks$name)
for(i in 1:length(index)) {
  if(is.na(index[i]))
  {
    r_data21[i,2:4] <- NA
  } else {
    r_data21[i,2] <- summarypost21[2, index[i]]
    r_data21[i,3] <- summarypost21[3, index[i]]
    r_data21[i,4] <- summarypost21[4, index[i]]
  }
}
names(r_data21) <- c("names", "expected", "lowerci", "upperci")
r_data21$expected <- as.numeric(r_data21$expected)
r_data21$lowerci <- as.numeric(r_data21$lowerci)
r_data21$upperci <- as.numeric(r_data21$upperci)

countyraster21 <- sf.ks
values(countyraster21) <- r_data21
squareraster21 <- countyraster21[countyraster21$names %in% counties,]
terra::plot(squareraster21, "lowerci",type = "continuous", main="Lowerci 21")
terra::plot(squareraster21, "upperci", type="continuous", main="Upperci 21")
terra::plot(squareraster21, "expected", type="continuous", main="Expected 21")


cropraster21 <- intersect(squareraster21, sf.lpc)
terra::plot(cropraster21, "expected", type="continuous", main="Expected 21")

sf.lpc <- rast(sf.lpc)
terra::plot(smallAreaRast)
plot(squareraster21)

##
##  2022 data analysis
##

post22 <- matrix(,length(county_HG22), n)
burn.in <- 1000;
for(i in 1:length(county_HG22)){
  post22[i,] <- MCMC(0, 1, county_HG22[[i]], n, 0.1)
  hist(post22[i,burn.in:n], freq = FALSE, main = paste(counties[i], "22"))
}

summarypost22 <- matrix(,4,length(counties))
summarypost22[1,] <- counties
summarypost22[2,] <- rowMeans(post22)
for(i in 1:length(counties)) {
  summarypost22[c(3, 4),i] <- quantile(post22[i,], probs=c(0.025, 0.975))
}

##
##  Raster stuff below
##

r_data22 <- rep(NA, ncell(sf.ks$name))
index <- matrix(,length(sf.ks$name),1)
index <- match(sf.ks$name, summarypost22[1,])
r_data22 <- data.frame(sf.ks$name)
for(i in 1:length(index)) {
  if(is.na(index[i]))
  {
    r_data22[i,2:4] <- NA
  } else {
    r_data22[i,2] <- summarypost22[2, index[i]]
    r_data22[i,3] <- summarypost22[3, index[i]]
    r_data22[i,4] <- summarypost22[4, index[i]]
  }
}

names(r_data22) <- c("names", "expected", "lowerci", "upperci")
r_data22$expected <- as.numeric(r_data22$expected)
r_data22$lowerci <- as.numeric(r_data22$lowerci)
r_data22$upperci <- as.numeric(r_data22$upperci)

countyraster22 <- sf.ks
values(countyraster22) <- r_data22
squareraster22 <- countyraster22[countyraster22$names %in% counties,]
terra::plot(squareraster22, "lowerci", type="continuous", main="Lowerci 22")
terra::plot(squareraster22, "upperci", type="continuous", main="Upperci 22")
terra::plot(squareraster22, "expected", type="continuous", main="Expected 22")

cropraster22 <- intersect(squareraster22, sf.lpc)
terra::plot(cropraster22, "expected", , type="continuous", main="Expected 22")

