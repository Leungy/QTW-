library(tidyverse)
library(magrittr)
library(lattice)
library(fields)
library(KernSmooth)
library(ggplot2)
library(class)

setwd("/home/yat-l/Documents/MSDS 7333 QTW/Wk1")

# Read in the raw "offline" text file
txt = readLines("/home/yat-l/Documents/MSDS 7333 QTW/Wk1/offline.final.trace.txt")

# Create a function to parse the data
processLine = function(x){
  tokens = strsplit(x, "[;=,]")[[1]]
  if (length(tokens) == 10) {
    return(NULL)
  }
  tmp = matrix(tokens[ - (1:10) ], ncol=4, byrow = TRUE)
  cbind(matrix(tokens[c(2, 4, 6:8, 10)], nrow(tmp), ncol=6,
               byrow = TRUE), tmp)
}

lines = txt[ substr(txt, 1, 1) != "#" ]
tmp = lapply(lines, processLine)

# Convert the offline data to a data frame
offline = as.data.frame(do.call("rbind", tmp),stringsAsFactors = FALSE)

#head(offline)

# Assign column names to the offline data frame
names(offline) = c("time", "scanMac", "posX", "posY", "posZ","orientation", "mac", "signal","channel", "type")

#column names
numVars = c("time", "posX", "posY", "posZ",
            "orientation", "signal")
#convert the position, signal, and time variables to numeric
offline[ numVars ] = lapply(offline[ numVars ], as.numeric)
offline = offline[ offline$type == "3", ]
offline = offline[ , "type" != names(offline) ]

#time conversion from POSIXt in sec to ms
offline$rawTime = offline$time
offline$time = offline$time/1000
class(offline$time) = c("POSIXt", "POSIXct")

#check variable type
unlist(lapply(offline, class))
summary(offline)

#filter out posZ column
offline = offline[ , !(names(offline) %in% c("scanMac", "posZ"))]
#summary(offline)

# unique orientation
length(unique(offline$orientation))

# #plot empirical CDF of offline orientation
# plot(ecdf(offline$orientation))

#create the rounded angles for ECDF
roundOrientation = function(angles) {
  refs = seq(0, by = 45, length = 9)
  q = sapply(angles, function(o) which.min(abs(o - refs)))
  c(refs[1:8], 0)[q]
}

offline$angle = roundOrientation(offline$orientation)

# #orientation boxplot
# with(offline, boxplot(orientation ~ angle,
#                       xlab = "nearest 45 degree angle",
#                       ylab="orientation"))

#unique addresses and channels
c(length(unique(offline$mac)), length(unique(offline$orientation)))
table(offline$mac)

#filter macs
subMacs = names(sort(table(offline$mac), decreasing = TRUE))[1:7]
offline = offline[ offline$mac %in% subMacs, ]

#count remaining MACxChannel
macChannel = with(offline, table(mac, channel))
apply(macChannel, 1, function(x) sum(x > 0))

#delete channel from offline
offline = offline[ , "channel" != names(offline)]

#combining X,Y together
locDF = with(offline, by(offline, list(posX, posY), function(x) x))
length(locDF)
#filter all n/a from locDF
sum(sapply(locDF, is.null))
locDF = locDF[ !sapply(locDF, is.null) ]
length(locDF)

#determine the number of observations recorded at each location
locCounts = sapply(locDF, nrow)

#keep the position information with the location
locCounts = sapply(locDF,function(df) c(df[1, c("posX", "posY")], count = nrow(df)))
#locCounts type and dimensions
class(locCounts)
dim(locCounts)

locCounts[ , 1:8]
# transpose location to column and build visualization
locCounts = t(locCounts)
plot(locCounts, type = "n", xlab = "", ylab = "")
text(locCounts, labels = locCounts[,3], cex = .8, srt = 45)

# 5-pt summary for offline signal
summary(offline$signal)

# density plot compare the distributions of signal strength for different angles
# and MAC addresses at the central location of x = 23 and y = 4
densityplot( ~ signal | mac + factor(angle), data = offline,
             subset = posX == 24 & posY == 4 &
               mac != "00:0f:a3:39:dd:cd",
             bw = 0.5, plot.points = FALSE)

#Unique PosX,PosY (PosXY) combinations
offline$posXY = paste(offline$posX, offline$posY, sep = "-")

#create list for each PosXY combination, angle and access pt
byLocAngleAP = with(offline,
                    by(offline, list(posXY, angle, mac),
                       function(x) x))
#summary stat for each PosXY combination, angle and access pt
signalSummary =
  lapply(byLocAngleAP,
         function(oneLoc) {
           ans = oneLoc[1, ]
           ans$medSignal = median(oneLoc$signal)
           ans$avgSignal = mean(oneLoc$signal)
           ans$num = length(oneLoc$signal)
           ans$sdSignal = sd(oneLoc$signal)
           ans$iqrSignal = IQR(oneLoc$signal)
           ans
         })
offlineSummary = do.call("rbind", signalSummary)

# #Std Dev Boxplot by subgroup signal mean
# breaks = seq(-90, -30, by = 5)
# bwplot(sdSignal ~ cut(avgSignal, breaks = breaks),
#        data = offlineSummary,
#        subset = mac != "00:0f:a3:39:dd:cd",
#        xlab = "Mean Signal", ylab = "SD Signal")

# #smooth scatter plot; Comparison of Mean and Median Signal Strength
# with(offlineSummary,
#      smoothScatter((avgSignal - medSignal) ~ num, xlab = "Number of Observations",
#                    ylab = "mean - median"))
# abline(h = 0, col = "#984ea3", lwd = 2)
# 
# #predicting the difference for each value of num
# lo.obj =  with(offlineSummary, loess(diff ~ num, data = data.frame(diff = (avgSignal - medSignal),num = num)))
# lo.obj.pr = predict(lo.obj, newdata = data.frame(num = (70:120)))
# lines(x = 70:120, y = lo.obj.pr, col = "#4daf4a", lwd = 2)

#filter to get the fith MAC and the Angle from offlineSummary df
oneAPAngle = subset(offlineSummary, mac == subMacs[5] & angle == 0)
smoothSS = Tps(oneAPAngle[, c("posX","posY")], oneAPAngle$avgSignal)
vizSmooth = predictSurface(smoothSS)
plot.surface(vizSmooth, type = "C")
points(oneAPAngle$posX, oneAPAngle$posY, pch=19, cex = 0.5)

#median signal maps for 2 AP+ 2Angles at 0 & 135
surfaceSS = function(data,mac=mac,angle=0){
  oneAPAngle = data[data$mac == mac & data$angle == angle,] 
  smoothSS = Tps(oneAPAngle[, c("posX","posY")], oneAPAngle$avgSignal)
  vizSmooth = predictSurface(smoothSS)
  plot.surface(vizSmooth, type = "C")
  points(oneAPAngle$posX, oneAPAngle$posY, pch=19, cex = 0.5)
}

# parCur = par(mfrow = c(2,2), mar = rep(1, 4))
# mapply(surfaceSS, mac = subMacs[ rep(c(5, 1), each = 2) ],
#        angle = rep(c(0, 135), 2),
#        data = list(data = offlineSummary))
# par(parCur)

# subset from offlineSummary to return all subMacs that are not the second subMac
offlineSummary = subset(offlineSummary, mac != subMacs[2])

#matrix with relevant positions for 6 access points
AP = matrix( c( 7.5, 6.3, 2.5, -.8, 12.8, -2.8,
                1, 14, 33.5, 9.3, 33.5, 2.8),
             ncol = 2, byrow = TRUE,
             dimnames = list(subMacs[ -2 ], c("x", "y") ))

# X,Y position differences between emitter(device) and receiver(AP)
diffs = offlineSummary[ , c("posX", "posY")] - AP[ offlineSummary$mac, ]

#euclidean distance between emitter and receiver
offlineSummary$dist = sqrt(diffs[ , 1]^2 + diffs[ , 2]^2)

#scatter plots for each AO and orientation
xyplot(signal ~ dist | factor(mac) + factor(angle),
       data = offlineSummary, pch = 19, cex = 0.3,
       xlab ="distance")

#readData function
readData = 
  function(filename = "/home/yat-l/Documents/MSDS 7333 QTW/Wk1/offline.final.trace.txt", 
           subMacs = c("00:0f:a3:39:e1:c0", "00:0f:a3:39:dd:cd", "00:14:bf:b1:97:8a",
                       "00:14:bf:3b:c7:c6", "00:14:bf:b1:97:90", "00:14:bf:b1:97:8d",
                       "00:14:bf:b1:97:81"))
  {
    txt = readLines(filename)
    lines = txt[ substr(txt, 1, 1) != "#" ]
    tmp = lapply(lines, processLine)
    offline = as.data.frame(do.call("rbind", tmp), stringsAsFactors= FALSE) 
    
    names(offline) = c("time", "scanMac", 
                       "posX", "posY", "posZ", "orientation", 
                       "mac", "signal", "channel", "type")
    
    # subset all type 3
    offline = offline[ offline$type == "3", ]
    
    # delete scanMac, posZ, channel, and type
    dropVars = c("scanMac", "posZ", "channel", "type")
    offline = offline[ , !( names(offline) %in% dropVars ) ]
    
    # subset out all macs that's not in subMacs
    offline = offline[ offline$mac %in% subMacs, ]
    
    # convert numeric values
    numVars = c("time", "posX", "posY", "orientation", "signal")
    offline[ numVars ] = lapply(offline[ numVars ], as.numeric)
    
    # convert time to POSIX
    offline$rawTime = offline$time
    offline$time = offline$time/1000
    class(offline$time) = c("POSIXt", "POSIXct")
    
    # round orientations to nearest 45
    offline$angle = roundOrientation(offline$orientation)
    
    # Combine X and Y positions into a X-Y column
    offline$posXY = paste(offline$posX, offline$posY, sep = "-")
    
    return(offline)
  }

#read online data 
macs = unique(offlineSummary$mac)
online = readData("/home/yat-l/Documents/MSDS 7333 QTW/Wk1/online.final.trace.txt", subMacs = macs)

# Tallying number of signal strengths at each location
tabonlineXYA = table(online$posXY, online$angle)
tabonlineXYA[1:6, ]

# reorganize signal strength into columns
keepVars = c("posXY", "posX","posY", "orientation", "angle")
byLoc = with(online,
             by(online, list(posXY),
                function(x) {
                  ans = x[1, keepVars]
                  avgSS = tapply(x$signal, x$mac, mean)
                  y = matrix(avgSS, nrow = 1, ncol = 6,
                             dimnames = list(ans$posXY, names(avgSS)))
                  cbind(ans, y)
                }))
# dimension [1] 60 11
onlineSummary = do.call("rbind", byLoc)
#names(onlineSummary)

# number of angle orientations to be included, angle at 45 deg increment 
# m = number of angles; angleNewObs = angle of the new observation
m = 3; angleNewObs = 230

refs = seq(0, by = 45, length = 8)
nearestAngle = roundOrientation(angleNewObs)
if (m %% 2 == 1) {
  angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
} else {
  m = m + 1
  angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
  if (sign(angleNewObs - nearestAngle) > -1)
    angles = angles[ -1 ]
  else
    angles = angles[ -m ]
}

angles = angles + nearestAngle
angles[angles < 0] = angles[ angles < 0 ] + 360
angles[angles > 360] = angles[ angles > 360 ] - 360

offlineSubset =  offlineSummary[ offlineSummary$angle %in% angles, ]

# aggregate signal strengths from angles and reshape data with reshapeSS()
reshapeSS = function(data, varSignal = "signal", 
                     keepVars = c("posXY", "posX","posY")) {
  byLocation =
    with(data, by(data, list(posXY), 
                  function(x) {
                    ans = x[1, keepVars]
                    avgSS = tapply(x[ , varSignal ], x$mac, mean)
                    y = matrix(avgSS, nrow = 1, ncol = 6,
                               dimnames = list(ans$posXY,
                                               names(avgSS)))
                    cbind(ans, y)
                  }))
  
  newDataSS = do.call("rbind", byLocation)
  return(newDataSS)
}

# build training signal strength from offlineSubset
trainSS = reshapeSS(offlineSubset, varSignal = "avgSignal")

# selectTrain function
# angleNewObs= the angle of the new observation; 
# signals = the training data, i.e., data in the format of offlineSummary;
# m = the number of angles to include from signals. can be 1 to 3 (3 is the total orientation angles they use)
selectTrain = function(angleNewObs, signals = NULL, m = 1){
    refs = seq(0, by = 45, length = 8)
    nearestAngle = roundOrientation(angleNewObs)
    if (m %% 2 == 1) {
      angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
    } else {
      m = m + 1
      angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
      if (sign(angleNewObs - nearestAngle) > -1)
        angles = angles[ -1 ]
      else
        angles = angles[ -m ]
    }
    
    angles = angles + nearestAngle
    angles[angles < 0] = angles[ angles < 0 ] + 360
    angles[angles > 360] = angles[ angles > 360 ] - 360
    
    offlineSubset =  offlineSummary[signals$angle %in% angles, ]  
    reshapeSS(offlineSubset, varSignal = "avgSignal")  
}

# # test angle observed at 130 using offlineSummary data with number of angles to include from signals column m=3
# train130 = selectTrain(130, offlineSummary, m = 3)

#head(train130)
#length(train130[[1]]) #166 observations

#find nearest neighbors, calculate distance of new source point to all observations in train set
findNN = function(newSignal, trainSubset) {
  diffs = apply(trainSubset[ , 4:9], 1,
                function(x) x - newSignal)
  dists = apply(diffs, 2, function(x) sqrt(sum(x^2)) )
  closest = order(dists)
  return(trainSubset[closest, 1:3 ])
}

# #estimate the XY pos, closeXY from findNN
# estXY = lapply(closeXY, function(x) sapply(x, function(x) mean(x[1:k])))

# predict XY position, default K neighbours =3 but can be overloaded to other integer values
predXY = function(newSignals, newAngles, trainData, numAngles = 1, k = 3){
  closeXY = list(length = nrow(newSignals))
  for (i in 1:nrow(newSignals)) {
    trainSS = selectTrain(newAngles[i], trainData, m = numAngles)
    closeXY[[i]] = findNN(newSignal = as.numeric(newSignals[i, ]), trainSS)
  }
  estXY = lapply(closeXY,function(x) sapply(x[ , 2:3], function(x) mean(x[1:k])))
  estXY = do.call("rbind", estXY)
  return(estXY)
}

#sum of squared errors
calcError = function(estXY, actualXY)  sum( rowSums( (estXY - actualXY)^2) )

# estXYk1 = predXY(newSignals = onlineSummary[ , 6:11],
#                  newAngles = onlineSummary[ , 4],
#                  offlineSummary, numAngles = 3, k = 1)
# estXYk3 = predXY(newSignals = onlineSummary[ , 6:11],
#                  newAngles = onlineSummary[ , 4],
#                  offlineSummary, numAngles = 3, k = 3)
# actualXY = onlineSummary[ , c("posX", "posY")]
# #compare 1-NN to 3-NN
# sapply(list(estXYk1, estXYk3), calcError, actualXY)


#k-fold Cross-Validation
v = 5
permuteLocs = sample(unique(offlineSummary$posXY))
permuteLocs = matrix(permuteLocs, ncol = v, nrow = floor(length(permuteLocs)/v))
onlineFold = subset(offlineSummary, posXY %in% permuteLocs[ , 1])


reshapeSS = function(data, varSignal = "signal", 
                     keepVars = c("posXY", "posX","posY"),
                     sampleAngle = FALSE, 
                     refs = seq(0, 315, by = 45)) {
  byLocation =  with(data, by(data, list(posXY), 
                  function(x) {
                    if (sampleAngle) {x = x[x$angle == sample(refs, size = 1), ]}
                    ans = x[1, keepVars]
                    avgSS = tapply(x[ , varSignal ], x$mac, mean)
                    y = matrix(avgSS, nrow = 1, ncol = 6, dimnames = list(ans$posXY,names(avgSS)))
                    cbind(ans, y)
                  }))
  
  newDataSS = do.call("rbind", byLocation)
  return(newDataSS)
}

# Exclude mac = 00:0f:a3:39:dd:cd ; error=259.4
offline1 = offline[offline$mac != "00:0f:a3:39:dd:cd",]
# variables to keep, building online CV set summary from offline
keepVars = c("posXY", "posX","posY", "orientation", "angle")
onlineCVSummary1 = reshapeSS(offline1, keepVars = keepVars, sampleAngle = TRUE)
# first online fold
onlineFold1 = subset(onlineCVSummary1,posXY %in% permuteLocs[ , 1])
# first offline fold
offlineFold1 = subset(offlineSummary,posXY %in% permuteLocs[ , -1])
# Use CV to estimate XY of online from offline
estFold1 = predXY(newSignals = onlineFold1[ , 6:11],
                 newAngles = onlineFold1[ , 4],
                 offlineFold1, numAngles = 3, k = k)

#estimate error in actual using predicted
actualFold1 = onlineFold1[ , c("posX", "posY")]
calcError(estFold1, actualFold1)

# Exclude mac = 00:0f:a3:39:e1:c0 ; error= 492.1
offline2 = offline[offline$mac != "00:0f:a3:39:e1:c0",]
# variables to keep, building online CV set summary from offline
keepVars = c("posXY", "posX","posY", "orientation", "angle")
onlineCVSummary2 = reshapeSS(offline2, keepVars = keepVars, sampleAngle = TRUE)
# first online fold
onlineFold2 = subset(onlineCVSummary2,posXY %in% permuteLocs[ , 1])
# first offline fold
offlineFold2 = subset(offlineSummary,posXY %in% permuteLocs[ , -1])
# Use CV to estimate XY of online from offline
estFold2 = predXY(newSignals = onlineFold2[ , 6:11],
                 newAngles = onlineFold2[ , 4],
                 offlineFold2, numAngles = 3, k = k)

#estimate error in actual using predicted
actualFold2 = onlineFold2[ , c("posX", "posY")]
calcError(estFold2, actualFold2)


# Approximate K-neighbors for KNN where K=1 to 20 aggregate over k-fold CV
K = 15

# SSE calculation. Masking one of the 00:0f:a3:39 MAC
err1 = rep(0, K)
for (j in 1:v) {
  onlineFold1 = subset(onlineCVSummary1,posXY %in% permuteLocs[ , j])
  offlineFold1 = subset(offlineSummary,posXY %in% permuteLocs[ , -j])
  actualFold1 = onlineFold1[ , c("posX", "posY")]
  for (k in 1:K) {
    estFold1 = predXY(newSignals = onlineFold1[ , 6:11],
                     newAngles = onlineFold1[ , 4],
                     offlineFold1, numAngles = 3, k = k)
    err1[k] = err1[k] + calcError(estFold1, actualFold1)
  }
}

# SSE calculation. Masking 00:0f:a3:39:e1:c0
err2 = rep(0, K)
for (j in 1:v) {
  onlineFold2 = subset(onlineCVSummary2,posXY %in% permuteLocs[ , j])
  offlineFold2 = subset(offlineSummary,posXY %in% permuteLocs[ , -j])
  actualFold2 = onlineFold2[ , c("posX", "posY")]
  for (k in 1:K) {
    estFold2 = predXY(newSignals = onlineFold2[ , 6:11],
                     newAngles = onlineFold2[ , 4],
                     offlineFold2, numAngles = 3, k = k)
    err2[k] = err2[k] + calcError(estFold2, actualFold2)
  }
}

# elbow plot
#SSE 1160 , KNN=7
plot(y = err1, x = (1:K),  type = "l", lwd= 2,
     ylim = c(900, 2200),
     xlab = "Number of Neighbors",
     ylab = "Sum of Square Errors")
title(main ="KNN Elbow Plot: 5-Fold Cross-Validated SSE values for K = 1 to 15", 
      sub="Masking MAC == 00:0f:a3:39:dd:cd ")
rmseMin = min(err1)
kMin = which(err1 == rmseMin)[1]
segments(x0 = 0, x1 = kMin, y0 = rmseMin, col = gray(0.4), lty = 2, lwd = 2)
segments(x0 = kMin, x1 = kMin, y0 = 1100,  y1 = rmseMin, col = grey(0.4), lty = 2, lwd = 2)
mtext(kMin, side = 1, line = 1, at = kMin, col = grey(0.4))
text(x = kMin - 2, y = rmseMin + 40, label = as.character(round(rmseMin)), col = grey(0.4))

# SSE after fixing line 433 from offline1 to offline2, 2710, KNN of 10
plot(y = err2, x = (1:K),  type = "l", lwd= 2,
     ylim = c(2000, 3900),
     xlab = "Number of Neighbors",
     ylab = "Sum of Square Errors")
title(main ="KNN Elbow Plot: 5-Fold Cross-Validated SSE values for K = 1 to 15", 
      sub="Masking MAC == 00:0f:a3:39:e1:c0 ")
rmseMin = min(err2)
kMin = which(err2 == rmseMin)[1]
segments(x0 = 0, x1 = kMin, y0 = rmseMin, col = gray(0.4), lty = 2, lwd = 2)
segments(x0 = kMin, x1 = kMin, y0 = 1100,  y1 = rmseMin, col = grey(0.4), lty = 2, lwd = 2)
mtext(kMin, side = 1, line = 1, at = kMin, col = grey(0.4))
text(x = kMin - 2, y = rmseMin + 40, label = as.character(round(rmseMin)), col = grey(0.4))

# KNN of 7 is the optimal solution regardless of which MAC to mask
# Calculated Actual-Predicted error. SSE= 273.7
estXYk7 = predXY(newSignals = onlineSummary[ , 6:11], 
                 newAngles = onlineSummary[ , 4], 
                 offlineSummary, numAngles = 3, k = 7)
actualXY = onlineSummary[ , c("posX", "posY")]
calcError(estXYk7, actualXY)

# find weighted dist average? 
# update the findNN function to return the distance (book p35)
# how close they actually are from the new observation’s signals
findNN = function(newSignal, trainSubset) {
  diffs = apply(trainSubset[ , 4:9], 1,
                function(x) x - newSignal)
  dists = apply(diffs, 2, function(x) sqrt(sum(x^2)) )
  closest = order(dists)
  return(list(trainSubset[closest, 1:3 ], dists[order(dists)]))
}

#find the weights (1/di) / summation of {i from 1 to k (1/di)}
# modify preXY function using the weighted distance
predXY = function(newSignals, newAngles, trainData, numAngles = 1, k = 3){
  closeXY = list(length = nrow(newSignals))
  closeDist = list(length = nrow(newSignals))
#findNN returns 2 variables, first one is XY and second on is the distance
    for (i in 1:nrow(newSignals)) {
    trainSS = selectTrain(newAngles[i], trainData, m = numAngles)
    results = findNN(newSignal = as.numeric(newSignals[i, ]), trainSS)
    closeXY[[i]] = results[[1]]
    closeDist[[i]] = results[[2]]
  }
#weighted distance
  # weightedDist returns a list of closeDist
  weightedDist = list(length = length(closeDist))
# find all weighted distance in the closeDist list
  for(i in 1:length(closeDist)){
    # k as in number of neighbours to use
    W = list(length(k))
    # for each j between 1 to k
    for (j in 1:k){
      # find weighted distance equation
      W[j] = (1/closeDist[[i]][j])/sum(1/closeDist[[i]][1:k])    
    }
    # return the weighted distance values to the weightedDist list
    weightedDist[[i]] = W
  }
# apply weighted distance to refine the estimated XY
  refinedXY = list(length = length(closeXY))
  # multiply weight to each closeXY
  for (i in 1:length(closeXY)){
    # not sure if this is right
    refinedXY[[i]]=as.matrix(closeXY[[i]][1:k,2:3])*unlist(weightedDist[[i]])
  }
  estXY = lapply(refinedXY,function(x) apply(x, 2, function(x)sum(x)))
  estXY = do.call("rbind", estXY)
  return(estXY)
}

# re-run the KNN masking 00:0f:a3:39:dd:cd
# Approximate K-neighbors for KNN where K=1 to 20 aggregate over k-fold CV
K = 15

# SSE calculation. Masking 00:0f:a3:39:dd:cd
err = rep(0, K)
for (j in 1:v) {
  onlineFold1 = subset(onlineCVSummary1,posXY %in% permuteLocs[ , j])
  offlineFold1 = subset(offlineSummary,posXY %in% permuteLocs[ , -j])
  actualFold1 = onlineFold1[ , c("posX", "posY")]
  for (k in 1:K) {
    estFold1 = predXY(newSignals = onlineFold1[ , 6:11],
                      newAngles = onlineFold1[ , 4],
                      offlineFold1, numAngles = 3, k = k)
    err[k] = err[k] + calcError(estFold1, actualFold1)
  }
}

plot(y = err, x = (1:K),  type = "l", lwd= 2,
     ylim = c(900, 2200),
     xlab = "Number of Neighbors",
     ylab = "Sum of Square Errors")
title(main ="KNN Elbow Plot: 5-Fold Cross-Validated SSE values for K = 1 to 15", 
      sub="Masking MAC == 00:0f:a3:39:dd:cd ; Modified findNN and predXY functions")
rmseMin = min(err)
kMin = which(err == rmseMin)[1]
segments(x0 = 0, x1 = kMin, y0 = rmseMin, col = gray(0.4), lty = 2, lwd = 2)
segments(x0 = kMin, x1 = kMin, y0 = 1100,  y1 = rmseMin, col = grey(0.4), lty = 2, lwd = 2)
mtext(kMin, side = 1, line = 1, at = kMin, col = grey(0.4))
text(x = kMin - 2, y = rmseMin + 40, label = as.character(round(rmseMin)), col = grey(0.4))

# SSE with weighted distance
# Calculated Actual-Predicted error. SSE= 267.4
estXYk8 = predXY(newSignals = onlineSummary[ , 6:11], 
                 newAngles = onlineSummary[ , 4], 
                 offlineSummary, numAngles = 3, k = kMin )
actualXY = onlineSummary[ , c("posX", "posY")]
calcError(estXYk8, actualXY)