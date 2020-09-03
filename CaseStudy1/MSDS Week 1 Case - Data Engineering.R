
setwd("/home/yat-l/Documents/MSDS 7333 QTW/Wk1")

# Read in the raw "offline" text file
txt = readLines("/home/yat-l/Documents/MSDS 7333 QTW/Wk1/offline.final.trace.txt")

# Create a function to parse the data
processLine = function(x)
{
  tokens = strsplit(x, "[;=,]")[[1]]
  if (length(tokens) == 10) {
    return(NULL)
  }
  tmp = matrix(tokens[ - (1:10) ], , 4, byrow = TRUE)
  cbind(matrix(tokens[c(2, 4, 6:8, 10)], nrow(tmp), 6,
               byrow = TRUE), tmp)
}

lines = txt[ substr(txt, 1, 1) != "#" ]
tmp = lapply(lines, processLine)

# Convert the offline data to a data frame
offline = as.data.frame(do.call("rbind", tmp),stringsAsFactors = FALSE)

head(offline)

# Assign column names to the offline data frame
names(offline) = c("time", "scanMac", "posX", "posY", "posZ","orientation", "mac", "signal","channel", "type")

#column names
numVars = c("time", "posX", "posY", "posZ",
            "orientation", "signal")
#convert the position, signal, and time variables to numeric
offline[ numVars ] = lapply(offline[ numVars ], as.numeric)
offline = offline[ offline$type == "3", ]
offline = offline[ , "type" != names(offline) ]

#dimension of the offline df
dim(offline)

#time conversion from POSIXt in sec to ms
offline$rawTime = offline$time
offline$time = offline$time/1000
class(offline$time) = c("POSIXt", "POSIXct")

#check variable type
unlist(lapply(offline, class))
summary(offline[, numVars])

#filter out posZ column
offline = offline[ , !(names(offline) %in% c("scanMac", "posZ"))]
summary(offline)

# unique orientation
length(unique(offline$orientation))

#plot empirical CDF of offline orientation
plot(ecdf(offline$orientation))

#create the rounded angles for ECDF
roundOrientation = function(angles) {
  refs = seq(0, by = 45, length = 9)
  q = sapply(angles, function(o) which.min(abs(o - refs)))
  c(refs[1:8], 0)[q]
}

offline$angle = roundOrientation(offline$orientation)

#orientation boxplot
with(offline, boxplot(orientation ~ angle,
                      xlab = "nearest 45 degree angle",
                      ylab="orientation"))

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
library(lattice)
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

#Std Dev Boxplot by subgroup signal mean
breaks = seq(-90, -30, by = 5)
bwplot(sdSignal ~ cut(avgSignal, breaks = breaks),
       data = offlineSummary,
       subset = mac != "00:0f:a3:39:dd:cd",
       xlab = "Mean Signal", ylab = "SD Signal")

#smooth scatter plot; Comparison of Mean and Median Signal Strength
library(KernSmooth)
with(offlineSummary,
     smoothScatter((avgSignal - medSignal) ~ num, xlab = "Number of Observations",
                   ylab = "mean - median"))
abline(h = 0, col = "#984ea3", lwd = 2)

#predicting the difference for each value of num
lo.obj =  with(offlineSummary, loess(diff ~ num, data = data.frame(diff = (avgSignal - medSignal),num = num)))
lo.obj.pr = predict(lo.obj, newdata = data.frame(num = (70:120)))
lines(x = 70:120, y = lo.obj.pr, col = "#4daf4a", lwd = 2)

#subset a MAC+ Angle from offlineSummary df
library(fields)
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

parCur = par(mfrow = c(2,2), mar = rep(1, 4))
mapply(surfaceSS, mac = subMacs[ rep(c(5, 1), each = 2) ],
       angle = rep(c(0, 135), 2),
       data = list(data = offlineSummary))
par(parCur)

# subset from offlineSummary to return all subMacs that are not the second subMac
offlineSummary = subset(offlineSummary, mac != subMacs[2])

#matrix with relevant positions for 6 access points
AP = matrix( c( 7.5, 6.3, 2.5, -.8, 12.8, -2.8,
                1, 14, 33.5, 9.3, 33.5, 2.8),
             ncol = 2, byrow = TRUE,
             dimnames = list(subMacs[ -2 ], c("x", "y") ))

AP

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
onlineSummary = do.call("rbind", byLoc)

#
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

# Observe one observation from the offline data frame
offline[offline$time==1139643118358,]

library(tidyverse)
library(magrittr)

# View the data frame exclusive of the scanMac, channel, and type columns
select(offline,-c(scanMac,channel,type))

# See if any observations have a posZ-value not equal to 0.0
offline[offline['posZ']!='0.0',]

# Get a full list of mac values and their types (we are interested in type == 3)
macTypeDF <- unique(select(offline, c(mac,type)))
macTypeDF[order(macTypeDF$type),]

# Filter mac values where type matches 3
type_3 <- filter(macTypeDF, type=='3')

# List the unique mac values
vals = data.frame(table(offline['mac']))
vals[order(-vals$Freq),]

valsUpdated <- vals %>% inner_join(macTypeDF, by = c("Var1" = "mac"))
valsUpdated <- valsUpdated[order(valsUpdated$type, -valsUpdated$Freq),]

# round all signals from decimals to integers
offline$signal %<>% as.integer

#filter all offline where type is not 1
offline <- offline[offline$type != 1, ]

#7 mac values to keep
keepMacs <- c('00:0f:a3:39:e1:c0',
              '00:0f:a3:39:dd:cd',
              '00:14:bf:b1:97:8a',
              '00:14:bf:3b:c7:c6',
              '00:14:bf:b1:97:90',
              '00:14:bf:b1:97:8d',
              '00:14:bf:b1:97:81'
              )
#filter offline df where the mac values match the values in keepMacs
offline <- offline[offline$mac %in% keepMacs ,]


# Pivot the data frame (or cast it; make it wider) but putting the mac values and their associated signals in the columns
offlineOut<-select(offline, -c(channel,scanMac,type)) %>% pivot_wider(names_from = mac,values_from = signal, values_fn = list(signal=mean))

#what is stored in nas column? how do you calculate that?
offlineOut$nas<-rowSums(is.na(offlineOut))

# View the final data frame
offlineOut

# Process the online data
onlineTxt = readLines("/home/yat-l/Documents/MSDS 7333 QTW/Wk1/online.final.trace.txt")

onlineLines = onlineTxt[ substr(onlineTxt, 1, 1) != "#" ]
onlineTmp = lapply(onlineLines, processLine)

# Convert the online data to a data frame
online = as.data.frame(do.call("rbind", onlineTmp),stringsAsFactors = FALSE)

head(online)

# there are 61 unique orientations in online df, 203 in offlineOut 
names(online) = c("time", "scanMac", "posX", "posY", "posZ",
                   "orientation", "mac", "signal",
                   "channel", "type")

head(online)

# View the data frame exclusive of the scanMac, channel, and type columns
select(online,-c(scanMac,channel,type))

# See if any observations have a posZ-value not equal to 0.0
online[online['posZ']!='0.0',]

# List the unique mac values and their frequencies in the offline data frame
onlineVals = data.frame(table(online['mac']))
onlineVals[order(-onlineVals$Freq),]

online$signal %<>% as.integer


# keepMacs <- c('00:0f:a3:39:e1:c0',
#               '00:0f:a3:39:dd:cd',
#               '00:14:bf:b1:97:8a',
#               '00:14:bf:3b:c7:c6',
#               '00:14:bf:b1:97:90',
#               '00:14:bf:b1:97:8d',
#               '00:14:bf:b1:97:81'
# )

#are we still using the keepMacs from earlier?
online <- online[online$mac %in% keepMacs ,]

# Pivot the data frame (or cast it; make it wider) but putting the mac values and their associated signals in the columns
onlineOut<-select(online, -c(channel,scanMac,type)) %>% pivot_wider(names_from = mac,values_from = signal, values_fn = list(signal=mean))

onlineOut$nas<-rowSums(is.na(onlineOut))

# View the final data frame
onlineOut


offlineOut$posXY <- paste(offlineOut$posX, offlineOut$posY, sep="-")
length(unique(offlineOut$posXY))

onlineOut$posXY <- paste(onlineOut$posX, onlineOut$posY, sep="-")
length(unique(onlineOut$posXY))







