
######## Basic setting ########

pathsnubh = "snubh_result"
omrsnubh = readRDS(file.path(pathsnubh, "outcomeModelReference.rds"))
mtsnubh <- readRDS(file.path(pathsnubh, omrsnubh$strataFile[omrsnubh$analysisId == 6 & omrsnubh$targetId == 2023 & omrsnubh$outcomeId == 2022]))

pathsnuh = "snuh_result"
omrsnuh = readRDS(file.path(pathsnuh, "outcomeModelReference.rds"))
mtsnuh <- readRDS(file.path(pathsnuh, omrsnuh$strataFile[omrsnuh$analysisId == 6 & omrsnuh$targetId == 182 & omrsnuh$outcomeId == 186]))

pathcmcs = "cmcs_result"
omrcmcs = readRDS(file.path(pathcmcs, "outcomeModelReference.rds"))
mtcmcs <- readRDS(file.path(pathcmcs, omrcmcs$strataFile[omrcmcs$analysisId == 6 & omrcmcs$targetId == 12 & omrcmcs$outcomeId == 10]))



#################### Differential Privacy for Kaplan-Meier curve ####################

##### SNUBH #####

snubh_comp = mtsnubh[mtsnubh$treatment == 0, ]
snubh_tar = mtsnubh[mtsnubh$treatment == 1, ]

snubh_timemax = max(mtsnubh$survivalTime)

### Comparator

snubh_comp_time = unique(sort(snubh_comp$survivalTime[which(snubh_comp$outcomeCount != 0)]))
snubh_comp_tk = length(snubh_comp_time)

snubh_comp_M = matrix(rep(NA, snubh_comp_tk*3), nrow = snubh_comp_tk, ncol = 3)

for (i in 1:snubh_comp_tk) {
  snubh_comp_M[i,1] = sum(snubh_comp$survivalTime >= snubh_comp_time[i])
  snubh_comp_M[i,2] = length(which(snubh_comp$survivalTime[which(snubh_comp$outcomeCount != 0)] == snubh_comp_time[i]))
  if(i >= 2) {
    snubh_comp_M[i-1,3] = snubh_comp_M[i-1,1] - snubh_comp_M[i,1] - snubh_comp_M[i-1,2]
  }
}

## Survival function
snubh_comp_S = rep(NA, snubh_comp_tk) 

for (i in snubh_comp_time) {
  j = which(snubh_comp_time == i)
  if (j == 1) {
    snubh_comp_S[j] = (snubh_comp_M[1, 1] - snubh_comp_M[1, 2]) / snubh_comp_M[1, 1]
  } else {
    snubh_comp_S[j] = snubh_comp_S[j-1] * (snubh_comp_M[j, 1] - snubh_comp_M[j, 2]) / snubh_comp_M[j, 1]
  }
}

snubh_comp_x = seq(1, snubh_comp_time[1])
snubh_comp_y = rep(1, snubh_comp_time[1])

for (i in 1:(snubh_comp_tk-1)) {
  snubh_comp_x = c(snubh_comp_x, seq(snubh_comp_time[i], snubh_comp_time[i+1]))
  snubh_comp_y = c(snubh_comp_y, rep(snubh_comp_S[i], snubh_comp_time[i+1]-snubh_comp_time[i]+1))
}

snubh_comp_x = c(snubh_comp_x, snubh_comp_time[snubh_comp_tk])
snubh_comp_y = c(snubh_comp_y, snubh_comp_S[snubh_comp_tk])


### Target

snubh_tar_time = unique(sort(snubh_tar$survivalTime[which(snubh_tar$outcomeCount != 0)]))
snubh_tar_tk = length(snubh_tar_time)

snubh_tar_M = matrix(rep(NA, snubh_tar_tk*3), nrow = snubh_tar_tk, ncol = 3)

for (i in 1:snubh_tar_tk) {
  snubh_tar_M[i,1] = sum(snubh_tar$survivalTime >= snubh_tar_time[i])
  snubh_tar_M[i,2] = length(which(snubh_tar$survivalTime[which(snubh_tar$outcomeCount != 0)] == snubh_tar_time[i]))
  if(i >= 2) {
    snubh_tar_M[i-1,3] = snubh_tar_M[i-1,1] - snubh_tar_M[i,1] - snubh_tar_M[i-1,2]
  }
}

snubh_tar_S = rep(NA, snubh_tar_tk)

for (i in snubh_tar_time) {
  j = which(snubh_tar_time == i)
  if (j == 1) {
    snubh_tar_S[j] = (snubh_tar_M[1, 1] - snubh_tar_M[1, 2]) / snubh_tar_M[1, 1]
  } else {
    snubh_tar_S[j] = snubh_tar_S[j-1] * (snubh_tar_M[j, 1] - snubh_tar_M[j, 2]) / snubh_tar_M[j, 1]
  }
}

snubh_tar_x = seq(1, snubh_tar_time[1])
snubh_tar_y = rep(1, snubh_tar_time[1])

for (i in 1:(snubh_tar_tk-1)) {
  snubh_tar_x = c(snubh_tar_x, seq(snubh_tar_time[i], snubh_tar_time[i+1]))
  snubh_tar_y = c(snubh_tar_y, rep(snubh_tar_S[i], snubh_tar_time[i+1]-snubh_tar_time[i]+1))
}

snubh_tar_x = c(snubh_tar_x, snubh_tar_time[snubh_tar_tk])
snubh_tar_y = c(snubh_tar_y, snubh_tar_S[snubh_tar_tk])

### Plot

plot(snubh_comp_x, snubh_comp_y, type = 'l', xlim = c(1, snubh_timemax), 
     ylim = c(min(snubh_comp_y, snubh_tar_y) - 0.01, 1), lwd = 3, main = "SNUBH(Original)")
points(snubh_tar_x, snubh_tar_y, type = 'l', lwd = 3, col ='red')
legend("topright", legend=c("Comparator", "Target"), col = c('black', 'red'), lwd = 3)

################ DP

set.seed(1)
S = 2 ## L1 Sensitivity
eps = c(1, 2, 3)

### eps = 1

### DP-comp(1)

snubh_comp_Mp_1 = snubh_comp_M

snubh_comp_noise_1 = rexp(nrow(snubh_comp_M)*2+1, eps[1]/S) * sign(runif(nrow(snubh_comp_M)*2+1, -1, 1))
snubh_comp_noise_1 = pmax(snubh_comp_noise_1, rep(0, nrow(snubh_comp_M)*2+1))

snubh_comp_Mp_1[1,1] = snubh_comp_Mp_1[1,1] + snubh_comp_noise_1[1]
snubh_comp_Mp_1[,2] = snubh_comp_M[,2] + snubh_comp_noise_1[2:(nrow(snubh_comp_M)+1)]
## snubh_comp_Mp_1[,3] = snubh_comp_M[,3] + snubh_comp_noise_1[-(1:(nrow(snubh_comp_M)+1))]

for (i in 2:nrow(snubh_comp_Mp_1)) {
  snubh_comp_Mp_1[i,1] = snubh_comp_Mp_1[i-1,1] - snubh_comp_Mp_1[i-1,2] - snubh_comp_Mp_1[i-1,3]
}

## Excluding negative at-risk
snubh_comp_Mp_1 = snubh_comp_Mp_1[which(snubh_comp_Mp_1[,1] > 0),]  
snubh_comp_tk_1 = sum(snubh_comp_Mp_1[,1] > 0)

## Survival function
snubh_comp_dp_1_S = rep(NA, snubh_comp_tk_1) 

for (i in snubh_comp_time[1:snubh_comp_tk_1]) {
  j = which(snubh_comp_time == i)
  if (j == 1) {
    snubh_comp_dp_1_S[j] = (snubh_comp_Mp_1[1, 1] - snubh_comp_Mp_1[1, 2]) / snubh_comp_Mp_1[1, 1]
  } else {
    snubh_comp_dp_1_S[j] = snubh_comp_dp_1_S[j-1] * (snubh_comp_Mp_1[j, 1] - snubh_comp_Mp_1[j, 2]) / snubh_comp_Mp_1[j, 1]
  }
}

snubh_comp_dp_1_x = seq(1, snubh_comp_time[1])
snubh_comp_dp_1_y = rep(1, snubh_comp_time[1])

for (i in 1:(snubh_comp_tk_1-1)) {
  snubh_comp_dp_1_x = c(snubh_comp_dp_1_x, seq(snubh_comp_time[i], snubh_comp_time[i+1]))
  snubh_comp_dp_1_y = c(snubh_comp_dp_1_y, rep(snubh_comp_dp_1_S[i], snubh_comp_time[i+1]-snubh_comp_time[i]+1))
}

snubh_comp_dp_1_x = c(snubh_comp_dp_1_x, snubh_comp_time[snubh_comp_tk_1])
snubh_comp_dp_1_y = c(snubh_comp_dp_1_y, snubh_comp_dp_1_S[snubh_comp_tk_1])


### DP-tar(1)

snubh_tar_Mp_1 = snubh_tar_M

snubh_tar_noise_1 = rexp(nrow(snubh_tar_M)*2+1, eps[1]/S) * sign(runif(nrow(snubh_tar_M)*2+1, -1, 1))
snubh_tar_noise_1 = pmax(snubh_tar_noise_1, rep(0, nrow(snubh_tar_M)*2+1))

snubh_tar_Mp_1[1,1] = snubh_tar_Mp_1[1,1] + snubh_tar_noise_1[1]
snubh_tar_Mp_1[,2] = snubh_tar_M[,2] + snubh_tar_noise_1[2:(nrow(snubh_tar_M)+1)]
## snubh_tar_Mp_1[,3] = snubh_tar_M[,3] + snubh_tar_noise_1[-(1:(nrow(snubh_tar_M)+1))]

for (i in 2:nrow(snubh_tar_Mp_1)) {
  snubh_tar_Mp_1[i,1] = snubh_tar_Mp_1[i-1,1] - snubh_tar_Mp_1[i-1,2] - snubh_tar_Mp_1[i-1,3]
}

## Excluding negative at-risk
snubh_tar_Mp_1 = snubh_tar_Mp_1[which(snubh_tar_Mp_1[,1] > 0),]  
snubh_tar_tk_1 = sum(snubh_tar_Mp_1[,1] > 0)

## Survival function
snubh_tar_dp_1_S = rep(NA, snubh_tar_tk_1) 

for (i in snubh_tar_time[1:snubh_tar_tk_1]) {
  j = which(snubh_tar_time == i)
  if (j == 1) {
    snubh_tar_dp_1_S[j] = (snubh_tar_Mp_1[1, 1] - snubh_tar_Mp_1[1, 2]) / snubh_tar_Mp_1[1, 1]
  } else {
    snubh_tar_dp_1_S[j] = snubh_tar_dp_1_S[j-1] * (snubh_tar_Mp_1[j, 1] - snubh_tar_Mp_1[j, 2]) / snubh_tar_Mp_1[j, 1]
  }
}

snubh_tar_dp_1_x = seq(1, snubh_tar_time[1])
snubh_tar_dp_1_y = rep(1, snubh_tar_time[1])

for (i in 1:(snubh_tar_tk_1-1)) {
  snubh_tar_dp_1_x = c(snubh_tar_dp_1_x, seq(snubh_tar_time[i], snubh_tar_time[i+1]))
  snubh_tar_dp_1_y = c(snubh_tar_dp_1_y, rep(snubh_tar_dp_1_S[i], snubh_tar_time[i+1]-snubh_tar_time[i]+1))
}

snubh_tar_dp_1_x = c(snubh_tar_dp_1_x, snubh_tar_time[snubh_tar_tk_1])
snubh_tar_dp_1_y = c(snubh_tar_dp_1_y, snubh_tar_dp_1_S[snubh_tar_tk_1])


plot(snubh_comp_dp_1_x, snubh_comp_dp_1_y, type = 'l', xlim = c(1, snubh_timemax), 
     ylim = c(min(snubh_comp_dp_1_y, snubh_tar_dp_1_y) - 0.01, 1), lwd = 3, main = "SNUBH(DP, eps=1)")
points(snubh_tar_dp_1_x, snubh_tar_dp_1_y, type = 'l', lwd = 3, col ='red')
legend("topright", legend = c("Comparator", "Target"), col = c('black', 'red'), lwd = 3)


### eps = 2

### DP-comp(2)

## Constructing M'
snubh_comp_Mp_2 = snubh_comp_M

snubh_comp_noise_2 = rexp(nrow(snubh_comp_M)*2+1, eps[1]/S) * sign(runif(nrow(snubh_comp_M)*2+1, -1, 1))
snubh_comp_noise_2 = pmax(snubh_comp_noise_2, rep(0, nrow(snubh_comp_M)*2+1))

snubh_comp_Mp_2[1,1] = snubh_comp_Mp_2[1,1] + snubh_comp_noise_2[1]
snubh_comp_Mp_2[,2] = snubh_comp_M[,2] + snubh_comp_noise_2[2:(nrow(snubh_comp_M)+1)]
# snubh_comp_Mp_2[,3] = snubh_comp_M[,3] + snubh_comp_noise_2[-(1:(nrow(snubh_comp_M)+1))]

for (i in 2:nrow(snubh_comp_Mp_2)) {
  snubh_comp_Mp_2[i,1] = snubh_comp_Mp_2[i-1,1] - snubh_comp_Mp_2[i-1,2] - snubh_comp_Mp_2[i-1,3]
}

## Excluding negative at-risk
snubh_comp_Mp_2 = snubh_comp_Mp_2[which(snubh_comp_Mp_2[,1] > 0),]  
snubh_comp_tk_2 = sum(snubh_comp_Mp_2[,1] > 0)

## Survival function
snubh_comp_dp_2_S = rep(NA, snubh_comp_tk_2)

for (i in snubh_comp_time[1:snubh_comp_tk_2]) {
  j = which(snubh_comp_time == i)
  if (j == 1) {
    snubh_comp_dp_2_S[j] = (snubh_comp_Mp_2[1, 1] - snubh_comp_Mp_2[1, 2]) / snubh_comp_Mp_2[1, 1]
  } else {
    snubh_comp_dp_2_S[j] = snubh_comp_dp_2_S[j-1] * (snubh_comp_Mp_2[j, 1] - snubh_comp_Mp_2[j, 2]) / snubh_comp_Mp_2[j, 1]
  }
}

snubh_comp_dp_2_x = seq(1, snubh_comp_time[1])
snubh_comp_dp_2_y = rep(1, snubh_comp_time[1])

for (i in 1:(snubh_comp_tk_2-1)) {
  snubh_comp_dp_2_x = c(snubh_comp_dp_2_x, seq(snubh_comp_time[i], snubh_comp_time[i+1]))
  snubh_comp_dp_2_y = c(snubh_comp_dp_2_y, rep(snubh_comp_dp_2_S[i], snubh_comp_time[i+1]-snubh_comp_time[i]+1))
}

snubh_comp_dp_2_x = c(snubh_comp_dp_2_x, snubh_comp_time[snubh_comp_tk_2])
snubh_comp_dp_2_y = c(snubh_comp_dp_2_y, snubh_comp_dp_2_S[snubh_comp_tk_2])

### DP-tar(2)

snubh_tar_Mp_2 = snubh_tar_M

snubh_tar_noise_2 = rexp(nrow(snubh_tar_M)*2+1, eps[1]/S) * sign(runif(nrow(snubh_tar_M)*2+1, -1, 1))
snubh_tar_noise_2 = pmax(snubh_tar_noise_2, rep(0, nrow(snubh_tar_M)*2+1))

snubh_tar_Mp_2[1,1] = snubh_tar_Mp_2[1,1] + snubh_tar_noise_2[1]
snubh_tar_Mp_2[,2] = snubh_tar_M[,2] + snubh_tar_noise_2[2:(nrow(snubh_tar_M)+1)]
## snubh_tar_Mp_2[,3] = snubh_tar_M[,3] + snubh_tar_noise_2[-(1:(nrow(snubh_tar_M)+1))]

for (i in 2:nrow(snubh_tar_Mp_2)) {
  snubh_tar_Mp_2[i,1] = snubh_tar_Mp_2[i-1,1] - snubh_tar_Mp_2[i-1,2] - snubh_tar_Mp_2[i-1,3]
}

## Excluding negative at-risk
snubh_tar_Mp_2 = snubh_tar_Mp_2[which(snubh_tar_Mp_2[,1] > 0),]  
snubh_tar_tk_2 = sum(snubh_tar_Mp_2[,1] > 0)

## Survival function
snubh_tar_dp_2_S = rep(NA, snubh_tar_tk_2)

for (i in snubh_tar_time[1:snubh_tar_tk_2]) {
  j = which(snubh_tar_time == i)
  if (j == 1) {
    snubh_tar_dp_2_S[j] = (snubh_tar_Mp_2[1, 1] - snubh_tar_Mp_2[1, 2]) / snubh_tar_Mp_2[1, 1]
  } else {
    snubh_tar_dp_2_S[j] = snubh_tar_dp_2_S[j-1] * (snubh_tar_Mp_2[j, 1] - snubh_tar_Mp_2[j, 2]) / snubh_tar_Mp_2[j, 1]
  }
}

snubh_tar_dp_2_x = seq(1, snubh_tar_time[1])
snubh_tar_dp_2_y = rep(1, snubh_tar_time[1])

for (i in 1:(snubh_tar_tk_2-1)) {
  snubh_tar_dp_2_x = c(snubh_tar_dp_2_x, seq(snubh_tar_time[i], snubh_tar_time[i+1]))
  snubh_tar_dp_2_y = c(snubh_tar_dp_2_y, rep(snubh_tar_dp_2_S[i], snubh_tar_time[i+1]-snubh_tar_time[i]+1))
}

snubh_tar_dp_2_x = c(snubh_tar_dp_2_x, snubh_tar_time[snubh_tar_tk_2])
snubh_tar_dp_2_y = c(snubh_tar_dp_2_y, snubh_tar_dp_2_S[snubh_tar_tk_2])


plot(snubh_comp_dp_2_x, snubh_comp_dp_2_y, type = 'l', xlim = c(1, snubh_timemax), 
     ylim = c(min(snubh_comp_dp_2_y, snubh_tar_dp_2_y) - 0.01, 1), lwd = 3, main = "SNUBH(DP, eps=2)")
points(snubh_tar_dp_2_x, snubh_tar_dp_2_y, type = 'l', lwd = 3, col ='red')
legend("topright", legend = c("Comparator", "Target"), col = c('black', 'red'), lwd = 3)


### eps = 3

### DP-comp(3)

snubh_comp_Mp_3 = snubh_comp_M

snubh_comp_noise_3 = rexp(nrow(snubh_comp_M)*2+1, eps[1]/S) * sign(runif(nrow(snubh_comp_M)*2+1, -1, 1))
snubh_comp_noise_3 = pmax(snubh_comp_noise_3, rep(0, nrow(snubh_comp_M)*2+1))

snubh_comp_Mp_3[1,1] = snubh_comp_Mp_3[1,1] + snubh_comp_noise_3[1]
snubh_comp_Mp_3[,2] = snubh_comp_M[,2] + snubh_comp_noise_3[2:(nrow(snubh_comp_M)+1)]
## snubh_comp_Mp_3[,3] = snubh_comp_M[,3] + snubh_comp_noise_3[-(1:(nrow(snubh_comp_M)+1))]

for (i in 2:nrow(snubh_comp_Mp_3)) {
  snubh_comp_Mp_3[i,1] = snubh_comp_Mp_3[i-1,1] - snubh_comp_Mp_3[i-1,2] - snubh_comp_Mp_3[i-1,3]
}

## Excluding negative at-risk
snubh_comp_Mp_3 = snubh_comp_Mp_3[which(snubh_comp_Mp_3[,1] > 0),]  
snubh_comp_tk_3 = sum(snubh_comp_Mp_3[,1] > 0)

## Survival function
snubh_comp_dp_3_S = rep(NA, snubh_comp_tk_3) 

for (i in snubh_comp_time[1:snubh_comp_tk_3]) {
  j = which(snubh_comp_time == i)
  if (j == 1) {
    snubh_comp_dp_3_S[j] = (snubh_comp_Mp_3[1, 1] - snubh_comp_Mp_3[1, 2]) / snubh_comp_Mp_3[1, 1]
  } else {
    snubh_comp_dp_3_S[j] = snubh_comp_dp_3_S[j-1] * (snubh_comp_Mp_3[j, 1] - snubh_comp_Mp_3[j, 2]) / snubh_comp_Mp_3[j, 1]
  }
}

snubh_comp_dp_3_x = seq(1, snubh_comp_time[1])
snubh_comp_dp_3_y = rep(1, snubh_comp_time[1])

for (i in 1:(snubh_comp_tk_3-1)) {
  snubh_comp_dp_3_x = c(snubh_comp_dp_3_x, seq(snubh_comp_time[i], snubh_comp_time[i+1]))
  snubh_comp_dp_3_y = c(snubh_comp_dp_3_y, rep(snubh_comp_dp_3_S[i], snubh_comp_time[i+1]-snubh_comp_time[i]+1))
}

snubh_comp_dp_3_x = c(snubh_comp_dp_3_x, snubh_comp_time[snubh_comp_tk_3])
snubh_comp_dp_3_y = c(snubh_comp_dp_3_y, snubh_comp_dp_3_S[snubh_comp_tk_3])


### DP-tar(3)

snubh_tar_Mp_3 = snubh_tar_M

snubh_tar_noise_3 = rexp(nrow(snubh_tar_M)*2+1, eps[1]/S) * sign(runif(nrow(snubh_tar_M)*2+1, -1, 1))
snubh_tar_noise_3 = pmax(snubh_tar_noise_3, rep(0, nrow(snubh_tar_M)*2+1))

snubh_tar_Mp_3[1,1] = snubh_tar_Mp_3[1,1] + snubh_tar_noise_3[1]
snubh_tar_Mp_3[,2] = snubh_tar_M[,2] + snubh_tar_noise_3[2:(nrow(snubh_tar_M)+1)]
## snubh_tar_Mp_3[,3] = snubh_tar_M[,3] + snubh_tar_noise_3[-(1:(nrow(snubh_tar_M)+1))]

for (i in 2:nrow(snubh_tar_Mp_3)) {
  snubh_tar_Mp_3[i,1] = snubh_tar_Mp_3[i-1,1] - snubh_tar_Mp_3[i-1,2] - snubh_tar_Mp_3[i-1,3]
}

## Excluding negative at-risk
snubh_tar_Mp_3 = snubh_tar_Mp_3[which(snubh_tar_Mp_3[,1] > 0),]  
snubh_tar_tk_3 = sum(snubh_tar_Mp_3[,1] > 0)

## Survival function
snubh_tar_dp_3_S = rep(NA, snubh_tar_tk_3) 

for (i in snubh_tar_time[1:snubh_tar_tk_3]) {
  j = which(snubh_tar_time == i)
  if (j == 1) {
    snubh_tar_dp_3_S[j] = (snubh_tar_Mp_3[1, 1] - snubh_tar_Mp_3[1, 2]) / snubh_tar_Mp_3[1, 1]
  } else {
    snubh_tar_dp_3_S[j] = snubh_tar_dp_3_S[j-1] * (snubh_tar_Mp_3[j, 1] - snubh_tar_Mp_3[j, 2]) / snubh_tar_Mp_3[j, 1]
  }
}

snubh_tar_dp_3_x = seq(1, snubh_tar_time[1])
snubh_tar_dp_3_y = rep(1, snubh_tar_time[1])

for (i in 1:(snubh_tar_tk_3-1)) {
  snubh_tar_dp_3_x = c(snubh_tar_dp_3_x, seq(snubh_tar_time[i], snubh_tar_time[i+1]))
  snubh_tar_dp_3_y = c(snubh_tar_dp_3_y, rep(snubh_tar_dp_3_S[i], snubh_tar_time[i+1]-snubh_tar_time[i]+1))
}

snubh_tar_dp_3_x = c(snubh_tar_dp_3_x, snubh_tar_time[snubh_tar_tk_3])
snubh_tar_dp_3_y = c(snubh_tar_dp_3_y, snubh_tar_dp_3_S[snubh_tar_tk_3])


plot(snubh_comp_dp_3_x, snubh_comp_dp_3_y, type = 'l', xlim = c(1, snubh_timemax), 
     ylim = c(min(snubh_comp_dp_3_y, snubh_tar_dp_3_y) - 0.01, 1), lwd = 3, main = "SNUBH(DP, eps=3)")
points(snubh_tar_dp_3_x, snubh_tar_dp_3_y, type = 'l', lwd = 3, col ='red')
legend("topright", legend = c("Comparator", "Target"), col = c('black', 'red'), lwd = 3)


### Compare by eps

# Comparator
plot(snubh_comp_x, snubh_comp_y, type = 'l', xlim = c(1, snubh_timemax), 
     ylim = c(min(snubh_comp_y, snubh_comp_dp_1_y, snubh_comp_dp_2_y, snubh_comp_dp_3_y) - 0.01, 1),
     lwd = 3, main = "SNUBH(Comparator, by eps)")
points(snubh_comp_dp_1_x, snubh_comp_dp_1_y, type = 'l', lwd = 3, col ='blue')
points(snubh_comp_dp_2_x, snubh_comp_dp_2_y, type = 'l', lwd = 3, col ='purple')
points(snubh_comp_dp_3_x, snubh_comp_dp_3_y, type = 'l', lwd = 3, col ='green')
legend("topright", legend = c("NA", "1", "2", "3"), col = c('black', 'blue', 'purple', 'green'), lwd = 3)

# Target
plot(snubh_tar_x, snubh_tar_y, type = 'l', xlim = c(1, snubh_timemax), 
     ylim = c(min(snubh_tar_y, snubh_tar_dp_1_y, snubh_tar_dp_2_y, snubh_tar_dp_3_y) - 0.01, 1),
     lwd = 3, main = "SNUBH(Target, by eps)")
points(snubh_tar_dp_1_x, snubh_tar_dp_1_y, type = 'l', lwd = 3, col ='blue')
points(snubh_tar_dp_2_x, snubh_tar_dp_2_y, type = 'l', lwd = 3, col ='purple')
points(snubh_tar_dp_3_x, snubh_tar_dp_3_y, type = 'l', lwd = 3, col ='green')
legend("topright", legend = c("NA", "1", "2", "3"), col = c('black', 'blue', 'purple', 'green'), lwd = 3)







##### SNUH #####

snuh_comp = mtsnuh[mtsnuh$treatment == 0, ]
snuh_tar = mtsnuh[mtsnuh$treatment == 1, ]

snuh_timemax = max(mtsnuh$survivalTime)

### Comparator

snuh_comp_time = unique(sort(snuh_comp$survivalTime[which(snuh_comp$outcomeCount != 0)]))
snuh_comp_tk = length(snuh_comp_time)

snuh_comp_M = matrix(rep(NA, snuh_comp_tk*3), nrow = snuh_comp_tk, ncol = 3)

for (i in 1:snuh_comp_tk) {
  snuh_comp_M[i,1] = sum(snuh_comp$survivalTime >= snuh_comp_time[i])
  snuh_comp_M[i,2] = length(which(snuh_comp$survivalTime[which(snuh_comp$outcomeCount != 0)] == snuh_comp_time[i]))
  if(i >= 2) {
    snuh_comp_M[i-1,3] = snuh_comp_M[i-1,1] - snuh_comp_M[i,1] - snuh_comp_M[i-1,2]
  }
}

## Survival function
snuh_comp_S = rep(NA, snuh_comp_tk) 

for (i in snuh_comp_time) {
  j = which(snuh_comp_time == i)
  if (j == 1) {
    snuh_comp_S[j] = (snuh_comp_M[1, 1] - snuh_comp_M[1, 2]) / snuh_comp_M[1, 1]
  } else {
    snuh_comp_S[j] = snuh_comp_S[j-1] * (snuh_comp_M[j, 1] - snuh_comp_M[j, 2]) / snuh_comp_M[j, 1]
  }
}

snuh_comp_x = seq(1, snuh_comp_time[1])
snuh_comp_y = rep(1, snuh_comp_time[1])

for (i in 1:(snuh_comp_tk-1)) {
  snuh_comp_x = c(snuh_comp_x, seq(snuh_comp_time[i], snuh_comp_time[i+1]))
  snuh_comp_y = c(snuh_comp_y, rep(snuh_comp_S[i], snuh_comp_time[i+1]-snuh_comp_time[i]+1))
}

snuh_comp_x = c(snuh_comp_x, snuh_comp_time[snuh_comp_tk])
snuh_comp_y = c(snuh_comp_y, snuh_comp_S[snuh_comp_tk])


### Target

snuh_tar_time = unique(sort(snuh_tar$survivalTime[which(snuh_tar$outcomeCount != 0)]))
snuh_tar_tk = length(snuh_tar_time)

snuh_tar_M = matrix(rep(NA, snuh_tar_tk*3), nrow = snuh_tar_tk, ncol = 3)

for (i in 1:snuh_tar_tk) {
  snuh_tar_M[i,1] = sum(snuh_tar$survivalTime >= snuh_tar_time[i])
  snuh_tar_M[i,2] = length(which(snuh_tar$survivalTime[which(snuh_tar$outcomeCount != 0)] == snuh_tar_time[i]))
  if(i >= 2) {
    snuh_tar_M[i-1,3] = snuh_tar_M[i-1,1] - snuh_tar_M[i,1] - snuh_tar_M[i-1,2]
  }
}

snuh_tar_S = rep(NA, snuh_tar_tk)

for (i in snuh_tar_time) {
  j = which(snuh_tar_time == i)
  if (j == 1) {
    snuh_tar_S[j] = (snuh_tar_M[1, 1] - snuh_tar_M[1, 2]) / snuh_tar_M[1, 1]
  } else {
    snuh_tar_S[j] = snuh_tar_S[j-1] * (snuh_tar_M[j, 1] - snuh_tar_M[j, 2]) / snuh_tar_M[j, 1]
  }
}

snuh_tar_x = seq(1, snuh_tar_time[1])
snuh_tar_y = rep(1, snuh_tar_time[1])

for (i in 1:(snuh_tar_tk-1)) {
  snuh_tar_x = c(snuh_tar_x, seq(snuh_tar_time[i], snuh_tar_time[i+1]))
  snuh_tar_y = c(snuh_tar_y, rep(snuh_tar_S[i], snuh_tar_time[i+1]-snuh_tar_time[i]+1))
}

snuh_tar_x = c(snuh_tar_x, snuh_tar_time[snuh_tar_tk])
snuh_tar_y = c(snuh_tar_y, snuh_tar_S[snuh_tar_tk])

### Plot

plot(snuh_comp_x, snuh_comp_y, type = 'l', xlim = c(1, snuh_timemax), 
     ylim = c(min(snuh_comp_y, snuh_tar_y) - 0.01, 1), lwd = 3, main = "SNUH(Original)")
points(snuh_tar_x, snuh_tar_y, type = 'l', lwd = 3, col ='red')
legend("topright", legend=c("Comparator", "Target"), col = c('black', 'red'), lwd = 3)

################ DP

set.seed(2)
S = 2 ## L1 Sensitivity
eps = c(1, 2, 3)

### eps = 1

### DP-comp(1)

snuh_comp_Mp_1 = snuh_comp_M

snuh_comp_noise_1 = rexp(nrow(snuh_comp_M)*2+1, eps[1]/S) * sign(runif(nrow(snuh_comp_M)*2+1, -1, 1))
snuh_comp_noise_1 = pmax(snuh_comp_noise_1, rep(0, nrow(snuh_comp_M)*2+1))

snuh_comp_Mp_1[1,1] = snuh_comp_Mp_1[1,1] + snuh_comp_noise_1[1]
snuh_comp_Mp_1[,2] = snuh_comp_M[,2] + snuh_comp_noise_1[2:(nrow(snuh_comp_M)+1)]
## snuh_comp_Mp_1[,3] = snuh_comp_M[,3] + snuh_comp_noise_1[-(1:(nrow(snuh_comp_M)+1))]

for (i in 2:nrow(snuh_comp_Mp_1)) {
  snuh_comp_Mp_1[i,1] = snuh_comp_Mp_1[i-1,1] - snuh_comp_Mp_1[i-1,2] - snuh_comp_Mp_1[i-1,3]
}

## Excluding negative at-risk
snuh_comp_Mp_1 = snuh_comp_Mp_1[which(snuh_comp_Mp_1[,1] > 0),]  
snuh_comp_tk_1 = sum(snuh_comp_Mp_1[,1] > 0)

## Survival function
snuh_comp_dp_1_S = rep(NA, snuh_comp_tk_1) 

for (i in snuh_comp_time[1:snuh_comp_tk_1]) {
  j = which(snuh_comp_time == i)
  if (j == 1) {
    snuh_comp_dp_1_S[j] = (snuh_comp_Mp_1[1, 1] - snuh_comp_Mp_1[1, 2]) / snuh_comp_Mp_1[1, 1]
  } else {
    snuh_comp_dp_1_S[j] = snuh_comp_dp_1_S[j-1] * (snuh_comp_Mp_1[j, 1] - snuh_comp_Mp_1[j, 2]) / snuh_comp_Mp_1[j, 1]
  }
}

snuh_comp_dp_1_x = seq(1, snuh_comp_time[1])
snuh_comp_dp_1_y = rep(1, snuh_comp_time[1])

for (i in 1:(snuh_comp_tk_1-1)) {
  snuh_comp_dp_1_x = c(snuh_comp_dp_1_x, seq(snuh_comp_time[i], snuh_comp_time[i+1]))
  snuh_comp_dp_1_y = c(snuh_comp_dp_1_y, rep(snuh_comp_dp_1_S[i], snuh_comp_time[i+1]-snuh_comp_time[i]+1))
}

snuh_comp_dp_1_x = c(snuh_comp_dp_1_x, snuh_comp_time[snuh_comp_tk_1])
snuh_comp_dp_1_y = c(snuh_comp_dp_1_y, snuh_comp_dp_1_S[snuh_comp_tk_1])


### DP-tar(1)

snuh_tar_Mp_1 = snuh_tar_M

snuh_tar_noise_1 = rexp(nrow(snuh_tar_M)*2+1, eps[1]/S) * sign(runif(nrow(snuh_tar_M)*2+1, -1, 1))
snuh_tar_noise_1 = pmax(snuh_tar_noise_1, rep(0, nrow(snuh_tar_M)*2+1))

snuh_tar_Mp_1[1,1] = snuh_tar_Mp_1[1,1] + snuh_tar_noise_1[1]
snuh_tar_Mp_1[,2] = snuh_tar_M[,2] + snuh_tar_noise_1[2:(nrow(snuh_tar_M)+1)]
## snuh_tar_Mp_1[,3] = snuh_tar_M[,3] + snuh_tar_noise_1[-(1:(nrow(snuh_tar_M)+1))]

for (i in 2:nrow(snuh_tar_Mp_1)) {
  snuh_tar_Mp_1[i,1] = snuh_tar_Mp_1[i-1,1] - snuh_tar_Mp_1[i-1,2] - snuh_tar_Mp_1[i-1,3]
}

## Excluding negative at-risk
snuh_tar_Mp_1 = snuh_tar_Mp_1[which(snuh_tar_Mp_1[,1] > 0),]  
snuh_tar_tk_1 = sum(snuh_tar_Mp_1[,1] > 0)

## Survival function
snuh_tar_dp_1_S = rep(NA, snuh_tar_tk_1) 

for (i in snuh_tar_time[1:snuh_tar_tk_1]) {
  j = which(snuh_tar_time == i)
  if (j == 1) {
    snuh_tar_dp_1_S[j] = (snuh_tar_Mp_1[1, 1] - snuh_tar_Mp_1[1, 2]) / snuh_tar_Mp_1[1, 1]
  } else {
    snuh_tar_dp_1_S[j] = snuh_tar_dp_1_S[j-1] * (snuh_tar_Mp_1[j, 1] - snuh_tar_Mp_1[j, 2]) / snuh_tar_Mp_1[j, 1]
  }
}

snuh_tar_dp_1_x = seq(1, snuh_tar_time[1])
snuh_tar_dp_1_y = rep(1, snuh_tar_time[1])

for (i in 1:(snuh_tar_tk_1-1)) {
  snuh_tar_dp_1_x = c(snuh_tar_dp_1_x, seq(snuh_tar_time[i], snuh_tar_time[i+1]))
  snuh_tar_dp_1_y = c(snuh_tar_dp_1_y, rep(snuh_tar_dp_1_S[i], snuh_tar_time[i+1]-snuh_tar_time[i]+1))
}

snuh_tar_dp_1_x = c(snuh_tar_dp_1_x, snuh_tar_time[snuh_tar_tk_1])
snuh_tar_dp_1_y = c(snuh_tar_dp_1_y, snuh_tar_dp_1_S[snuh_tar_tk_1])


plot(snuh_comp_dp_1_x, snuh_comp_dp_1_y, type = 'l', xlim = c(1, snuh_timemax), 
     ylim = c(min(snuh_comp_dp_1_y, snuh_tar_dp_1_y) - 0.01, 1), lwd = 3, main = "SNUH(DP, eps=1)")
points(snuh_tar_dp_1_x, snuh_tar_dp_1_y, type = 'l', lwd = 3, col ='red')
legend("topright", legend = c("Comparator", "Target"), col = c('black', 'red'), lwd = 3)


### eps = 2

### DP-comp(2)

## Constructing M'
snuh_comp_Mp_2 = snuh_comp_M

snuh_comp_noise_2 = rexp(nrow(snuh_comp_M)*2+1, eps[1]/S) * sign(runif(nrow(snuh_comp_M)*2+1, -1, 1))
snuh_comp_noise_2 = pmax(snuh_comp_noise_2, rep(0, nrow(snuh_comp_M)*2+1))

snuh_comp_Mp_2[1,1] = snuh_comp_Mp_2[1,1] + snuh_comp_noise_2[1]
snuh_comp_Mp_2[,2] = snuh_comp_M[,2] + snuh_comp_noise_2[2:(nrow(snuh_comp_M)+1)]
# snuh_comp_Mp_2[,3] = snuh_comp_M[,3] + snuh_comp_noise_2[-(1:(nrow(snuh_comp_M)+1))]

for (i in 2:nrow(snuh_comp_Mp_2)) {
  snuh_comp_Mp_2[i,1] = snuh_comp_Mp_2[i-1,1] - snuh_comp_Mp_2[i-1,2] - snuh_comp_Mp_2[i-1,3]
}

## Excluding negative at-risk
snuh_comp_Mp_2 = snuh_comp_Mp_2[which(snuh_comp_Mp_2[,1] > 0),]  
snuh_comp_tk_2 = sum(snuh_comp_Mp_2[,1] > 0)

## Survival function
snuh_comp_dp_2_S = rep(NA, snuh_comp_tk_2)

for (i in snuh_comp_time[1:snuh_comp_tk_2]) {
  j = which(snuh_comp_time == i)
  if (j == 1) {
    snuh_comp_dp_2_S[j] = (snuh_comp_Mp_2[1, 1] - snuh_comp_Mp_2[1, 2]) / snuh_comp_Mp_2[1, 1]
  } else {
    snuh_comp_dp_2_S[j] = snuh_comp_dp_2_S[j-1] * (snuh_comp_Mp_2[j, 1] - snuh_comp_Mp_2[j, 2]) / snuh_comp_Mp_2[j, 1]
  }
}

snuh_comp_dp_2_x = seq(1, snuh_comp_time[1])
snuh_comp_dp_2_y = rep(1, snuh_comp_time[1])

for (i in 1:(snuh_comp_tk_2-1)) {
  snuh_comp_dp_2_x = c(snuh_comp_dp_2_x, seq(snuh_comp_time[i], snuh_comp_time[i+1]))
  snuh_comp_dp_2_y = c(snuh_comp_dp_2_y, rep(snuh_comp_dp_2_S[i], snuh_comp_time[i+1]-snuh_comp_time[i]+1))
}

snuh_comp_dp_2_x = c(snuh_comp_dp_2_x, snuh_comp_time[snuh_comp_tk_2])
snuh_comp_dp_2_y = c(snuh_comp_dp_2_y, snuh_comp_dp_2_S[snuh_comp_tk_2])

### DP-tar(2)

snuh_tar_Mp_2 = snuh_tar_M

snuh_tar_noise_2 = rexp(nrow(snuh_tar_M)*2+1, eps[1]/S) * sign(runif(nrow(snuh_tar_M)*2+1, -1, 1))
snuh_tar_noise_2 = pmax(snuh_tar_noise_2, rep(0, nrow(snuh_tar_M)*2+1))

snuh_tar_Mp_2[1,1] = snuh_tar_Mp_2[1,1] + snuh_tar_noise_2[1]
snuh_tar_Mp_2[,2] = snuh_tar_M[,2] + snuh_tar_noise_2[2:(nrow(snuh_tar_M)+1)]
## snuh_tar_Mp_2[,3] = snuh_tar_M[,3] + snuh_tar_noise_2[-(1:(nrow(snuh_tar_M)+1))]

for (i in 2:nrow(snuh_tar_Mp_2)) {
  snuh_tar_Mp_2[i,1] = snuh_tar_Mp_2[i-1,1] - snuh_tar_Mp_2[i-1,2] - snuh_tar_Mp_2[i-1,3]
}

## Excluding negative at-risk
snuh_tar_Mp_2 = snuh_tar_Mp_2[which(snuh_tar_Mp_2[,1] > 0),]  
snuh_tar_tk_2 = sum(snuh_tar_Mp_2[,1] > 0)

## Survival function
snuh_tar_dp_2_S = rep(NA, snuh_tar_tk_2)

for (i in snuh_tar_time[1:snuh_tar_tk_2]) {
  j = which(snuh_tar_time == i)
  if (j == 1) {
    snuh_tar_dp_2_S[j] = (snuh_tar_Mp_2[1, 1] - snuh_tar_Mp_2[1, 2]) / snuh_tar_Mp_2[1, 1]
  } else {
    snuh_tar_dp_2_S[j] = snuh_tar_dp_2_S[j-1] * (snuh_tar_Mp_2[j, 1] - snuh_tar_Mp_2[j, 2]) / snuh_tar_Mp_2[j, 1]
  }
}

snuh_tar_dp_2_x = seq(1, snuh_tar_time[1])
snuh_tar_dp_2_y = rep(1, snuh_tar_time[1])

for (i in 1:(snuh_tar_tk_2-1)) {
  snuh_tar_dp_2_x = c(snuh_tar_dp_2_x, seq(snuh_tar_time[i], snuh_tar_time[i+1]))
  snuh_tar_dp_2_y = c(snuh_tar_dp_2_y, rep(snuh_tar_dp_2_S[i], snuh_tar_time[i+1]-snuh_tar_time[i]+1))
}

snuh_tar_dp_2_x = c(snuh_tar_dp_2_x, snuh_tar_time[snuh_tar_tk_2])
snuh_tar_dp_2_y = c(snuh_tar_dp_2_y, snuh_tar_dp_2_S[snuh_tar_tk_2])


plot(snuh_comp_dp_2_x, snuh_comp_dp_2_y, type = 'l', xlim = c(1, snuh_timemax), 
     ylim = c(min(snuh_comp_dp_2_y, snuh_tar_dp_2_y) - 0.01, 1), lwd = 3, main = "SNUH(DP, eps=2)")
points(snuh_tar_dp_2_x, snuh_tar_dp_2_y, type = 'l', lwd = 3, col ='red')
legend("topright", legend = c("Comparator", "Target"), col = c('black', 'red'), lwd = 3)


### eps = 3

### DP-comp(3)

snuh_comp_Mp_3 = snuh_comp_M

snuh_comp_noise_3 = rexp(nrow(snuh_comp_M)*2+1, eps[1]/S) * sign(runif(nrow(snuh_comp_M)*2+1, -1, 1))
snuh_comp_noise_3 = pmax(snuh_comp_noise_3, rep(0, nrow(snuh_comp_M)*2+1))

snuh_comp_Mp_3[1,1] = snuh_comp_Mp_3[1,1] + snuh_comp_noise_3[1]
snuh_comp_Mp_3[,2] = snuh_comp_M[,2] + snuh_comp_noise_3[2:(nrow(snuh_comp_M)+1)]
## snuh_comp_Mp_3[,3] = snuh_comp_M[,3] + snuh_comp_noise_3[-(1:(nrow(snuh_comp_M)+1))]

for (i in 2:nrow(snuh_comp_Mp_3)) {
  snuh_comp_Mp_3[i,1] = snuh_comp_Mp_3[i-1,1] - snuh_comp_Mp_3[i-1,2] - snuh_comp_Mp_3[i-1,3]
}

## Excluding negative at-risk
snuh_comp_Mp_3 = snuh_comp_Mp_3[which(snuh_comp_Mp_3[,1] > 0),]  
snuh_comp_tk_3 = sum(snuh_comp_Mp_3[,1] > 0)

## Survival function
snuh_comp_dp_3_S = rep(NA, snuh_comp_tk_3) 

for (i in snuh_comp_time[1:snuh_comp_tk_3]) {
  j = which(snuh_comp_time == i)
  if (j == 1) {
    snuh_comp_dp_3_S[j] = (snuh_comp_Mp_3[1, 1] - snuh_comp_Mp_3[1, 2]) / snuh_comp_Mp_3[1, 1]
  } else {
    snuh_comp_dp_3_S[j] = snuh_comp_dp_3_S[j-1] * (snuh_comp_Mp_3[j, 1] - snuh_comp_Mp_3[j, 2]) / snuh_comp_Mp_3[j, 1]
  }
}

snuh_comp_dp_3_x = seq(1, snuh_comp_time[1])
snuh_comp_dp_3_y = rep(1, snuh_comp_time[1])

for (i in 1:(snuh_comp_tk_3-1)) {
  snuh_comp_dp_3_x = c(snuh_comp_dp_3_x, seq(snuh_comp_time[i], snuh_comp_time[i+1]))
  snuh_comp_dp_3_y = c(snuh_comp_dp_3_y, rep(snuh_comp_dp_3_S[i], snuh_comp_time[i+1]-snuh_comp_time[i]+1))
}

snuh_comp_dp_3_x = c(snuh_comp_dp_3_x, snuh_comp_time[snuh_comp_tk_3])
snuh_comp_dp_3_y = c(snuh_comp_dp_3_y, snuh_comp_dp_3_S[snuh_comp_tk_3])


### DP-tar(3)

snuh_tar_Mp_3 = snuh_tar_M

snuh_tar_noise_3 = rexp(nrow(snuh_tar_M)*2+1, eps[1]/S) * sign(runif(nrow(snuh_tar_M)*2+1, -1, 1))
snuh_tar_noise_3 = pmax(snuh_tar_noise_3, rep(0, nrow(snuh_tar_M)*2+1))

snuh_tar_Mp_3[1,1] = snuh_tar_Mp_3[1,1] + snuh_tar_noise_3[1]
snuh_tar_Mp_3[,2] = snuh_tar_M[,2] + snuh_tar_noise_3[2:(nrow(snuh_tar_M)+1)]
## snuh_tar_Mp_3[,3] = snuh_tar_M[,3] + snuh_tar_noise_3[-(1:(nrow(snuh_tar_M)+1))]

for (i in 2:nrow(snuh_tar_Mp_3)) {
  snuh_tar_Mp_3[i,1] = snuh_tar_Mp_3[i-1,1] - snuh_tar_Mp_3[i-1,2] - snuh_tar_Mp_3[i-1,3]
}

## Excluding negative at-risk
snuh_tar_Mp_3 = snuh_tar_Mp_3[which(snuh_tar_Mp_3[,1] > 0),]  
snuh_tar_tk_3 = sum(snuh_tar_Mp_3[,1] > 0)

## Survival function
snuh_tar_dp_3_S = rep(NA, snuh_tar_tk_3) 

for (i in snuh_tar_time[1:snuh_tar_tk_3]) {
  j = which(snuh_tar_time == i)
  if (j == 1) {
    snuh_tar_dp_3_S[j] = (snuh_tar_Mp_3[1, 1] - snuh_tar_Mp_3[1, 2]) / snuh_tar_Mp_3[1, 1]
  } else {
    snuh_tar_dp_3_S[j] = snuh_tar_dp_3_S[j-1] * (snuh_tar_Mp_3[j, 1] - snuh_tar_Mp_3[j, 2]) / snuh_tar_Mp_3[j, 1]
  }
}

snuh_tar_dp_3_x = seq(1, snuh_tar_time[1])
snuh_tar_dp_3_y = rep(1, snuh_tar_time[1])

for (i in 1:(snuh_tar_tk_3-1)) {
  snuh_tar_dp_3_x = c(snuh_tar_dp_3_x, seq(snuh_tar_time[i], snuh_tar_time[i+1]))
  snuh_tar_dp_3_y = c(snuh_tar_dp_3_y, rep(snuh_tar_dp_3_S[i], snuh_tar_time[i+1]-snuh_tar_time[i]+1))
}

snuh_tar_dp_3_x = c(snuh_tar_dp_3_x, snuh_tar_time[snuh_tar_tk_3])
snuh_tar_dp_3_y = c(snuh_tar_dp_3_y, snuh_tar_dp_3_S[snuh_tar_tk_3])


plot(snuh_comp_dp_3_x, snuh_comp_dp_3_y, type = 'l', xlim = c(1, snuh_timemax), 
     ylim = c(min(snuh_comp_dp_3_y, snuh_tar_dp_3_y) - 0.01, 1), lwd = 3, main = "SNUH(DP, eps=3)")
points(snuh_tar_dp_3_x, snuh_tar_dp_3_y, type = 'l', lwd = 3, col ='red')
legend("topright", legend = c("Comparator", "Target"), col = c('black', 'red'), lwd = 3)


### Compare by eps

# Comparator
plot(snuh_comp_x, snuh_comp_y, type = 'l', xlim = c(1, snuh_timemax), 
     ylim = c(min(snuh_comp_y, snuh_comp_dp_1_y, snuh_comp_dp_2_y, snuh_comp_dp_3_y) - 0.01, 1),
     lwd = 3, main = "SNUH(Comparator, by eps)")
points(snuh_comp_dp_1_x, snuh_comp_dp_1_y, type = 'l', lwd = 3, col ='blue')
points(snuh_comp_dp_2_x, snuh_comp_dp_2_y, type = 'l', lwd = 3, col ='purple')
points(snuh_comp_dp_3_x, snuh_comp_dp_3_y, type = 'l', lwd = 3, col ='green')
legend("topright", legend = c("NA", "1", "2", "3"), col = c('black', 'blue', 'purple', 'green'), lwd = 3)

# Target
plot(snuh_tar_x, snuh_tar_y, type = 'l', xlim = c(1, snuh_timemax), 
     ylim = c(min(snuh_tar_y, snuh_tar_dp_1_y, snuh_tar_dp_2_y, snuh_tar_dp_3_y) - 0.01, 1),
     lwd = 3, main = "SNUH(Target, by eps)")
points(snuh_tar_dp_1_x, snuh_tar_dp_1_y, type = 'l', lwd = 3, col ='blue')
points(snuh_tar_dp_2_x, snuh_tar_dp_2_y, type = 'l', lwd = 3, col ='purple')
points(snuh_tar_dp_3_x, snuh_tar_dp_3_y, type = 'l', lwd = 3, col ='green')
legend("topright", legend = c("NA", "1", "2", "3"), col = c('black', 'blue', 'purple', 'green'), lwd = 3)









##### CMCS #####

cmcs_comp = mtcmcs[mtcmcs$treatment == 0, ]
cmcs_tar = mtcmcs[mtcmcs$treatment == 1, ]

cmcs_timemax = max(mtcmcs$survivalTime)

### Comparator

cmcs_comp_time = unique(sort(cmcs_comp$survivalTime[which(cmcs_comp$outcomeCount != 0)]))
cmcs_comp_tk = length(cmcs_comp_time)

cmcs_comp_M = matrix(rep(NA, cmcs_comp_tk*3), nrow = cmcs_comp_tk, ncol = 3)

for (i in 1:cmcs_comp_tk) {
  cmcs_comp_M[i,1] = sum(cmcs_comp$survivalTime >= cmcs_comp_time[i])
  cmcs_comp_M[i,2] = length(which(cmcs_comp$survivalTime[which(cmcs_comp$outcomeCount != 0)] == cmcs_comp_time[i]))
  if(i >= 2) {
    cmcs_comp_M[i-1,3] = cmcs_comp_M[i-1,1] - cmcs_comp_M[i,1] - cmcs_comp_M[i-1,2]
  }
}

## Survival function
cmcs_comp_S = rep(NA, cmcs_comp_tk) 

for (i in cmcs_comp_time) {
  j = which(cmcs_comp_time == i)
  if (j == 1) {
    cmcs_comp_S[j] = (cmcs_comp_M[1, 1] - cmcs_comp_M[1, 2]) / cmcs_comp_M[1, 1]
  } else {
    cmcs_comp_S[j] = cmcs_comp_S[j-1] * (cmcs_comp_M[j, 1] - cmcs_comp_M[j, 2]) / cmcs_comp_M[j, 1]
  }
}

cmcs_comp_x = seq(1, cmcs_comp_time[1])
cmcs_comp_y = rep(1, cmcs_comp_time[1])

for (i in 1:(cmcs_comp_tk-1)) {
  cmcs_comp_x = c(cmcs_comp_x, seq(cmcs_comp_time[i], cmcs_comp_time[i+1]))
  cmcs_comp_y = c(cmcs_comp_y, rep(cmcs_comp_S[i], cmcs_comp_time[i+1]-cmcs_comp_time[i]+1))
}

cmcs_comp_x = c(cmcs_comp_x, cmcs_comp_time[cmcs_comp_tk])
cmcs_comp_y = c(cmcs_comp_y, cmcs_comp_S[cmcs_comp_tk])


### Target

cmcs_tar_time = unique(sort(cmcs_tar$survivalTime[which(cmcs_tar$outcomeCount != 0)]))
cmcs_tar_tk = length(cmcs_tar_time)

cmcs_tar_M = matrix(rep(NA, cmcs_tar_tk*3), nrow = cmcs_tar_tk, ncol = 3)

for (i in 1:cmcs_tar_tk) {
  cmcs_tar_M[i,1] = sum(cmcs_tar$survivalTime >= cmcs_tar_time[i])
  cmcs_tar_M[i,2] = length(which(cmcs_tar$survivalTime[which(cmcs_tar$outcomeCount != 0)] == cmcs_tar_time[i]))
  if(i >= 2) {
    cmcs_tar_M[i-1,3] = cmcs_tar_M[i-1,1] - cmcs_tar_M[i,1] - cmcs_tar_M[i-1,2]
  }
}

cmcs_tar_S = rep(NA, cmcs_tar_tk)

for (i in cmcs_tar_time) {
  j = which(cmcs_tar_time == i)
  if (j == 1) {
    cmcs_tar_S[j] = (cmcs_tar_M[1, 1] - cmcs_tar_M[1, 2]) / cmcs_tar_M[1, 1]
  } else {
    cmcs_tar_S[j] = cmcs_tar_S[j-1] * (cmcs_tar_M[j, 1] - cmcs_tar_M[j, 2]) / cmcs_tar_M[j, 1]
  }
}

cmcs_tar_x = seq(1, cmcs_tar_time[1])
cmcs_tar_y = rep(1, cmcs_tar_time[1])

for (i in 1:(cmcs_tar_tk-1)) {
  cmcs_tar_x = c(cmcs_tar_x, seq(cmcs_tar_time[i], cmcs_tar_time[i+1]))
  cmcs_tar_y = c(cmcs_tar_y, rep(cmcs_tar_S[i], cmcs_tar_time[i+1]-cmcs_tar_time[i]+1))
}

cmcs_tar_x = c(cmcs_tar_x, cmcs_tar_time[cmcs_tar_tk])
cmcs_tar_y = c(cmcs_tar_y, cmcs_tar_S[cmcs_tar_tk])

### Plot

plot(cmcs_comp_x, cmcs_comp_y, type = 'l', xlim = c(1, cmcs_timemax), 
     ylim = c(min(cmcs_comp_y, cmcs_tar_y) - 0.01, 1), lwd = 3, main = "CMCS(Original)")
points(cmcs_tar_x, cmcs_tar_y, type = 'l', lwd = 3, col ='red')
legend("topright", legend=c("Comparator", "Target"), col = c('black', 'red'), lwd = 3)

################ DP

set.seed(3)
S = 2 ## L1 Sensitivity
eps = c(1, 2, 3)

### eps = 1

### DP-comp(1)

cmcs_comp_Mp_1 = cmcs_comp_M

cmcs_comp_noise_1 = rexp(nrow(cmcs_comp_M)*2+1, eps[1]/S) * sign(runif(nrow(cmcs_comp_M)*2+1, -1, 1))
cmcs_comp_noise_1 = pmax(cmcs_comp_noise_1, rep(0, nrow(cmcs_comp_M)*2+1))

cmcs_comp_Mp_1[1,1] = cmcs_comp_Mp_1[1,1] + cmcs_comp_noise_1[1]
cmcs_comp_Mp_1[,2] = cmcs_comp_M[,2] + cmcs_comp_noise_1[2:(nrow(cmcs_comp_M)+1)]
## cmcs_comp_Mp_1[,3] = cmcs_comp_M[,3] + cmcs_comp_noise_1[-(1:(nrow(cmcs_comp_M)+1))]

for (i in 2:nrow(cmcs_comp_Mp_1)) {
  cmcs_comp_Mp_1[i,1] = cmcs_comp_Mp_1[i-1,1] - cmcs_comp_Mp_1[i-1,2] - cmcs_comp_Mp_1[i-1,3]
}

## Excluding negative at-risk
cmcs_comp_Mp_1 = cmcs_comp_Mp_1[which(cmcs_comp_Mp_1[,1] > 0),]  
cmcs_comp_tk_1 = sum(cmcs_comp_Mp_1[,1] > 0)

## Survival function
cmcs_comp_dp_1_S = rep(NA, cmcs_comp_tk_1) 

for (i in cmcs_comp_time[1:cmcs_comp_tk_1]) {
  j = which(cmcs_comp_time == i)
  if (j == 1) {
    cmcs_comp_dp_1_S[j] = (cmcs_comp_Mp_1[1, 1] - cmcs_comp_Mp_1[1, 2]) / cmcs_comp_Mp_1[1, 1]
  } else {
    cmcs_comp_dp_1_S[j] = cmcs_comp_dp_1_S[j-1] * (cmcs_comp_Mp_1[j, 1] - cmcs_comp_Mp_1[j, 2]) / cmcs_comp_Mp_1[j, 1]
  }
}

cmcs_comp_dp_1_x = seq(1, cmcs_comp_time[1])
cmcs_comp_dp_1_y = rep(1, cmcs_comp_time[1])

for (i in 1:(cmcs_comp_tk_1-1)) {
  cmcs_comp_dp_1_x = c(cmcs_comp_dp_1_x, seq(cmcs_comp_time[i], cmcs_comp_time[i+1]))
  cmcs_comp_dp_1_y = c(cmcs_comp_dp_1_y, rep(cmcs_comp_dp_1_S[i], cmcs_comp_time[i+1]-cmcs_comp_time[i]+1))
}

cmcs_comp_dp_1_x = c(cmcs_comp_dp_1_x, cmcs_comp_time[cmcs_comp_tk_1])
cmcs_comp_dp_1_y = c(cmcs_comp_dp_1_y, cmcs_comp_dp_1_S[cmcs_comp_tk_1])


### DP-tar(1)

cmcs_tar_Mp_1 = cmcs_tar_M

cmcs_tar_noise_1 = rexp(nrow(cmcs_tar_M)*2+1, eps[1]/S) * sign(runif(nrow(cmcs_tar_M)*2+1, -1, 1))
cmcs_tar_noise_1 = pmax(cmcs_tar_noise_1, rep(0, nrow(cmcs_tar_M)*2+1))

cmcs_tar_Mp_1[1,1] = cmcs_tar_Mp_1[1,1] + cmcs_tar_noise_1[1]
cmcs_tar_Mp_1[,2] = cmcs_tar_M[,2] + cmcs_tar_noise_1[2:(nrow(cmcs_tar_M)+1)]
## cmcs_tar_Mp_1[,3] = cmcs_tar_M[,3] + cmcs_tar_noise_1[-(1:(nrow(cmcs_tar_M)+1))]

for (i in 2:nrow(cmcs_tar_Mp_1)) {
  cmcs_tar_Mp_1[i,1] = cmcs_tar_Mp_1[i-1,1] - cmcs_tar_Mp_1[i-1,2] - cmcs_tar_Mp_1[i-1,3]
}

## Excluding negative at-risk
cmcs_tar_Mp_1 = cmcs_tar_Mp_1[which(cmcs_tar_Mp_1[,1] > 0),]  
cmcs_tar_tk_1 = sum(cmcs_tar_Mp_1[,1] > 0)

## Survival function
cmcs_tar_dp_1_S = rep(NA, cmcs_tar_tk_1) 

for (i in cmcs_tar_time[1:cmcs_tar_tk_1]) {
  j = which(cmcs_tar_time == i)
  if (j == 1) {
    cmcs_tar_dp_1_S[j] = (cmcs_tar_Mp_1[1, 1] - cmcs_tar_Mp_1[1, 2]) / cmcs_tar_Mp_1[1, 1]
  } else {
    cmcs_tar_dp_1_S[j] = cmcs_tar_dp_1_S[j-1] * (cmcs_tar_Mp_1[j, 1] - cmcs_tar_Mp_1[j, 2]) / cmcs_tar_Mp_1[j, 1]
  }
}

cmcs_tar_dp_1_x = seq(1, cmcs_tar_time[1])
cmcs_tar_dp_1_y = rep(1, cmcs_tar_time[1])

for (i in 1:(cmcs_tar_tk_1-1)) {
  cmcs_tar_dp_1_x = c(cmcs_tar_dp_1_x, seq(cmcs_tar_time[i], cmcs_tar_time[i+1]))
  cmcs_tar_dp_1_y = c(cmcs_tar_dp_1_y, rep(cmcs_tar_dp_1_S[i], cmcs_tar_time[i+1]-cmcs_tar_time[i]+1))
}

cmcs_tar_dp_1_x = c(cmcs_tar_dp_1_x, cmcs_tar_time[cmcs_tar_tk_1])
cmcs_tar_dp_1_y = c(cmcs_tar_dp_1_y, cmcs_tar_dp_1_S[cmcs_tar_tk_1])


plot(cmcs_comp_dp_1_x, cmcs_comp_dp_1_y, type = 'l', xlim = c(1, cmcs_timemax), 
     ylim = c(min(cmcs_comp_dp_1_y, cmcs_tar_dp_1_y) - 0.01, 1), lwd = 3, main = "CMCS(DP, eps=1)")
points(cmcs_tar_dp_1_x, cmcs_tar_dp_1_y, type = 'l', lwd = 3, col ='red')
legend("topright", legend = c("Comparator", "Target"), col = c('black', 'red'), lwd = 3)


### eps = 2

### DP-comp(2)

## Constructing M'
cmcs_comp_Mp_2 = cmcs_comp_M

cmcs_comp_noise_2 = rexp(nrow(cmcs_comp_M)*2+1, eps[1]/S) * sign(runif(nrow(cmcs_comp_M)*2+1, -1, 1))
cmcs_comp_noise_2 = pmax(cmcs_comp_noise_2, rep(0, nrow(cmcs_comp_M)*2+1))

cmcs_comp_Mp_2[1,1] = cmcs_comp_Mp_2[1,1] + cmcs_comp_noise_2[1]
cmcs_comp_Mp_2[,2] = cmcs_comp_M[,2] + cmcs_comp_noise_2[2:(nrow(cmcs_comp_M)+1)]
# cmcs_comp_Mp_2[,3] = cmcs_comp_M[,3] + cmcs_comp_noise_2[-(1:(nrow(cmcs_comp_M)+1))]

for (i in 2:nrow(cmcs_comp_Mp_2)) {
  cmcs_comp_Mp_2[i,1] = cmcs_comp_Mp_2[i-1,1] - cmcs_comp_Mp_2[i-1,2] - cmcs_comp_Mp_2[i-1,3]
}

## Excluding negative at-risk
cmcs_comp_Mp_2 = cmcs_comp_Mp_2[which(cmcs_comp_Mp_2[,1] > 0),]  
cmcs_comp_tk_2 = sum(cmcs_comp_Mp_2[,1] > 0)

## Survival function
cmcs_comp_dp_2_S = rep(NA, cmcs_comp_tk_2)

for (i in cmcs_comp_time[1:cmcs_comp_tk_2]) {
  j = which(cmcs_comp_time == i)
  if (j == 1) {
    cmcs_comp_dp_2_S[j] = (cmcs_comp_Mp_2[1, 1] - cmcs_comp_Mp_2[1, 2]) / cmcs_comp_Mp_2[1, 1]
  } else {
    cmcs_comp_dp_2_S[j] = cmcs_comp_dp_2_S[j-1] * (cmcs_comp_Mp_2[j, 1] - cmcs_comp_Mp_2[j, 2]) / cmcs_comp_Mp_2[j, 1]
  }
}

cmcs_comp_dp_2_x = seq(1, cmcs_comp_time[1])
cmcs_comp_dp_2_y = rep(1, cmcs_comp_time[1])

for (i in 1:(cmcs_comp_tk_2-1)) {
  cmcs_comp_dp_2_x = c(cmcs_comp_dp_2_x, seq(cmcs_comp_time[i], cmcs_comp_time[i+1]))
  cmcs_comp_dp_2_y = c(cmcs_comp_dp_2_y, rep(cmcs_comp_dp_2_S[i], cmcs_comp_time[i+1]-cmcs_comp_time[i]+1))
}

cmcs_comp_dp_2_x = c(cmcs_comp_dp_2_x, cmcs_comp_time[cmcs_comp_tk_2])
cmcs_comp_dp_2_y = c(cmcs_comp_dp_2_y, cmcs_comp_dp_2_S[cmcs_comp_tk_2])

### DP-tar(2)

cmcs_tar_Mp_2 = cmcs_tar_M

cmcs_tar_noise_2 = rexp(nrow(cmcs_tar_M)*2+1, eps[1]/S) * sign(runif(nrow(cmcs_tar_M)*2+1, -1, 1))
cmcs_tar_noise_2 = pmax(cmcs_tar_noise_2, rep(0, nrow(cmcs_tar_M)*2+1))

cmcs_tar_Mp_2[1,1] = cmcs_tar_Mp_2[1,1] + cmcs_tar_noise_2[1]
cmcs_tar_Mp_2[,2] = cmcs_tar_M[,2] + cmcs_tar_noise_2[2:(nrow(cmcs_tar_M)+1)]
## cmcs_tar_Mp_2[,3] = cmcs_tar_M[,3] + cmcs_tar_noise_2[-(1:(nrow(cmcs_tar_M)+1))]

for (i in 2:nrow(cmcs_tar_Mp_2)) {
  cmcs_tar_Mp_2[i,1] = cmcs_tar_Mp_2[i-1,1] - cmcs_tar_Mp_2[i-1,2] - cmcs_tar_Mp_2[i-1,3]
}

## Excluding negative at-risk
cmcs_tar_Mp_2 = cmcs_tar_Mp_2[which(cmcs_tar_Mp_2[,1] > 0),]  
cmcs_tar_tk_2 = sum(cmcs_tar_Mp_2[,1] > 0)

## Survival function
cmcs_tar_dp_2_S = rep(NA, cmcs_tar_tk_2)

for (i in cmcs_tar_time[1:cmcs_tar_tk_2]) {
  j = which(cmcs_tar_time == i)
  if (j == 1) {
    cmcs_tar_dp_2_S[j] = (cmcs_tar_Mp_2[1, 1] - cmcs_tar_Mp_2[1, 2]) / cmcs_tar_Mp_2[1, 1]
  } else {
    cmcs_tar_dp_2_S[j] = cmcs_tar_dp_2_S[j-1] * (cmcs_tar_Mp_2[j, 1] - cmcs_tar_Mp_2[j, 2]) / cmcs_tar_Mp_2[j, 1]
  }
}

cmcs_tar_dp_2_x = seq(1, cmcs_tar_time[1])
cmcs_tar_dp_2_y = rep(1, cmcs_tar_time[1])

for (i in 1:(cmcs_tar_tk_2-1)) {
  cmcs_tar_dp_2_x = c(cmcs_tar_dp_2_x, seq(cmcs_tar_time[i], cmcs_tar_time[i+1]))
  cmcs_tar_dp_2_y = c(cmcs_tar_dp_2_y, rep(cmcs_tar_dp_2_S[i], cmcs_tar_time[i+1]-cmcs_tar_time[i]+1))
}

cmcs_tar_dp_2_x = c(cmcs_tar_dp_2_x, cmcs_tar_time[cmcs_tar_tk_2])
cmcs_tar_dp_2_y = c(cmcs_tar_dp_2_y, cmcs_tar_dp_2_S[cmcs_tar_tk_2])


plot(cmcs_comp_dp_2_x, cmcs_comp_dp_2_y, type = 'l', xlim = c(1, cmcs_timemax), 
     ylim = c(min(cmcs_comp_dp_2_y, cmcs_tar_dp_2_y) - 0.01, 1), lwd = 3, main = "CMCS(DP, eps=2)")
points(cmcs_tar_dp_2_x, cmcs_tar_dp_2_y, type = 'l', lwd = 3, col ='red')
legend("topright", legend = c("Comparator", "Target"), col = c('black', 'red'), lwd = 3)


### eps = 3

### DP-comp(3)

cmcs_comp_Mp_3 = cmcs_comp_M

cmcs_comp_noise_3 = rexp(nrow(cmcs_comp_M)*2+1, eps[1]/S) * sign(runif(nrow(cmcs_comp_M)*2+1, -1, 1))
cmcs_comp_noise_3 = pmax(cmcs_comp_noise_3, rep(0, nrow(cmcs_comp_M)*2+1))

cmcs_comp_Mp_3[1,1] = cmcs_comp_Mp_3[1,1] + cmcs_comp_noise_3[1]
cmcs_comp_Mp_3[,2] = cmcs_comp_M[,2] + cmcs_comp_noise_3[2:(nrow(cmcs_comp_M)+1)]
## cmcs_comp_Mp_3[,3] = cmcs_comp_M[,3] + cmcs_comp_noise_3[-(1:(nrow(cmcs_comp_M)+1))]

for (i in 2:nrow(cmcs_comp_Mp_3)) {
  cmcs_comp_Mp_3[i,1] = cmcs_comp_Mp_3[i-1,1] - cmcs_comp_Mp_3[i-1,2] - cmcs_comp_Mp_3[i-1,3]
}

## Excluding negative at-risk
cmcs_comp_Mp_3 = cmcs_comp_Mp_3[which(cmcs_comp_Mp_3[,1] > 0),]  
cmcs_comp_tk_3 = sum(cmcs_comp_Mp_3[,1] > 0)

## Survival function
cmcs_comp_dp_3_S = rep(NA, cmcs_comp_tk_3) 

for (i in cmcs_comp_time[1:cmcs_comp_tk_3]) {
  j = which(cmcs_comp_time == i)
  if (j == 1) {
    cmcs_comp_dp_3_S[j] = (cmcs_comp_Mp_3[1, 1] - cmcs_comp_Mp_3[1, 2]) / cmcs_comp_Mp_3[1, 1]
  } else {
    cmcs_comp_dp_3_S[j] = cmcs_comp_dp_3_S[j-1] * (cmcs_comp_Mp_3[j, 1] - cmcs_comp_Mp_3[j, 2]) / cmcs_comp_Mp_3[j, 1]
  }
}

cmcs_comp_dp_3_x = seq(1, cmcs_comp_time[1])
cmcs_comp_dp_3_y = rep(1, cmcs_comp_time[1])

for (i in 1:(cmcs_comp_tk_3-1)) {
  cmcs_comp_dp_3_x = c(cmcs_comp_dp_3_x, seq(cmcs_comp_time[i], cmcs_comp_time[i+1]))
  cmcs_comp_dp_3_y = c(cmcs_comp_dp_3_y, rep(cmcs_comp_dp_3_S[i], cmcs_comp_time[i+1]-cmcs_comp_time[i]+1))
}

cmcs_comp_dp_3_x = c(cmcs_comp_dp_3_x, cmcs_comp_time[cmcs_comp_tk_3])
cmcs_comp_dp_3_y = c(cmcs_comp_dp_3_y, cmcs_comp_dp_3_S[cmcs_comp_tk_3])


### DP-tar(3)

cmcs_tar_Mp_3 = cmcs_tar_M

cmcs_tar_noise_3 = rexp(nrow(cmcs_tar_M)*2+1, eps[1]/S) * sign(runif(nrow(cmcs_tar_M)*2+1, -1, 1))
cmcs_tar_noise_3 = pmax(cmcs_tar_noise_3, rep(0, nrow(cmcs_tar_M)*2+1))

cmcs_tar_Mp_3[1,1] = cmcs_tar_Mp_3[1,1] + cmcs_tar_noise_3[1]
cmcs_tar_Mp_3[,2] = cmcs_tar_M[,2] + cmcs_tar_noise_3[2:(nrow(cmcs_tar_M)+1)]
## cmcs_tar_Mp_3[,3] = cmcs_tar_M[,3] + cmcs_tar_noise_3[-(1:(nrow(cmcs_tar_M)+1))]

for (i in 2:nrow(cmcs_tar_Mp_3)) {
  cmcs_tar_Mp_3[i,1] = cmcs_tar_Mp_3[i-1,1] - cmcs_tar_Mp_3[i-1,2] - cmcs_tar_Mp_3[i-1,3]
}

## Excluding negative at-risk
cmcs_tar_Mp_3 = cmcs_tar_Mp_3[which(cmcs_tar_Mp_3[,1] > 0),]  
cmcs_tar_tk_3 = sum(cmcs_tar_Mp_3[,1] > 0)

## Survival function
cmcs_tar_dp_3_S = rep(NA, cmcs_tar_tk_3) 

for (i in cmcs_tar_time[1:cmcs_tar_tk_3]) {
  j = which(cmcs_tar_time == i)
  if (j == 1) {
    cmcs_tar_dp_3_S[j] = (cmcs_tar_Mp_3[1, 1] - cmcs_tar_Mp_3[1, 2]) / cmcs_tar_Mp_3[1, 1]
  } else {
    cmcs_tar_dp_3_S[j] = cmcs_tar_dp_3_S[j-1] * (cmcs_tar_Mp_3[j, 1] - cmcs_tar_Mp_3[j, 2]) / cmcs_tar_Mp_3[j, 1]
  }
}

cmcs_tar_dp_3_x = seq(1, cmcs_tar_time[1])
cmcs_tar_dp_3_y = rep(1, cmcs_tar_time[1])

for (i in 1:(cmcs_tar_tk_3-1)) {
  cmcs_tar_dp_3_x = c(cmcs_tar_dp_3_x, seq(cmcs_tar_time[i], cmcs_tar_time[i+1]))
  cmcs_tar_dp_3_y = c(cmcs_tar_dp_3_y, rep(cmcs_tar_dp_3_S[i], cmcs_tar_time[i+1]-cmcs_tar_time[i]+1))
}

cmcs_tar_dp_3_x = c(cmcs_tar_dp_3_x, cmcs_tar_time[cmcs_tar_tk_3])
cmcs_tar_dp_3_y = c(cmcs_tar_dp_3_y, cmcs_tar_dp_3_S[cmcs_tar_tk_3])


plot(cmcs_comp_dp_3_x, cmcs_comp_dp_3_y, type = 'l', xlim = c(1, cmcs_timemax), 
     ylim = c(min(cmcs_comp_dp_3_y, cmcs_tar_dp_3_y) - 0.01, 1), lwd = 3, main = "CMCS(DP, eps=3)")
points(cmcs_tar_dp_3_x, cmcs_tar_dp_3_y, type = 'l', lwd = 3, col ='red')
legend("topright", legend = c("Comparator", "Target"), col = c('black', 'red'), lwd = 3)


### Compare by eps

# Comparator
plot(cmcs_comp_x, cmcs_comp_y, type = 'l', xlim = c(1, cmcs_timemax), 
     ylim = c(min(cmcs_comp_y, cmcs_comp_dp_1_y, cmcs_comp_dp_2_y, cmcs_comp_dp_3_y) - 0.01, 1),
     lwd = 3, main = "CMCS(Comparator, by eps)")
points(cmcs_comp_dp_1_x, cmcs_comp_dp_1_y, type = 'l', lwd = 3, col ='blue')
points(cmcs_comp_dp_2_x, cmcs_comp_dp_2_y, type = 'l', lwd = 3, col ='purple')
points(cmcs_comp_dp_3_x, cmcs_comp_dp_3_y, type = 'l', lwd = 3, col ='green')
legend("topright", legend = c("NA", "1", "2", "3"), col = c('black', 'blue', 'purple', 'green'), lwd = 3)

# Target
plot(cmcs_tar_x, cmcs_tar_y, type = 'l', xlim = c(1, cmcs_timemax), 
     ylim = c(min(cmcs_tar_y, cmcs_tar_dp_1_y, cmcs_tar_dp_2_y, cmcs_tar_dp_3_y) - 0.01, 1),
     lwd = 3, main = "CMCS(Target, by eps)")
points(cmcs_tar_dp_1_x, cmcs_tar_dp_1_y, type = 'l', lwd = 3, col ='blue')
points(cmcs_tar_dp_2_x, cmcs_tar_dp_2_y, type = 'l', lwd = 3, col ='purple')
points(cmcs_tar_dp_3_x, cmcs_tar_dp_3_y, type = 'l', lwd = 3, col ='green')
legend("topright", legend = c("NA", "1", "2", "3"), col = c('black', 'blue', 'purple', 'green'), lwd = 3)





