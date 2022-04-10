
eps_p_1 = matrix(rep(NA, 30), ncol = 3)
colnames(eps_p_1) = c("SNUBH", "SNUH", "CMCS")
row.names(eps_p_1) = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
eps_p_2 = matrix(rep(NA, 30), ncol = 3)
colnames(eps_p_2) = c("SNUBH", "SNUH", "CMCS")
row.names(eps_p_2) = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
eps_p_3 = matrix(rep(NA, 30), ncol = 3)
colnames(eps_p_3) = c("SNUBH", "SNUH", "CMCS")
row.names(eps_p_3) = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

true_p = rep(NA, 3)

S = 2
eps = c(1,2,3)

### SNUBH

snubh_M = unique(sort(mtsnubh$survivalTime[mtsnubh$outcomeCount != 0]))
snubh_M = cbind(snubh_M, matrix(rep(NA, length(snubh_M)*9), ncol=9))
colnames(snubh_M) = c("Time", "r1", "d1", "c1", "r2", "d2", "c2", "E1", "E2", "V")

for (i in 1:nrow(snubh_M)) {
  snubh_M[i,2] = sum(mtsnubh$survivalTime[which(mtsnubh$treatment == 0)] >= snubh_M[i,1])
  snubh_M[i,5] = sum(mtsnubh$survivalTime[which(mtsnubh$treatment == 1)] >= snubh_M[i,1])
  snubh_M[i,3] = sum(mtsnubh$survivalTime[intersect(which(mtsnubh$treatment == 0), which(mtsnubh$outcomeCount != 0))] == snubh_M[i,1])
  snubh_M[i,6] = sum(mtsnubh$survivalTime[intersect(which(mtsnubh$treatment == 1), which(mtsnubh$outcomeCount != 0))] == snubh_M[i,1])
  snubh_M[i,8] = (snubh_M[i,3] + snubh_M[i,6]) * snubh_M[i,2] / (snubh_M[i,2] + snubh_M[i,5])
  snubh_M[i,9] = (snubh_M[i,3] + snubh_M[i,6]) * snubh_M[i,5] / (snubh_M[i,2] + snubh_M[i,5])
  snubh_M[i,10] = snubh_M[i,2] * snubh_M[i,5] * (snubh_M[i,3] + snubh_M[i,6]) *
    (snubh_M[i,2] + snubh_M[i,5] - (snubh_M[i,3] + snubh_M[i,6])) / ((snubh_M[i,2] + snubh_M[i,5])^2 * (snubh_M[i,2] + snubh_M[i,5] - 1))
  if(i >= 2) {
    snubh_M[i-1,4] = snubh_M[i-1,2] - snubh_M[i,2] - snubh_M[i,3]
    snubh_M[i-1,7] = snubh_M[i-1,5] - snubh_M[i,5] - snubh_M[i,6]
  }
}

snubh_Z = sum(snubh_M[,3] - snubh_M[,8])^2 / sum(snubh_M[,10])
snubh_Z_p = 1 - pchisq(snubh_Z, 1)

true_p[1] = snubh_Z_p

for(j in 1:10) {
  set.seed(j)
  
  ### eps = 1
  
  snubh_Mp_1 = snubh_M
  # snubh_noise_1 = rexp(nrow(snubh_M)*2+2, eps[1]/S) * sign(runif(nrow(snubh_M)*2+2, -1, 1))
  snubh_noise_1 = rexp(nrow(snubh_M)*4+2, eps[1]/S) * sign(runif(nrow(snubh_M)*4+2, -1, 1))
  snubh_noise_1 = pmax(snubh_noise_1, 0)
  
  snubh_Mp_1[1,2] = snubh_Mp_1[1,2] + snubh_noise_1[1]
  snubh_Mp_1[1,5] = snubh_Mp_1[1,5] + snubh_noise_1[2]
  #snubh_Mp_1[,3] = snubh_Mp_1[,3] + snubh_noise_1[3:(nrow(snubh_M)+2)]
  #snubh_Mp_1[,6] = snubh_Mp_1[,6] + snubh_noise_1[-c(1:(nrow(snubh_M)+2))]
  
  snubh_Mp_1[,3] = snubh_Mp_1[,3] + snubh_noise_1[3:(nrow(snubh_M)+2)]
  snubh_Mp_1[,6] = snubh_Mp_1[,6] + snubh_noise_1[(nrow(snubh_M)+3):(nrow(snubh_M)*2+2)]
  snubh_Mp_1[,4] = snubh_Mp_1[,4] + snubh_noise_1[(nrow(snubh_M)*2+3):(nrow(snubh_M)*3+2)]
  snubh_Mp_1[,7] = snubh_Mp_1[,7] + snubh_noise_1[(nrow(snubh_M)*3+3):(nrow(snubh_M)*4+2)]
  
  
  for(i in 2:nrow(snubh_Mp_1)) {
    snubh_Mp_1[i,2] = snubh_Mp_1[i-1,2] - snubh_Mp_1[i-1,3] - snubh_Mp_1[i-1,4]
    snubh_Mp_1[i,5] = snubh_Mp_1[i-1,5] - snubh_Mp_1[i-1,6] - snubh_Mp_1[i-1,7]
  }
  
  snubh_Mp_1[,8] = (snubh_Mp_1[,3] + snubh_Mp_1[,6]) * snubh_Mp_1[,2] / (snubh_Mp_1[,2] + snubh_Mp_1[,5])
  snubh_Mp_1[,9] = (snubh_Mp_1[,3] + snubh_Mp_1[,6]) * snubh_Mp_1[,5] / (snubh_Mp_1[,2] + snubh_Mp_1[,5])
  
  snubh_Mp_1_lim = min(max(which(snubh_Mp_1[,2] > 0)), max(which(snubh_Mp_1[,5] > 0)))
  
  snubh_Mp_1 = snubh_Mp_1[1:snubh_Mp_1_lim,]
  
  snubh_Z_1 = sum(snubh_Mp_1[,3] - snubh_Mp_1[,8])^2 / sum(snubh_Mp_1[,10])
  snubh_Z_1_p = 1 - pchisq(snubh_Z_1, 1)
  
  eps_p_1[j,1] = snubh_Z_1_p
  
  ### eps = 2
  
  snubh_Mp_2 = snubh_M
  # snubh_noise_2 = rexp(nrow(snubh_M)*2+2, eps[1]/S) * sign(runif(nrow(snubh_M)*2+2, -1, 1))
  snubh_noise_2 = rexp(nrow(snubh_M)*4+2, eps[1]/S) * sign(runif(nrow(snubh_M)*4+2, -1, 1))
  snubh_noise_2 = pmax(snubh_noise_2, 0)
  
  snubh_Mp_2[1,2] = snubh_Mp_2[1,2] + snubh_noise_2[1]
  snubh_Mp_2[1,5] = snubh_Mp_2[1,5] + snubh_noise_2[2]
  #snubh_Mp_2[,3] = snubh_Mp_2[,3] + snubh_noise_2[3:(nrow(snubh_M)+2)]
  #snubh_Mp_2[,6] = snubh_Mp_2[,6] + snubh_noise_2[-c(1:(nrow(snubh_M)+2))]
  
  snubh_Mp_2[,3] = snubh_Mp_2[,3] + snubh_noise_2[3:(nrow(snubh_M)+2)]
  snubh_Mp_2[,6] = snubh_Mp_2[,6] + snubh_noise_2[(nrow(snubh_M)+3):(nrow(snubh_M)*2+2)]
  snubh_Mp_2[,4] = snubh_Mp_2[,4] + snubh_noise_2[(nrow(snubh_M)*2+3):(nrow(snubh_M)*3+2)]
  snubh_Mp_2[,7] = snubh_Mp_2[,7] + snubh_noise_2[(nrow(snubh_M)*3+3):(nrow(snubh_M)*4+2)]
  
  
  for(i in 2:nrow(snubh_Mp_2)) {
    snubh_Mp_2[i,2] = snubh_Mp_2[i-1,2] - snubh_Mp_2[i-1,3] - snubh_Mp_2[i-1,4]
    snubh_Mp_2[i,5] = snubh_Mp_2[i-1,5] - snubh_Mp_2[i-1,6] - snubh_Mp_2[i-1,7]
  }
  
  snubh_Mp_2[,8] = (snubh_Mp_2[,3] + snubh_Mp_2[,6]) * snubh_Mp_2[,2] / (snubh_Mp_2[,2] + snubh_Mp_2[,5])
  snubh_Mp_2[,9] = (snubh_Mp_2[,3] + snubh_Mp_2[,6]) * snubh_Mp_2[,5] / (snubh_Mp_2[,2] + snubh_Mp_2[,5])
  
  snubh_Mp_2_lim = min(max(which(snubh_Mp_2[,2] > 0)), max(which(snubh_Mp_2[,5] > 0)))
  
  snubh_Mp_2 = snubh_Mp_2[1:snubh_Mp_2_lim,]
  
  snubh_Z_2 = sum(snubh_Mp_2[,3] - snubh_Mp_2[,8])^2 / sum(snubh_Mp_2[,10])
  snubh_Z_2_p = 1 - pchisq(snubh_Z_2, 1)
  
  eps_p_2[j,1] = snubh_Z_2_p
  
  ### eps = 3
  
  snubh_Mp_3 = snubh_M
  # snubh_noise_3 = rexp(nrow(snubh_M)*2+2, eps[1]/S) * sign(runif(nrow(snubh_M)*2+2, -1, 1))
  snubh_noise_3 = rexp(nrow(snubh_M)*4+2, eps[1]/S) * sign(runif(nrow(snubh_M)*4+2, -1, 1))
  snubh_noise_3 = pmax(snubh_noise_3, 0)
  
  snubh_Mp_3[1,2] = snubh_Mp_3[1,2] + snubh_noise_3[1]
  snubh_Mp_3[1,5] = snubh_Mp_3[1,5] + snubh_noise_3[2]
  #snubh_Mp_3[,3] = snubh_Mp_3[,3] + snubh_noise_3[3:(nrow(snubh_M)+2)]
  #snubh_Mp_3[,6] = snubh_Mp_3[,6] + snubh_noise_3[-c(1:(nrow(snubh_M)+2))]
  
  snubh_Mp_3[,3] = snubh_Mp_3[,3] + snubh_noise_3[3:(nrow(snubh_M)+2)]
  snubh_Mp_3[,6] = snubh_Mp_3[,6] + snubh_noise_3[(nrow(snubh_M)+3):(nrow(snubh_M)*2+2)]
  snubh_Mp_3[,4] = snubh_Mp_3[,4] + snubh_noise_3[(nrow(snubh_M)*2+3):(nrow(snubh_M)*3+2)]
  snubh_Mp_3[,7] = snubh_Mp_3[,7] + snubh_noise_3[(nrow(snubh_M)*3+3):(nrow(snubh_M)*4+2)]
  
  
  for(i in 2:nrow(snubh_Mp_3)) {
    snubh_Mp_3[i,2] = snubh_Mp_3[i-1,2] - snubh_Mp_3[i-1,3] - snubh_Mp_3[i-1,4]
    snubh_Mp_3[i,5] = snubh_Mp_3[i-1,5] - snubh_Mp_3[i-1,6] - snubh_Mp_3[i-1,7]
  }
  
  snubh_Mp_3[,8] = (snubh_Mp_3[,3] + snubh_Mp_3[,6]) * snubh_Mp_3[,2] / (snubh_Mp_3[,2] + snubh_Mp_3[,5])
  snubh_Mp_3[,9] = (snubh_Mp_3[,3] + snubh_Mp_3[,6]) * snubh_Mp_3[,5] / (snubh_Mp_3[,2] + snubh_Mp_3[,5])
  
  snubh_Mp_3_lim = min(max(which(snubh_Mp_3[,2] > 0)), max(which(snubh_Mp_3[,5] > 0)))
  
  snubh_Mp_3 = snubh_Mp_3[1:snubh_Mp_3_lim,]
  
  snubh_Z_3 = sum(snubh_Mp_3[,3] - snubh_Mp_3[,8])^2 / sum(snubh_Mp_3[,10])
  snubh_Z_3_p = 1 - pchisq(snubh_Z_3, 1)
  
  eps_p_3[j,1] = snubh_Z_3_p
}


### SNUH

snuh_M = unique(sort(mtsnuh$survivalTime[mtsnuh$outcomeCount != 0]))
snuh_M = cbind(snuh_M, matrix(rep(NA, length(snuh_M)*9), ncol=9))
colnames(snuh_M) = c("Time", "r1", "d1", "c1", "r2", "d2", "c2", "E1", "E2", "V")

for (i in 1:nrow(snuh_M)) {
  snuh_M[i,2] = sum(mtsnuh$survivalTime[which(mtsnuh$treatment == 0)] >= snuh_M[i,1])
  snuh_M[i,5] = sum(mtsnuh$survivalTime[which(mtsnuh$treatment == 1)] >= snuh_M[i,1])
  snuh_M[i,3] = sum(mtsnuh$survivalTime[intersect(which(mtsnuh$treatment == 0), which(mtsnuh$outcomeCount != 0))] == snuh_M[i,1])
  snuh_M[i,6] = sum(mtsnuh$survivalTime[intersect(which(mtsnuh$treatment == 1), which(mtsnuh$outcomeCount != 0))] == snuh_M[i,1])
  snuh_M[i,8] = (snuh_M[i,3] + snuh_M[i,6]) * snuh_M[i,2] / (snuh_M[i,2] + snuh_M[i,5])
  snuh_M[i,9] = (snuh_M[i,3] + snuh_M[i,6]) * snuh_M[i,5] / (snuh_M[i,2] + snuh_M[i,5])
  snuh_M[i,10] = snuh_M[i,2] * snuh_M[i,5] * (snuh_M[i,3] + snuh_M[i,6]) *
    (snuh_M[i,2] + snuh_M[i,5] - (snuh_M[i,3] + snuh_M[i,6])) / ((snuh_M[i,2] + snuh_M[i,5])^2 * (snuh_M[i,2] + snuh_M[i,5] - 1))
  if(i >= 2) {
    snuh_M[i-1,4] = snuh_M[i-1,2] - snuh_M[i,2] - snuh_M[i,3]
    snuh_M[i-1,7] = snuh_M[i-1,5] - snuh_M[i,5] - snuh_M[i,6]
  }
}

snuh_Z = sum(snuh_M[,3] - snuh_M[,8])^2 / sum(snuh_M[,10])
snuh_Z_p = 1 - pchisq(snuh_Z, 1)

true_p[2] = snuh_Z_p

for(j in 11:20) {
  set.seed(j)
  
  ### eps = 1
  
  snuh_Mp_1 = snuh_M
  # snuh_noise_1 = rexp(nrow(snuh_M)*2+2, eps[1]/S) * sign(runif(nrow(snuh_M)*2+2, -1, 1))
  snuh_noise_1 = rexp(nrow(snuh_M)*4+2, eps[1]/S) * sign(runif(nrow(snuh_M)*4+2, -1, 1))
  snuh_noise_1 = pmax(snuh_noise_1, 0)
  
  snuh_Mp_1[1,2] = snuh_Mp_1[1,2] + snuh_noise_1[1]
  snuh_Mp_1[1,5] = snuh_Mp_1[1,5] + snuh_noise_1[2]
  #snuh_Mp_1[,3] = snuh_Mp_1[,3] + snuh_noise_1[3:(nrow(snuh_M)+2)]
  #snuh_Mp_1[,6] = snuh_Mp_1[,6] + snuh_noise_1[-c(1:(nrow(snuh_M)+2))]
  
  snuh_Mp_1[,3] = snuh_Mp_1[,3] + snuh_noise_1[3:(nrow(snuh_M)+2)]
  snuh_Mp_1[,6] = snuh_Mp_1[,6] + snuh_noise_1[(nrow(snuh_M)+3):(nrow(snuh_M)*2+2)]
  snuh_Mp_1[,4] = snuh_Mp_1[,4] + snuh_noise_1[(nrow(snuh_M)*2+3):(nrow(snuh_M)*3+2)]
  snuh_Mp_1[,7] = snuh_Mp_1[,7] + snuh_noise_1[(nrow(snuh_M)*3+3):(nrow(snuh_M)*4+2)]
  
  
  for(i in 2:nrow(snuh_Mp_1)) {
    snuh_Mp_1[i,2] = snuh_Mp_1[i-1,2] - snuh_Mp_1[i-1,3] - snuh_Mp_1[i-1,4]
    snuh_Mp_1[i,5] = snuh_Mp_1[i-1,5] - snuh_Mp_1[i-1,6] - snuh_Mp_1[i-1,7]
  }
  
  snuh_Mp_1[,8] = (snuh_Mp_1[,3] + snuh_Mp_1[,6]) * snuh_Mp_1[,2] / (snuh_Mp_1[,2] + snuh_Mp_1[,5])
  snuh_Mp_1[,9] = (snuh_Mp_1[,3] + snuh_Mp_1[,6]) * snuh_Mp_1[,5] / (snuh_Mp_1[,2] + snuh_Mp_1[,5])
  
  snuh_Mp_1_lim = min(max(which(snuh_Mp_1[,2] > 0)), max(which(snuh_Mp_1[,5] > 0)))
  
  snuh_Mp_1 = snuh_Mp_1[1:snuh_Mp_1_lim,]
  
  snuh_Z_1 = sum(snuh_Mp_1[,3] - snuh_Mp_1[,8])^2 / sum(snuh_Mp_1[,10])
  snuh_Z_1_p = 1 - pchisq(snuh_Z_1, 1)
  
  eps_p_1[j-10,2] = snuh_Z_1_p
  
  ### eps = 2
  
  snuh_Mp_2 = snuh_M
  # snuh_noise_2 = rexp(nrow(snuh_M)*2+2, eps[1]/S) * sign(runif(nrow(snuh_M)*2+2, -1, 1))
  snuh_noise_2 = rexp(nrow(snuh_M)*4+2, eps[1]/S) * sign(runif(nrow(snuh_M)*4+2, -1, 1))
  snuh_noise_2 = pmax(snuh_noise_2, 0)
  
  snuh_Mp_2[1,2] = snuh_Mp_2[1,2] + snuh_noise_2[1]
  snuh_Mp_2[1,5] = snuh_Mp_2[1,5] + snuh_noise_2[2]
  #snuh_Mp_2[,3] = snuh_Mp_2[,3] + snuh_noise_2[3:(nrow(snuh_M)+2)]
  #snuh_Mp_2[,6] = snuh_Mp_2[,6] + snuh_noise_2[-c(1:(nrow(snuh_M)+2))]
  
  snuh_Mp_2[,3] = snuh_Mp_2[,3] + snuh_noise_2[3:(nrow(snuh_M)+2)]
  snuh_Mp_2[,6] = snuh_Mp_2[,6] + snuh_noise_2[(nrow(snuh_M)+3):(nrow(snuh_M)*2+2)]
  snuh_Mp_2[,4] = snuh_Mp_2[,4] + snuh_noise_2[(nrow(snuh_M)*2+3):(nrow(snuh_M)*3+2)]
  snuh_Mp_2[,7] = snuh_Mp_2[,7] + snuh_noise_2[(nrow(snuh_M)*3+3):(nrow(snuh_M)*4+2)]
  
  
  for(i in 2:nrow(snuh_Mp_2)) {
    snuh_Mp_2[i,2] = snuh_Mp_2[i-1,2] - snuh_Mp_2[i-1,3] - snuh_Mp_2[i-1,4]
    snuh_Mp_2[i,5] = snuh_Mp_2[i-1,5] - snuh_Mp_2[i-1,6] - snuh_Mp_2[i-1,7]
  }
  
  snuh_Mp_2[,8] = (snuh_Mp_2[,3] + snuh_Mp_2[,6]) * snuh_Mp_2[,2] / (snuh_Mp_2[,2] + snuh_Mp_2[,5])
  snuh_Mp_2[,9] = (snuh_Mp_2[,3] + snuh_Mp_2[,6]) * snuh_Mp_2[,5] / (snuh_Mp_2[,2] + snuh_Mp_2[,5])
  
  snuh_Mp_2_lim = min(max(which(snuh_Mp_2[,2] > 0)), max(which(snuh_Mp_2[,5] > 0)))
  
  snuh_Mp_2 = snuh_Mp_2[1:snuh_Mp_2_lim,]
  
  snuh_Z_2 = sum(snuh_Mp_2[,3] - snuh_Mp_2[,8])^2 / sum(snuh_Mp_2[,10])
  snuh_Z_2_p = 1 - pchisq(snuh_Z_2, 1)
  
  eps_p_2[j-10,2] = snuh_Z_2_p
  
  ### eps = 3
  
  snuh_Mp_3 = snuh_M
  # snuh_noise_3 = rexp(nrow(snuh_M)*2+2, eps[1]/S) * sign(runif(nrow(snuh_M)*2+2, -1, 1))
  snuh_noise_3 = rexp(nrow(snuh_M)*4+2, eps[1]/S) * sign(runif(nrow(snuh_M)*4+2, -1, 1))
  snuh_noise_3 = pmax(snuh_noise_3, 0)
  
  snuh_Mp_3[1,2] = snuh_Mp_3[1,2] + snuh_noise_3[1]
  snuh_Mp_3[1,5] = snuh_Mp_3[1,5] + snuh_noise_3[2]
  #snuh_Mp_3[,3] = snuh_Mp_3[,3] + snuh_noise_3[3:(nrow(snuh_M)+2)]
  #snuh_Mp_3[,6] = snuh_Mp_3[,6] + snuh_noise_3[-c(1:(nrow(snuh_M)+2))]
  
  snuh_Mp_3[,3] = snuh_Mp_3[,3] + snuh_noise_3[3:(nrow(snuh_M)+2)]
  snuh_Mp_3[,6] = snuh_Mp_3[,6] + snuh_noise_3[(nrow(snuh_M)+3):(nrow(snuh_M)*2+2)]
  snuh_Mp_3[,4] = snuh_Mp_3[,4] + snuh_noise_3[(nrow(snuh_M)*2+3):(nrow(snuh_M)*3+2)]
  snuh_Mp_3[,7] = snuh_Mp_3[,7] + snuh_noise_3[(nrow(snuh_M)*3+3):(nrow(snuh_M)*4+2)]
  
  
  for(i in 2:nrow(snuh_Mp_3)) {
    snuh_Mp_3[i,2] = snuh_Mp_3[i-1,2] - snuh_Mp_3[i-1,3] - snuh_Mp_3[i-1,4]
    snuh_Mp_3[i,5] = snuh_Mp_3[i-1,5] - snuh_Mp_3[i-1,6] - snuh_Mp_3[i-1,7]
  }
  
  snuh_Mp_3[,8] = (snuh_Mp_3[,3] + snuh_Mp_3[,6]) * snuh_Mp_3[,2] / (snuh_Mp_3[,2] + snuh_Mp_3[,5])
  snuh_Mp_3[,9] = (snuh_Mp_3[,3] + snuh_Mp_3[,6]) * snuh_Mp_3[,5] / (snuh_Mp_3[,2] + snuh_Mp_3[,5])
  
  snuh_Mp_3_lim = min(max(which(snuh_Mp_3[,2] > 0)), max(which(snuh_Mp_3[,5] > 0)))
  
  snuh_Mp_3 = snuh_Mp_3[1:snuh_Mp_3_lim,]
  
  snuh_Z_3 = sum(snuh_Mp_3[,3] - snuh_Mp_3[,8])^2 / sum(snuh_Mp_3[,10])
  snuh_Z_3_p = 1 - pchisq(snuh_Z_3, 1)
  
  eps_p_3[j-10,2] = snuh_Z_3_p
}



### CMCS

cmcs_M = unique(sort(mtcmcs$survivalTime[mtcmcs$outcomeCount != 0]))
cmcs_M = cbind(cmcs_M, matrix(rep(NA, length(cmcs_M)*9), ncol=9))
colnames(cmcs_M) = c("Time", "r1", "d1", "c1", "r2", "d2", "c2", "E1", "E2", "V")

for (i in 1:nrow(cmcs_M)) {
  cmcs_M[i,2] = sum(mtcmcs$survivalTime[which(mtcmcs$treatment == 0)] >= cmcs_M[i,1])
  cmcs_M[i,5] = sum(mtcmcs$survivalTime[which(mtcmcs$treatment == 1)] >= cmcs_M[i,1])
  cmcs_M[i,3] = sum(mtcmcs$survivalTime[intersect(which(mtcmcs$treatment == 0), which(mtcmcs$outcomeCount != 0))] == cmcs_M[i,1])
  cmcs_M[i,6] = sum(mtcmcs$survivalTime[intersect(which(mtcmcs$treatment == 1), which(mtcmcs$outcomeCount != 0))] == cmcs_M[i,1])
  cmcs_M[i,8] = (cmcs_M[i,3] + cmcs_M[i,6]) * cmcs_M[i,2] / (cmcs_M[i,2] + cmcs_M[i,5])
  cmcs_M[i,9] = (cmcs_M[i,3] + cmcs_M[i,6]) * cmcs_M[i,5] / (cmcs_M[i,2] + cmcs_M[i,5])
  cmcs_M[i,10] = cmcs_M[i,2] * cmcs_M[i,5] * (cmcs_M[i,3] + cmcs_M[i,6]) *
    (cmcs_M[i,2] + cmcs_M[i,5] - (cmcs_M[i,3] + cmcs_M[i,6])) / ((cmcs_M[i,2] + cmcs_M[i,5])^2 * (cmcs_M[i,2] + cmcs_M[i,5] - 1))
  if(i >= 2) {
    cmcs_M[i-1,4] = cmcs_M[i-1,2] - cmcs_M[i,2] - cmcs_M[i,3]
    cmcs_M[i-1,7] = cmcs_M[i-1,5] - cmcs_M[i,5] - cmcs_M[i,6]
  }
}

cmcs_Z = sum(cmcs_M[,3] - cmcs_M[,8])^2 / sum(cmcs_M[,10])
cmcs_Z_p = 1 - pchisq(cmcs_Z, 1)

true_p[3] = cmcs_Z_p

for(j in 21:30) {
  set.seed(j)
  
  ### eps = 1
  
  cmcs_Mp_1 = cmcs_M
  # cmcs_noise_1 = rexp(nrow(cmcs_M)*2+2, eps[1]/S) * sign(runif(nrow(cmcs_M)*2+2, -1, 1))
  cmcs_noise_1 = rexp(nrow(cmcs_M)*4+2, eps[1]/S) * sign(runif(nrow(cmcs_M)*4+2, -1, 1))
  cmcs_noise_1 = pmax(cmcs_noise_1, 0)
  
  cmcs_Mp_1[1,2] = cmcs_Mp_1[1,2] + cmcs_noise_1[1]
  cmcs_Mp_1[1,5] = cmcs_Mp_1[1,5] + cmcs_noise_1[2]
  #cmcs_Mp_1[,3] = cmcs_Mp_1[,3] + cmcs_noise_1[3:(nrow(cmcs_M)+2)]
  #cmcs_Mp_1[,6] = cmcs_Mp_1[,6] + cmcs_noise_1[-c(1:(nrow(cmcs_M)+2))]
  
  cmcs_Mp_1[,3] = cmcs_Mp_1[,3] + cmcs_noise_1[3:(nrow(cmcs_M)+2)]
  cmcs_Mp_1[,6] = cmcs_Mp_1[,6] + cmcs_noise_1[(nrow(cmcs_M)+3):(nrow(cmcs_M)*2+2)]
  cmcs_Mp_1[,4] = cmcs_Mp_1[,4] + cmcs_noise_1[(nrow(cmcs_M)*2+3):(nrow(cmcs_M)*3+2)]
  cmcs_Mp_1[,7] = cmcs_Mp_1[,7] + cmcs_noise_1[(nrow(cmcs_M)*3+3):(nrow(cmcs_M)*4+2)]
  
  
  for(i in 2:nrow(cmcs_Mp_1)) {
    cmcs_Mp_1[i,2] = cmcs_Mp_1[i-1,2] - cmcs_Mp_1[i-1,3] - cmcs_Mp_1[i-1,4]
    cmcs_Mp_1[i,5] = cmcs_Mp_1[i-1,5] - cmcs_Mp_1[i-1,6] - cmcs_Mp_1[i-1,7]
  }
  
  cmcs_Mp_1[,8] = (cmcs_Mp_1[,3] + cmcs_Mp_1[,6]) * cmcs_Mp_1[,2] / (cmcs_Mp_1[,2] + cmcs_Mp_1[,5])
  cmcs_Mp_1[,9] = (cmcs_Mp_1[,3] + cmcs_Mp_1[,6]) * cmcs_Mp_1[,5] / (cmcs_Mp_1[,2] + cmcs_Mp_1[,5])
  
  cmcs_Mp_1_lim = min(max(which(cmcs_Mp_1[,2] > 0)), max(which(cmcs_Mp_1[,5] > 0)))
  
  cmcs_Mp_1 = cmcs_Mp_1[1:cmcs_Mp_1_lim,]
  
  cmcs_Z_1 = sum(cmcs_Mp_1[,3] - cmcs_Mp_1[,8])^2 / sum(cmcs_Mp_1[,10])
  cmcs_Z_1_p = 1 - pchisq(cmcs_Z_1, 1)
  
  eps_p_1[j-20,3] = cmcs_Z_1_p
  
  ### eps = 2
  
  cmcs_Mp_2 = cmcs_M
  # cmcs_noise_2 = rexp(nrow(cmcs_M)*2+2, eps[1]/S) * sign(runif(nrow(cmcs_M)*2+2, -1, 1))
  cmcs_noise_2 = rexp(nrow(cmcs_M)*4+2, eps[1]/S) * sign(runif(nrow(cmcs_M)*4+2, -1, 1))
  cmcs_noise_2 = pmax(cmcs_noise_2, 0)
  
  cmcs_Mp_2[1,2] = cmcs_Mp_2[1,2] + cmcs_noise_2[1]
  cmcs_Mp_2[1,5] = cmcs_Mp_2[1,5] + cmcs_noise_2[2]
  #cmcs_Mp_2[,3] = cmcs_Mp_2[,3] + cmcs_noise_2[3:(nrow(cmcs_M)+2)]
  #cmcs_Mp_2[,6] = cmcs_Mp_2[,6] + cmcs_noise_2[-c(1:(nrow(cmcs_M)+2))]
  
  cmcs_Mp_2[,3] = cmcs_Mp_2[,3] + cmcs_noise_2[3:(nrow(cmcs_M)+2)]
  cmcs_Mp_2[,6] = cmcs_Mp_2[,6] + cmcs_noise_2[(nrow(cmcs_M)+3):(nrow(cmcs_M)*2+2)]
  cmcs_Mp_2[,4] = cmcs_Mp_2[,4] + cmcs_noise_2[(nrow(cmcs_M)*2+3):(nrow(cmcs_M)*3+2)]
  cmcs_Mp_2[,7] = cmcs_Mp_2[,7] + cmcs_noise_2[(nrow(cmcs_M)*3+3):(nrow(cmcs_M)*4+2)]
  
  
  for(i in 2:nrow(cmcs_Mp_2)) {
    cmcs_Mp_2[i,2] = cmcs_Mp_2[i-1,2] - cmcs_Mp_2[i-1,3] - cmcs_Mp_2[i-1,4]
    cmcs_Mp_2[i,5] = cmcs_Mp_2[i-1,5] - cmcs_Mp_2[i-1,6] - cmcs_Mp_2[i-1,7]
  }
  
  cmcs_Mp_2[,8] = (cmcs_Mp_2[,3] + cmcs_Mp_2[,6]) * cmcs_Mp_2[,2] / (cmcs_Mp_2[,2] + cmcs_Mp_2[,5])
  cmcs_Mp_2[,9] = (cmcs_Mp_2[,3] + cmcs_Mp_2[,6]) * cmcs_Mp_2[,5] / (cmcs_Mp_2[,2] + cmcs_Mp_2[,5])
  
  cmcs_Mp_2_lim = min(max(which(cmcs_Mp_2[,2] > 0)), max(which(cmcs_Mp_2[,5] > 0)))
  
  cmcs_Mp_2 = cmcs_Mp_2[1:cmcs_Mp_2_lim,]
  
  cmcs_Z_2 = sum(cmcs_Mp_2[,3] - cmcs_Mp_2[,8])^2 / sum(cmcs_Mp_2[,10])
  cmcs_Z_2_p = 1 - pchisq(cmcs_Z_2, 1)
  
  eps_p_2[j-20,3] = cmcs_Z_2_p
  
  ### eps = 3
  
  cmcs_Mp_3 = cmcs_M
  # cmcs_noise_3 = rexp(nrow(cmcs_M)*2+2, eps[1]/S) * sign(runif(nrow(cmcs_M)*2+2, -1, 1))
  cmcs_noise_3 = rexp(nrow(cmcs_M)*4+2, eps[1]/S) * sign(runif(nrow(cmcs_M)*4+2, -1, 1))
  cmcs_noise_3 = pmax(cmcs_noise_3, 0)
  
  cmcs_Mp_3[1,2] = cmcs_Mp_3[1,2] + cmcs_noise_3[1]
  cmcs_Mp_3[1,5] = cmcs_Mp_3[1,5] + cmcs_noise_3[2]
  #cmcs_Mp_3[,3] = cmcs_Mp_3[,3] + cmcs_noise_3[3:(nrow(cmcs_M)+2)]
  #cmcs_Mp_3[,6] = cmcs_Mp_3[,6] + cmcs_noise_3[-c(1:(nrow(cmcs_M)+2))]
  
  cmcs_Mp_3[,3] = cmcs_Mp_3[,3] + cmcs_noise_3[3:(nrow(cmcs_M)+2)]
  cmcs_Mp_3[,6] = cmcs_Mp_3[,6] + cmcs_noise_3[(nrow(cmcs_M)+3):(nrow(cmcs_M)*2+2)]
  cmcs_Mp_3[,4] = cmcs_Mp_3[,4] + cmcs_noise_3[(nrow(cmcs_M)*2+3):(nrow(cmcs_M)*3+2)]
  cmcs_Mp_3[,7] = cmcs_Mp_3[,7] + cmcs_noise_3[(nrow(cmcs_M)*3+3):(nrow(cmcs_M)*4+2)]
  
  
  for(i in 2:nrow(cmcs_Mp_3)) {
    cmcs_Mp_3[i,2] = cmcs_Mp_3[i-1,2] - cmcs_Mp_3[i-1,3] - cmcs_Mp_3[i-1,4]
    cmcs_Mp_3[i,5] = cmcs_Mp_3[i-1,5] - cmcs_Mp_3[i-1,6] - cmcs_Mp_3[i-1,7]
  }
  
  cmcs_Mp_3[,8] = (cmcs_Mp_3[,3] + cmcs_Mp_3[,6]) * cmcs_Mp_3[,2] / (cmcs_Mp_3[,2] + cmcs_Mp_3[,5])
  cmcs_Mp_3[,9] = (cmcs_Mp_3[,3] + cmcs_Mp_3[,6]) * cmcs_Mp_3[,5] / (cmcs_Mp_3[,2] + cmcs_Mp_3[,5])
  
  cmcs_Mp_3_lim = min(max(which(cmcs_Mp_3[,2] > 0)), max(which(cmcs_Mp_3[,5] > 0)))
  
  cmcs_Mp_3 = cmcs_Mp_3[1:cmcs_Mp_3_lim,]
  
  cmcs_Z_3 = sum(cmcs_Mp_3[,3] - cmcs_Mp_3[,8])^2 / sum(cmcs_Mp_3[,10])
  cmcs_Z_3_p = 1 - pchisq(cmcs_Z_3, 1)
  
  eps_p_3[j-20,3] = cmcs_Z_3_p
}
