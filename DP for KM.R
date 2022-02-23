

####################   Differential Privacy for Kaplan-Meier curve    ####################

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
  if(j == 1) {
    snubh_comp_S[j] = (snubh_comp_M[1,1] - snubh_comp_M[1,2]) / snubh_comp_M[1,1]
  } else {
    snubh_comp_S[j] = snubh_comp_S[j-1] * (snubh_comp_M[j,1] - snubh_comp_M[j,2]) / snubh_comp_M[j,1]
  }
}





