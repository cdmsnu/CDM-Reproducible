library(meta); library(grid);
library(metaBMA)
library(dplyr)
library(ReplicationSuccess)
library(metap)
library(magrittr)

#------ data load ------
data = read.csv("C:/Users/dy/Desktop/snu_bh/Meta analysis/R/1st_thyroid/1st_thyroid_table1.csv", header = TRUE, encoding = 'UTF-8') %>% na.omit
#data$TE = c(log(data[data$Analysis == "snubh",]$HR), log(data[data$Analysis == "snuh",]$HR))
data = data %>% subset(select= c(Time.at.risk, Database, rr, ci95lb, ci95ub))
data['Time.at.risk'] = c(rep('1y5y', 4), rep('1y10y', 4), rep('1yobs', 4))
#data = rbind(data, c("1yobs", "Taiwan", 1.279, 0.92, 1.776))
  

data2 = read.csv("C:/Users/dy/Desktop/snu_bh/Meta analysis/R/1st_thyroid/1st_thyroid_table2.csv", header = TRUE, encoding = 'UTF-8') %>% na.omit
data2 = data2 %>% subset(select= c(Time.at.risk, Database, HR, LL, UL))
colnames(data2)= colnames(data)
data2['Time.at.risk'] = c(rep('1y5y_', 4), rep('5yobs', 3))

data = rbind(data, data2)
data = transform(data, rr = as.numeric(rr), ci95lb = as.numeric(ci95lb), ci95ub = as.numeric(ci95ub))

meta_data=data %>%
  subset(select= c(Time.at.risk, Database, rr, ci95lb, ci95ub)) %>%
  mutate(se= (ci95ub- ci95lb)/(2*qnorm(0.975))) %>%
  mutate(logHR= log(rr)) %>%
  mutate(se_logHR = (log(ci95ub)-log(ci95lb))/(2*qnorm(0.975))) %>%
  mutate(p_onesided = ci2p(ci95lb, ci95ub, ratio= TRUE, alternative = "one.sided"))

#colnames(meta_data)[1]='databseId'
#meta_data['databaseId']= c(rep('1y5y', 4), rep('1y10y', 4), rep('1yobs', 4))
#meta_data_res = meta_data %>% group_by(analysisId) %>% mutate(combined_p = allmetap(p_onesided, method= "all"))

meta_data_1y5y  = meta_data %>% subset(Time.at.risk == '1y5y')
meta_data_1y10y = meta_data %>% subset(Time.at.risk == '1y10y')
meta_data_1yobs = meta_data %>% subset(Time.at.risk == '1yobs')

meta_data_1y5y_   = meta_data %>% subset(Time.at.risk == '1y5y_')
meta_data_5yobs  = meta_data %>% subset(Time.at.risk == '5yobs')

p_1y5y  = allmetap(meta_data_1y5y$p_onesided, method = "all")[c(1,6,8),]
p_1y10y = allmetap(meta_data_1y10y$p_onesided, method = "all")[c(1,6,8),]
p_1yobs = allmetap(meta_data_1yobs$p_onesided, method = "all")[c(1,6,8),]

p_1y5y_  = allmetap(meta_data_1y5y_$p_onesided, method = "all")[c(1,6,8),]
p_5yobs  = allmetap(meta_data_5yobs$p_onesided, method = "all")[c(1,6,8),]

#p_all = allmetap(meta_data$p_onesided, method="all")[c(1,6,8),]
#p_all

p_all = rbind(p_1y5y, p_1y10y, p_1yobs, p_1y5y_, p_5yobs)
p_all['Time.at.risk'] = c(rep('1y5y', 3), rep('1y10y', 3), rep('1yobs', 3), rep('1y5y_', 3), rep('5yobs', 3))
p_all

res_1y5y  = metagen(TE = logHR, seTE = se_logHR, sm = "HR", studlab = Database, hakn = TRUE, method.tau= "PM", data = meta_data_1y5y)
res_1y10y = metagen(TE = logHR, seTE = se_logHR, sm = "HR", studlab = Database, hakn = TRUE, method.tau= "PM", data = meta_data_1y10y)
res_1yobs = metagen(TE = logHR, seTE = se_logHR, sm = "HR", studlab = Database, hakn = TRUE, method.tau= "PM", data = meta_data_1yobs)

res_1y5y_  = metagen(TE = logHR, seTE = se_logHR, sm = "HR", studlab = Database, hakn = TRUE, method.tau= "PM", data = meta_data_1y5y_)
res_5yobs = metagen(TE = logHR, seTE = se_logHR, sm = "HR", studlab = Database, hakn = TRUE, method.tau= "PM", data = meta_data_5yobs)
#res = metagen(TE = logHR, seTE = se_logHR, sm = "HR", studlab = Database, hakn = TRUE, data = meta_data)

png(file="C:/Users/dy/Desktop/snu_bh/Meta analysis/R/1st_thyroid/forest_plot/forest_1y5y.png",width=600, height=350)
forest(res_1y5y)
grid.text("Time at risk : 1y~5y", .5, .8, gp=gpar(cex=1.5))
dev.off()

png(file="C:/Users/dy/Desktop/snu_bh/Meta analysis/R/1st_thyroid/forest_plot/forest_1y10y.png",width=600, height=350)
forest(res_1y10y)
grid.text("Time at risk : 1y~10y", .5, .8, gp=gpar(cex=1.5))
dev.off()

png(file="C:/Users/dy/Desktop/snu_bh/Meta analysis/R/1st_thyroid/forest_plot/forest_1yobs.png",width=600, height=350)
forest(res_1yobs)
grid.text("Time at risk : 1y~obs", .5, .8, gp=gpar(cex=1.5))
dev.off()

png(file="C:/Users/dy/Desktop/snu_bh/Meta analysis/R/1st_thyroid/forest_plot/forest_1y5y_.png",width=600, height=350)
forest(res_1y5y_)
dev.off()

png(file="C:/Users/dy/Desktop/snu_bh/Meta analysis/R/1st_thyroid/forest_plot/forest_5yobs.png",width=600, height=350)
forest.meta(res_5yobs)
grid.text("Time at risk : 5y~obs", .5, .8, gp=gpar(cex=1.5))
dev.off()


bma_1y5y = meta_bma(logHR, se_logHR, Database, data= meta_data_1y5y,
                    d = prior("norm", c(mean=0, sd=0.3)),  ### Half-Normal Distribution
                    tau = prior("t", c(location=0, scale=0.3, nu=1), lower=0), ### Half-Cauchy Distribution
                    control = list(adapt_delta = 0.999), # To avoid "There were 3 divergent transitions after warmup"
                    prior = c(1,1,1,1), silent_stan = FALSE) # in the order c(fixed_H0,fixed_H1,random_H0,random_H1)

bma_1y10y = meta_bma(logHR, se_logHR, Database, data= meta_data_1y10y,
                     d = prior("norm", c(mean=0, sd=0.3)),  ### Half-Normal Distribution
                     tau = prior("t", c(location=0, scale=0.3, nu=1), lower=0), ### Half-Cauchy Distribution
                     control = list(adapt_delta = 0.999), # To avoid "There were 3 divergent transitions after warmup"
                     prior = c(1,1,1,1), silent_stan = FALSE) # in the order c(fixed_H0,fixed_H1,random_H0,random_H1)

bma_1yobs = meta_bma(logHR, se_logHR, Database, data= meta_data_1yobs,
                    d = prior("norm", c(mean=0, sd=0.3)),  ### Half-Normal Distribution
                    tau = prior("t", c(location=0, scale=0.3, nu=1), lower=0), ### Half-Cauchy Distribution
                    control = list(adapt_delta = 0.999), # To avoid "There were 3 divergent transitions after warmup"
                    prior = c(1,1,1,1), silent_stan = FALSE) # in the order c(fixed_H0,fixed_H1,random_H0,random_H1)


bma_1y5y_ = meta_bma(logHR, se_logHR, Database, data= meta_data_1y5y_,
                    d = prior("norm", c(mean=0, sd=0.3)),  ### Half-Normal Distribution
                    tau = prior("t", c(location=0, scale=0.3, nu=1), lower=0), ### Half-Cauchy Distribution
                    control = list(adapt_delta = 0.999), # To avoid "There were 3 divergent transitions after warmup"
                    prior = c(1,1,1,1), silent_stan = FALSE) # in the order c(fixed_H0,fixed_H1,random_H0,random_H1)

bma_5yobs = meta_bma(logHR, se_logHR, Database, data= meta_data_5yobs,
                     d = prior("norm", c(mean=0, sd=0.3)),  ### Half-Normal Distribution
                     tau = prior("t", c(location=0, scale=0.3, nu=1), lower=0), ### Half-Cauchy Distribution
                     control = list(adapt_delta = 0.99999), # To avoid "There were 3 divergent transitions after warmup"
                     prior = c(1,1,1,1), silent_stan = FALSE) # in the order c(fixed_H0,fixed_H1,random_H0,random_H1)



# bma = meta_bma(logHR, se_logHR, Database, data= meta_data,
#                d = prior("norm", c(mean=0, sd=0.3)),  ### Half-Normal Distribution
#                tau = prior("t", c(location=0, scale=0.3, nu=1), lower=0), ### Half-Cauchy Distribution
#                control = list(adapt_delta = 0.999), # To avoid "There were 3 divergent transitions after warmup"
#                prior = c(1,1,1,1), silent_stan = FALSE) # in the order c(fixed_H0,fixed_H1,random_H0,random_H1)


############# 1y5y #############
bma_df= data.frame("time.at.risk"=c("1y5y", "1y10y", "1yobs", "1y5y_", "5yobs"), "Prob"=rep(0,5))
cdf=cumsum(bma_1y5y$posterior_d(seq(-2, 2, 0.01))*0.01)
exp(seq(-2, 2, 0.01)[which(0.025 <= cdf & cdf <= 0.975)])
plot_posterior(bma_1y5y, parameter = "d", from = -0.4, to= 0.8, main= paste("1y5y,P(HR > 1 | X)=", round(prob, 4)))
prob=sum(bma_1y5y$posterior_d(seq(0, 2, 0.01))*0.01)
bma_df[1,2]= prob

############# 1y10y #############
cdf=cumsum(bma_1y10y$posterior_d(seq(-2, 2, 0.01))*0.01)
exp(seq(-2, 2, 0.01)[which(0.025 <= cdf & cdf <= 0.975)])
plot_posterior(bma_1y10y, parameter = "d", from = -0.4, to= 0.6, main=paste("1y10y, P(HR > 1 | X)=", round(prob, 4)))
prob=sum(bma_1y10y$posterior_d(seq(0, 2, 0.01))*0.01)
bma_df[2,2]= prob


############# 1yobs #############

cdf=cumsum(bma_1y10y$posterior_d(seq(-2, 2, 0.01))*0.01)
exp(seq(-2, 2, 0.01)[which(0.025 <= cdf & cdf <= 0.975)])
plot_posterior(bma_1yobs, parameter = "d", from = -0.5, to= 1.4, main=paste("1yobs, P(HR > 1 | X)=", round(prob, 4)))
prob=sum(bma_1yobs$posterior_d(seq(0, 2, 0.01))*0.01)
bma_df[3,2]= prob

############# 1y5y_ #############

cdf=cumsum(bma_1y10y$posterior_d(seq(-2, 2, 0.01))*0.01)
exp(seq(-2, 2, 0.01)[which(0.025 <= cdf & cdf <= 0.975)])
plot_posterior(bma_1y5y_, parameter = "d", from = -0.4, to= 0.8, main= paste("1y5y_,P(HR > 1 | X)=", round(prob, 4)))
prob=sum(bma_1y5y_$posterior_d(seq(0, 2, 0.01))*0.01)
bma_df[4,2]= prob

############# 5yobs #############

cdf=cumsum(bma_5yobs$posterior_d(seq(-2, 2, 0.01))*0.01)
exp(seq(-2, 2, 0.01)[which(0.025 <= cdf & cdf <= 0.975)])
prob=sum(bma_5yobs$posterior_d(seq(0, 2, 0.01))*0.01)
plot_posterior(bma_5yobs, parameter = "d", from = -0.4, to= 0.6, main=paste("5yobs, P(HR > 1 | X)=", round(prob, 4)))
bma_df[5,2]= prob


bma_df[,2]=round(bma_df[,2],4)
bma_df
#prob=sum(bma$posterior_d(seq(0, 2, 0.01))*0.01)
#plot_posterior(bma, parameter = "d", from = -0.4, to= 0.8, main= paste("P(HR > 1 | X)=", round(prob, 4)))

