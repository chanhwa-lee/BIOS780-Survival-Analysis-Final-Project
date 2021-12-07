library(survival)
library(dplyr)
# install.packages("survminer")
library(survminer)
library(cowplot)

setwd("C:/Users/Chanhwa/OneDrive - University of North Carolina at Chapel Hill/2021/FA21/BIOS780/Final Project")

#--------------------#
#                    #
#    Introduction    #
#                    #
#--------------------#


# Load data
data <- read.table("colon.txt", header = F)
colnames(data) <- c("patno", 
                    
                    "treat", 
                    # 1: Obs, 2: Lev, 3: 5-FU + Lev
                    
                    "dayreg", "mthreg", "yrreg", 
                    # Date of registration
                    
                    "daylc", "mthlc", "yrlc",  
                    # Date of last contact or death
                    
                    "status",
                    # 0: Alive (censored), 1: Dead
                    
                    "dayprgrel", "mthprgrel", "yrprgrel",
                    # Date of progress/relapse
                    
                    "stratno",
                    # stratification level
                    
                    "sex",
                    # 0: Female, 1: Male
                    
                    "ptlsite", "histype", "cdifintg", "extloc",
                    # categorical variables
                    
                    "obstrc", "perfor", "adhere", "regimp", 
                    # binary variables
                    
                    "ceapreop",
                    # pre-operative CEA level
                    
                    "noposnd", 
                    # Nodal involvement (number of positive nodes)
                    
                    "dayres", "mthres", "yrres",
                    # Date of tumor resection
                    
                    "dayrx", "mthrx", "yrrx",
                    # Date of start of treatment
                    
                    "prgrel",
                    # 0: Progress or relapse (worse) 1: no (good)
                    
                    "age", 
                    
                    "days",
                    # Days from tumor resection to start of treatment
                    
                    "time", 
                    # Days from registration to death or last contact
                    
                    "group"
                    # 1: Data from SWOG 8591, 2: Data from NCCTG
                    )


# Compute time to event
data$dayprgrel <- as.numeric(data$dayprgrel)
data$mthprgrel <- as.numeric(data$mthprgrel)
data$yrprgrel  <- as.numeric(data$yrprgrel)

data$dayres <- as.numeric(data$dayres)
data$mthres <- as.numeric(data$mthres)
data$yrres  <- as.numeric(data$yrres)

data$dayrx <- as.numeric(data$dayrx)
data$mthrx <- as.numeric(data$mthrx)
data$yrrx  <- as.numeric(data$yrrx)

data$days <- as.numeric(data$days)

data = data %>% 
  mutate(datereg = as.Date(paste(mthreg, dayreg, yrreg-1900, sep = "/"), format = "%m/%d/%y"),
         datelc  = as.Date(paste(mthlc,  daylc,  yrlc -1900, sep = "/"), format = "%m/%d/%y"),
         dateprgrel = as.Date(paste(mthprgrel, dayprgrel, yrprgrel-1900, sep = "/"), format = "%m/%d/%y"),
         dateres = as.Date(paste(mthres, dayres, yrres-1900, sep = "/"), format = "%m/%d/%y"),
         daterx  = as.Date(paste(mthrx,  dayrx,  yrrx -1900, sep = "/"), format = "%m/%d/%y")
         ) %>%
  select(-c(mthreg, dayreg, yrreg, mthlc,  daylc,  yrlc, mthprgrel, dayprgrel, yrprgrel, mthres, dayres, yrres,mthrx,  dayrx,  yrrx))

# Check whether manual computation is the same with original data
data = data %>% mutate(days_manual = as.numeric(daterx - dateres), 
                       time_manual = as.numeric(datelc - datereg))
table(data$days - data$days_manual)
table(data$time - data$time_manual)

# Check the delay between dates
table(data$datereg - data$dateres)
table(data$daterx  - data$datereg)
table(data$dateprgrel  - data$daterx)
table(data$datelc - data$dateprgrel)

# Since datereg and datelc are the only zero-missing data, we use this dates
# Time to death = Days from reg to last contact = data$time
# Time to disease recurrence = Days from reg to prgrel (if prgrel != NA) / last contact (if prgrel == NA)
data = data %>% mutate(
  time_dea = time,
  time_rec = ifelse(prgrel == 1, dateprgrel, datelc) - as.numeric(datereg),
  status_dea = status,
  status_rec = prgrel
)

# Finally organize columns
data$noposnd <- as.numeric(data$noposnd)
data$ceapreop <- as.numeric(data$ceapreop)
data$treat <- as.character(data$treat)
data$sex <- as.character(data$sex)
data$obstrc <- as.character(data$obstrc)
data = data %>% select(-c(patno, stratno, days, time, datereg, datelc, dateprgrel, 
                          dateres, daterx, days_manual, time_manual, status, prgrel))

# How many records are missing?
data.1 = data %>% filter(group == 1) %>% select(-group)
data.2 = data %>% filter(group == 2) %>% select(-group)

apply(X = data, MARGIN = 2, function(x) sum(x %in% c(".",NA)))
apply(X = data.1, MARGIN = 2, function(x) sum(x %in% c(".",NA)))
apply(X = data.2, MARGIN = 2, function(x) sum(x %in% c(".",NA)))

# ptlsite, histype, cdifintg, extloc, regimp are missing for the same 243 subjects
# perfor, adhere are missing for 7 subjects
# noposnd are missing for 31 subjects
# Group2 subjects are missing ptlsite, histype, cdifintg, extloc, regimp => Set this as validation set?


#---------------------#
#                     #
#  TRT on Recurrence  #
#                     #
#---------------------#

# KM plots #
fit.rec <- survfit(Surv(time_rec, as.numeric(status_rec)) ~ treat, data = data.1)
plot(fit.rec, col = 1:3, xlab = "Days", ylab = "Survival Probability")
legend(0, .3, c("Observation", "Levamisole", "5-FU + Levamisole"), lty = 1, col = 1:3) 
# title("Kaplan-Meier Curves of Disease recurrence by Treatment \nfor Colon Cancer study") 

# Cox PH PMLE
fit <- coxph(Surv(time_rec, status_rec) ~ treat, data = data.1, method = 'breslow')
summary(fit)

# Log rank tests #
Lev.Obs.fit.rec <- survdiff(Surv(time_rec, status_rec) ~ treat, data = data.1 %>% filter(treat !=3))
Lev.Obs.fit.rec

FLe.Obs.fit.rec <- survdiff(Surv(time_rec, status_rec) ~ treat, data = data.1 %>% filter(treat !=2))
FLe.Obs.fit.rec

Lev.FLe.fit.rec <- survdiff(Surv(time_rec, status_rec) ~ treat, data = data.1 %>% filter(treat !=1))
Lev.FLe.fit.rec



#---------------------#
#                     #
#  TRT on Surv time   #
#                     #
#---------------------#

# KM plots #
fit.dea <- survfit(Surv(time_dea, status_dea) ~ treat, data = data.1)
plot(fit.dea, col = 1:3, xlab = "Days", ylab = "Survival Probability")
legend(0, .3, c("Observation", "Levamisole", "5-FU + Levamisole"), lty = 1, col = 1:3) 
# title("Kaplan-Meier Curves of Survival time by Treatment \nfor Colon Cancer study") 

# Cox PH PMLE
fit <- coxph(Surv(time_dea, status_dea) ~ treat, data = data.1, method = 'breslow')
summary(fit)

# Log rank tests #
Lev.Obs.fit.dea <- survdiff(Surv(time_dea, status_dea) ~ treat, data = data.1 %>% filter(treat !=3))
Lev.Obs.fit.dea

FLe.Obs.fit.dea <- survdiff(Surv(time_dea, status_dea) ~ treat, data = data.1 %>% filter(treat !=2))
FLe.Obs.fit.dea

Lev.FLe.fit.dea <- survdiff(Surv(time_dea, status_dea) ~ treat, data = data.1 %>% filter(treat !=1))
Lev.FLe.fit.dea


#----------------------------------------#
#                                        #
#   Natural history model on Surv time   #
#                                        #
#----------------------------------------#

#----------------------------------------#
#                                        #
#        Univariate statistics           #
#                                        #
#----------------------------------------#

for(covar in colnames(data.1)){
  if(!(covar %in% c("time_dea", "time_rec", "status_dea", "status_rec"))){
    print(covar)
    
    # Summary table
    if(class(data.1[[covar]]) == "character"){
      print(table(data.1[[covar]], useNA = "always"))
    }
    else{
      print(summary(data.1[[covar]]))
    }
    
    # Univariate coxph model statistics
    fit <- coxph(Surv(data.1$time_dea, data.1$status_dea) ~ data.1[[covar]], method = 'breslow')
    print(summary(fit))
  }
}


#----------------------------------------#
#                                        #
#     Step-down variable selection       #
#                                        #
#----------------------------------------#

# -- variable "age" should be included even if it is statistically insignificant predictor
# -- variable "trt" should be included

# Step1
data.vs <- data.1 %>% select(-c(ptlsite, histype, perfor, ceapreop, time_rec, status_rec))
fit.1 <- coxph(Surv(time_dea, status_dea) ~ ., data = data.vs, method = 'breslow')
fit.1.df = length(fit.1$coefficients)

result.1 <- data.frame(variable=colnames(data.vs %>% select(-c(time_dea, status_dea))), stat = 0, df = 0, pval = 0)

  # Exclude sex
  fit.1.sex <- coxph(Surv(time_dea, status_dea) ~ .-sex, data = data.vs, method = 'breslow')
  fit.1.sex.df = length(fit.1.sex$coefficients)
  stat = 2*(fit.1$loglik[2] - fit.1.sex$loglik[2])
  df = fit.1.df - fit.1.sex.df
  pval = 1 - pchisq(q = stat, df = df)
  result.1$stat[result.1$variable == "sex"] = stat
  result.1$df  [result.1$variable == "sex"] = df
  result.1$pval[result.1$variable == "sex"] = pval

  # Exclude cdifintg
  fit.1.cdifintg <- coxph(Surv(time_dea, status_dea) ~ .-cdifintg, data = data.vs, method = 'breslow')
  fit.1.cdifintg.df = length(fit.1.cdifintg$coefficients)
  stat = 2*(fit.1$loglik[2] - fit.1.cdifintg$loglik[2])
  df = fit.1.df - fit.1.cdifintg.df
  pval = 1 - pchisq(q = stat, df = df)
  result.1$stat[result.1$variable == "cdifintg"] = stat
  result.1$df  [result.1$variable == "cdifintg"] = df
  result.1$pval[result.1$variable == "cdifintg"] = pval
  
  # Exclude extloc
  fit.1.extloc <- coxph(Surv(time_dea, status_dea) ~ .-extloc, data = data.vs, method = 'breslow')
  fit.1.extloc.df = length(fit.1.extloc$coefficients)
  stat = 2*(fit.1$loglik[2] - fit.1.extloc$loglik[2])
  df = fit.1.df - fit.1.extloc.df
  pval = 1 - pchisq(q = stat, df = df)
  result.1$stat[result.1$variable == "extloc"] = stat
  result.1$df  [result.1$variable == "extloc"] = df
  result.1$pval[result.1$variable == "extloc"] = pval
  
  # Exclude obstrc
  fit.1.obstrc <- coxph(Surv(time_dea, status_dea) ~ .-obstrc, data = data.vs, method = 'breslow')
  fit.1.obstrc.df = length(fit.1.obstrc$coefficients)
  stat = 2*(fit.1$loglik[2] - fit.1.obstrc$loglik[2])
  df = fit.1.df - fit.1.obstrc.df
  pval = 1 - pchisq(q = stat, df = df)
  result.1$stat[result.1$variable == "obstrc"] = stat
  result.1$df  [result.1$variable == "obstrc"] = df
  result.1$pval[result.1$variable == "obstrc"] = pval
  
  # Exclude adhere
  fit.1.adhere <- coxph(Surv(time_dea, status_dea) ~ .-adhere, data = data.vs, method = 'breslow')
  fit.1.adhere.df = length(fit.1.adhere$coefficients)
  stat = 2*(fit.1$loglik[2] - fit.1.adhere$loglik[2])
  df = fit.1.df - fit.1.adhere.df
  pval = 1 - pchisq(q = stat, df = df)
  result.1$stat[result.1$variable == "adhere"] = stat
  result.1$df  [result.1$variable == "adhere"] = df
  result.1$pval[result.1$variable == "adhere"] = pval

  # Exclude regimp
  fit.1.regimp <- coxph(Surv(time_dea, status_dea) ~ .-regimp, data = data.vs, method = 'breslow')
  fit.1.regimp.df = length(fit.1.regimp$coefficients)
  stat = 2*(fit.1$loglik[2] - fit.1.regimp$loglik[2])
  df = fit.1.df - fit.1.regimp.df
  pval = 1 - pchisq(q = stat, df = df)
  result.1$stat[result.1$variable == "regimp"] = stat
  result.1$df  [result.1$variable == "regimp"] = df
  result.1$pval[result.1$variable == "regimp"] = pval
  
  # Exclude noposnd
  fit.1.noposnd <- coxph(Surv(time_dea, status_dea) ~ .-noposnd, data = data.vs, method = 'breslow')
  fit.1.noposnd.df = length(fit.1.noposnd$coefficients)
  stat = 2*(fit.1$loglik[2] - fit.1.noposnd$loglik[2])
  df = fit.1.df - fit.1.noposnd.df
  pval = 1 - pchisq(q = stat, df = df)
  result.1$stat[result.1$variable == "noposnd"] = stat
  result.1$df  [result.1$variable == "noposnd"] = df
  result.1$pval[result.1$variable == "noposnd"] = pval
  

  result.1[which.max(result.1$pval),]
  # sex is excluded
  

# Step2
  data.vs <- data.vs %>% select(-sex)
  fit.2 <- coxph(Surv(time_dea, status_dea) ~ ., data = data.vs, method = 'breslow')
  fit.2.df = length(fit.2$coefficients)
  
  result.2 <- data.frame(variable=colnames(data.vs %>% select(-c(time_dea, status_dea))), stat = 0, df = 0, pval = 0)
  
  # Exclude cdifintg
  fit.2.cdifintg <- coxph(Surv(time_dea, status_dea) ~ .-cdifintg, data = data.vs, method = 'breslow')
  fit.2.cdifintg.df = length(fit.2.cdifintg$coefficients)
  stat = 2*(fit.2$loglik[2] - fit.2.cdifintg$loglik[2])
  df = fit.2.df - fit.2.cdifintg.df
  pval = 1 - pchisq(q = stat, df = df)
  result.2$stat[result.2$variable == "cdifintg"] = stat
  result.2$df  [result.2$variable == "cdifintg"] = df
  result.2$pval[result.2$variable == "cdifintg"] = pval
  
  # Exclude extloc
  fit.2.extloc <- coxph(Surv(time_dea, status_dea) ~ .-extloc, data = data.vs, method = 'breslow')
  fit.2.extloc.df = length(fit.2.extloc$coefficients)
  stat = 2*(fit.2$loglik[2] - fit.2.extloc$loglik[2])
  df = fit.2.df - fit.2.extloc.df
  pval = 1 - pchisq(q = stat, df = df)
  result.2$stat[result.2$variable == "extloc"] = stat
  result.2$df  [result.2$variable == "extloc"] = df
  result.2$pval[result.2$variable == "extloc"] = pval
  
  # Exclude obstrc
  fit.2.obstrc <- coxph(Surv(time_dea, status_dea) ~ .-obstrc, data = data.vs, method = 'breslow')
  fit.2.obstrc.df = length(fit.2.obstrc$coefficients)
  stat = 2*(fit.2$loglik[2] - fit.2.obstrc$loglik[2])
  df = fit.2.df - fit.2.obstrc.df
  pval = 1 - pchisq(q = stat, df = df)
  result.2$stat[result.2$variable == "obstrc"] = stat
  result.2$df  [result.2$variable == "obstrc"] = df
  result.2$pval[result.2$variable == "obstrc"] = pval
  
  # Exclude adhere
  fit.2.adhere <- coxph(Surv(time_dea, status_dea) ~ .-adhere, data = data.vs, method = 'breslow')
  fit.2.adhere.df = length(fit.2.adhere$coefficients)
  stat = 2*(fit.2$loglik[2] - fit.2.adhere$loglik[2])
  df = fit.2.df - fit.2.adhere.df
  pval = 1 - pchisq(q = stat, df = df)
  result.2$stat[result.2$variable == "adhere"] = stat
  result.2$df  [result.2$variable == "adhere"] = df
  result.2$pval[result.2$variable == "adhere"] = pval
  
  # Exclude regimp
  fit.2.regimp <- coxph(Surv(time_dea, status_dea) ~ .-regimp, data = data.vs, method = 'breslow')
  fit.2.regimp.df = length(fit.2.regimp$coefficients)
  stat = 2*(fit.2$loglik[2] - fit.2.regimp$loglik[2])
  df = fit.2.df - fit.2.regimp.df
  pval = 1 - pchisq(q = stat, df = df)
  result.2$stat[result.2$variable == "regimp"] = stat
  result.2$df  [result.2$variable == "regimp"] = df
  result.2$pval[result.2$variable == "regimp"] = pval
  
  # Exclude noposnd
  fit.2.noposnd <- coxph(Surv(time_dea, status_dea) ~ .-noposnd, data = data.vs, method = 'breslow')
  fit.2.noposnd.df = length(fit.2.noposnd$coefficients)
  stat = 2*(fit.2$loglik[2] - fit.2.noposnd$loglik[2])
  df = fit.2.df - fit.2.noposnd.df
  pval = 1 - pchisq(q = stat, df = df)
  result.2$stat[result.2$variable == "noposnd"] = stat
  result.2$df  [result.2$variable == "noposnd"] = df
  result.2$pval[result.2$variable == "noposnd"] = pval
  
  # Final Decision
  result.2[which.max(result.2$pval),]
  # cdifintg is excluded
  
  
# Step3
  data.vs <- data.vs %>% select(-cdifintg)
  fit.3 <- coxph(Surv(time_dea, status_dea) ~ ., data = data.vs, method = 'breslow')
  fit.3.df = length(fit.3$coefficients)
  
  result.3 <- data.frame(variable=colnames(data.vs %>% select(-c(time_dea, status_dea))), stat = 0, df = 0, pval = 0)
  
  # Exclude extloc
  fit.3.extloc <- coxph(Surv(time_dea, status_dea) ~ .-extloc, data = data.vs, method = 'breslow')
  fit.3.extloc.df = length(fit.3.extloc$coefficients)
  stat = 2*(fit.3$loglik[2] - fit.3.extloc$loglik[2])
  df = fit.3.df - fit.3.extloc.df
  pval = 1 - pchisq(q = stat, df = df)
  result.3$stat[result.3$variable == "extloc"] = stat
  result.3$df  [result.3$variable == "extloc"] = df
  result.3$pval[result.3$variable == "extloc"] = pval
  
  # Exclude obstrc
  fit.3.obstrc <- coxph(Surv(time_dea, status_dea) ~ .-obstrc, data = data.vs, method = 'breslow')
  fit.3.obstrc.df = length(fit.3.obstrc$coefficients)
  stat = 2*(fit.3$loglik[2] - fit.3.obstrc$loglik[2])
  df = fit.3.df - fit.3.obstrc.df
  pval = 1 - pchisq(q = stat, df = df)
  result.3$stat[result.3$variable == "obstrc"] = stat
  result.3$df  [result.3$variable == "obstrc"] = df
  result.3$pval[result.3$variable == "obstrc"] = pval
  
  # Exclude adhere
  fit.3.adhere <- coxph(Surv(time_dea, status_dea) ~ .-adhere, data = data.vs, method = 'breslow')
  fit.3.adhere.df = length(fit.3.adhere$coefficients)
  stat = 2*(fit.3$loglik[2] - fit.3.adhere$loglik[2])
  df = fit.3.df - fit.3.adhere.df
  pval = 1 - pchisq(q = stat, df = df)
  result.3$stat[result.3$variable == "adhere"] = stat
  result.3$df  [result.3$variable == "adhere"] = df
  result.3$pval[result.3$variable == "adhere"] = pval
  
  # Exclude regimp
  fit.3.regimp <- coxph(Surv(time_dea, status_dea) ~ .-regimp, data = data.vs, method = 'breslow')
  fit.3.regimp.df = length(fit.3.regimp$coefficients)
  stat = 2*(fit.3$loglik[2] - fit.3.regimp$loglik[2])
  df = fit.3.df - fit.3.regimp.df
  pval = 1 - pchisq(q = stat, df = df)
  result.3$stat[result.3$variable == "regimp"] = stat
  result.3$df  [result.3$variable == "regimp"] = df
  result.3$pval[result.3$variable == "regimp"] = pval
  
  # Exclude noposnd
  fit.3.noposnd <- coxph(Surv(time_dea, status_dea) ~ .-noposnd, data = data.vs, method = 'breslow')
  fit.3.noposnd.df = length(fit.3.noposnd$coefficients)
  stat = 2*(fit.3$loglik[2] - fit.3.noposnd$loglik[2])
  df = fit.3.df - fit.3.noposnd.df
  pval = 1 - pchisq(q = stat, df = df)
  result.3$stat[result.3$variable == "noposnd"] = stat
  result.3$df  [result.3$variable == "noposnd"] = df
  result.3$pval[result.3$variable == "noposnd"] = pval
  
  # Final Decision
  result.3[which.max(result.3$pval),]
  # adhere is excluded
  fit.3

# Step4
  data.vs <- data.vs %>% select(-adhere)
  fit.4 <- coxph(Surv(time_dea, status_dea) ~ ., data = data.vs, method = 'breslow')
  fit.4.df = length(fit.4$coefficients)
  
  result.4 <- data.frame(variable=colnames(data.vs %>% select(-c(time_dea, status_dea))), stat = 0, df = 0, pval = 0)
  
  # Exclude extloc
  fit.4.extloc <- coxph(Surv(time_dea, status_dea) ~ .-extloc, data = data.vs, method = 'breslow')
  fit.4.extloc.df = length(fit.4.extloc$coefficients)
  stat = 2*(fit.4$loglik[2] - fit.4.extloc$loglik[2])
  df = fit.4.df - fit.4.extloc.df
  pval = 1 - pchisq(q = stat, df = df)
  result.4$stat[result.4$variable == "extloc"] = stat
  result.4$df  [result.4$variable == "extloc"] = df
  result.4$pval[result.4$variable == "extloc"] = pval
  
  # Exclude obstrc
  fit.4.obstrc <- coxph(Surv(time_dea, status_dea) ~ .-obstrc, data = data.vs, method = 'breslow')
  fit.4.obstrc.df = length(fit.4.obstrc$coefficients)
  stat = 2*(fit.4$loglik[2] - fit.4.obstrc$loglik[2])
  df = fit.4.df - fit.4.obstrc.df
  pval = 1 - pchisq(q = stat, df = df)
  result.4$stat[result.4$variable == "obstrc"] = stat
  result.4$df  [result.4$variable == "obstrc"] = df
  result.4$pval[result.4$variable == "obstrc"] = pval
  
  # Exclude regimp
  fit.4.regimp <- coxph(Surv(time_dea, status_dea) ~ .-regimp, data = data.vs, method = 'breslow')
  fit.4.regimp.df = length(fit.4.regimp$coefficients)
  stat = 2*(fit.4$loglik[2] - fit.4.regimp$loglik[2])
  df = fit.4.df - fit.4.regimp.df
  pval = 1 - pchisq(q = stat, df = df)
  result.4$stat[result.4$variable == "regimp"] = stat
  result.4$df  [result.4$variable == "regimp"] = df
  result.4$pval[result.4$variable == "regimp"] = pval
  
  # Exclude noposnd
  fit.4.noposnd <- coxph(Surv(time_dea, status_dea) ~ .-noposnd, data = data.vs, method = 'breslow')
  fit.4.noposnd.df = length(fit.4.noposnd$coefficients)
  stat = 2*(fit.4$loglik[2] - fit.4.noposnd$loglik[2])
  df = fit.4.df - fit.4.noposnd.df
  pval = 1 - pchisq(q = stat, df = df)
  result.4$stat[result.4$variable == "noposnd"] = stat
  result.4$df  [result.4$variable == "noposnd"] = df
  result.4$pval[result.4$variable == "noposnd"] = pval
  
  # Final Decision
  result.4[which.max(result.4$pval),]
  # obstrc is excluded (p-val = 0.064295)
  fit.4

# Step5
  data.vs <- data.vs %>% select(-obstrc)
  fit.5 <- coxph(Surv(time_dea, status_dea) ~ ., data = data.vs, method = 'breslow')
  fit.5.df = length(fit.5$coefficients)
  
  result.5 <- data.frame(variable=colnames(data.vs %>% select(-c(time_dea, status_dea))), stat = 0, df = 0, pval = 0)
  
  # Exclude extloc
  fit.5.extloc <- coxph(Surv(time_dea, status_dea) ~ .-extloc, data = data.vs, method = 'breslow')
  fit.5.extloc.df = length(fit.5.extloc$coefficients)
  stat = 2*(fit.5$loglik[2] - fit.5.extloc$loglik[2])
  df = fit.5.df - fit.5.extloc.df
  pval = 1 - pchisq(q = stat, df = df)
  result.5$stat[result.5$variable == "extloc"] = stat
  result.5$df  [result.5$variable == "extloc"] = df
  result.5$pval[result.5$variable == "extloc"] = pval
  
  # Exclude regimp
  fit.5.regimp <- coxph(Surv(time_dea, status_dea) ~ .-regimp, data = data.vs, method = 'breslow')
  fit.5.regimp.df = length(fit.5.regimp$coefficients)
  stat = 2*(fit.5$loglik[2] - fit.5.regimp$loglik[2])
  df = fit.5.df - fit.5.regimp.df
  pval = 1 - pchisq(q = stat, df = df)
  result.5$stat[result.5$variable == "regimp"] = stat
  result.5$df  [result.5$variable == "regimp"] = df
  result.5$pval[result.5$variable == "regimp"] = pval
  
  # Exclude adhere
  fit.5.adhere <- coxph(Surv(time_dea, status_dea) ~ .-adhere, data = data.vs, method = 'breslow')
  fit.5.adhere.df = length(fit.5.adhere$coefficients)
  stat = 2*(fit.5$loglik[2] - fit.5.adhere$loglik[2])
  df = fit.5.df - fit.5.adhere.df
  pval = 1 - pchisq(q = stat, df = df)
  result.5$stat[result.5$variable == "adhere"] = stat
  result.5$df  [result.5$variable == "adhere"] = df
  result.5$pval[result.5$variable == "adhere"] = pval
  
  # Exclude noposnd
  fit.5.noposnd <- coxph(Surv(time_dea, status_dea) ~ .-noposnd, data = data.vs, method = 'breslow')
  fit.5.noposnd.df = length(fit.5.noposnd$coefficients)
  stat = 2*(fit.5$loglik[2] - fit.5.noposnd$loglik[2])
  df = fit.5.df - fit.5.noposnd.df
  pval = 1 - pchisq(q = stat, df = df)
  result.5$stat[result.5$variable == "noposnd"] = stat
  result.5$df  [result.5$variable == "noposnd"] = df
  result.5$pval[result.5$variable == "noposnd"] = pval
  
  # Final Decision
  result.5[which.max(result.5$pval),]
  # NO variable is excluded because all pvalues < 0.05
  result.5
  fit.5

# Final model: treat + extloc + regimp + noposnd + age
fit.5 <- coxph(Surv(time_dea, status_dea) ~ ., data = data.vs, method = 'breslow')
result.5
fit.5
fit.5.df
fit.5$loglik

# Original model: all variables  
result.1
fit.1
fit.1.df
fit.1$loglik

# Difference between original model and reduced model
2*(fit.1$loglik[2] - fit.5$loglik[2])
1-pchisq(2*(fit.1$loglik[2] - fit.5$loglik[2]), df = fit.1.df - fit.5.df)



#----------------------------------------#
#                                        #
#     Functional form of covariates      #
#                                        #
#----------------------------------------#

# treat + extloc + regimp + noposnd + log(noposnd) + age + log(age)
table(data.vs$noposnd) # ignore data with noposnd == 0 for log transformation
data.vs = data.vs %>% mutate(log_noposnd = log(pmax(noposnd,1)))
fit.5.1 <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + noposnd + log_noposnd+ age + log(age), data = data.vs, method = 'breslow')
fit.5.1
fit.5

fit.5.1$loglik
fit.5$loglik
-2*(fit.5$loglik[2] - fit.5.1$loglik[2])

# treat + extloc + regimp + log(noposnd) + age
fit.5.2 <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age, data = data.vs, method = 'breslow')
fit.5.2
fit.5.2$loglik[2]

# treat + extloc + regimp + log(noposnd) + age + age^2
fit.5.3 <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + I(age**2), data = data.vs, method = 'breslow')
fit.5.3
fit.5.3$loglik[2]

# All 2 way interactions treat * extloc * regimp * log(noposnd) * age
fit.5.4 <- coxph(Surv(time_dea, status_dea) ~ treat:extloc + treat:regimp + treat:log_noposnd + treat:age
                 + extloc:regimp + extloc:log_noposnd + extloc:age
                 + regimp:log_noposnd + regimp:age
                 + log_noposnd:age, data = data.vs, method = 'breslow')
fit.5.4
fit.5.4$loglik[2]
# No significant 2 way interactions if considering all interactions simultaneously

# treat + extloc + regimp + log(noposnd) + age + treat:extloc
fit.5.5 <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + treat:extloc, data = data.vs, method = 'breslow')
fit.5.5
fit.5.5$loglik[2]

# treat + extloc + regimp + log(noposnd) + age + treat:regimp
fit.5.6 <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + treat:regimp, data = data.vs, method = 'breslow')
fit.5.6
fit.5.6$loglik[2]

# treat + extloc + regimp + log(noposnd) + age + treat:log(noposnd)
fit.5.7 <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + treat:log_noposnd, data = data.vs, method = 'breslow')
fit.5.7
fit.5.7$loglik[2]

# treat + extloc + regimp + log(noposnd) + age + treat:age
fit.5.8 <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + treat:age, data = data.vs, method = 'breslow')
fit.5.8
fit.5.8$loglik[2]
# Treat3(5-FU + Lev) and Age interaction pvalue = 0.089 (but then Treat3 main effect insignificant)

# treat + extloc + regimp + log(noposnd) + age + extloc:regimp
fit.5.9 <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + extloc:regimp, data = data.vs, method = 'breslow')
fit.5.9
fit.5.9$loglik[2]

# treat + extloc + regimp + log(noposnd) + age + extloc:log(noposnd)
fit.5.10 <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + extloc:log_noposnd, data = data.vs, method = 'breslow')
fit.5.10
fit.5.10$loglik[2]

# treat + extloc + regimp + log(noposnd) + age + extloc:age
fit.5.11 <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + extloc:age, data = data.vs, method = 'breslow')
fit.5.11
fit.5.11$loglik[2]
# extloc:age is higly significant

# treat + extloc + regimp + log(noposnd) + age + regimp:log(noposnd)
fit.5.12 <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + regimp:log_noposnd, data = data.vs, method = 'breslow')
fit.5.12
fit.5.12$loglik[2]

# treat + extloc + regimp + log(noposnd) + age + regimp:age
fit.5.13 <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + regimp:age, data = data.vs, method = 'breslow')
fit.5.13
fit.5.13$loglik[2]

# treat + extloc + regimp + log(noposnd) + age + log(noposnd):age
fit.5.14 <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + log_noposnd:age, data = data.vs, method = 'breslow')
fit.5.14
fit.5.14$loglik[2]

# treat + extloc + regimp + log(noposnd) + age + extloc:age + log(noposnd):age
fit.5.15 <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + extloc:age + log_noposnd:age, data = data.vs, method = 'breslow')
fit.5.15
fit.5.15$loglik[2]


# Final model : treat + extloc + regimp + log(noposnd) + age + extloc:age
fit.final <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + extloc:age, data = data.vs, method = 'breslow')
fit.final
fit.final$loglik[2]



#----------------------------------------#
#                                        #
#    Prediction of survival function     #
#                                        #
#----------------------------------------#

# median risk: 9.241
med.risk = median(fit.final$linear.predictors) 

# baseline hazard function
breslow <- function(t){
  points = basehaz(fit.final)$time
  basehaz(fit.final)$hazard[length(points[points <= t])]
}

# 1 year survival rate of median risk
exp(-breslow(365)) ^ exp(med.risk)

# 5 year survival rate of median risk
exp(-breslow(365*5)) ^ exp(med.risk)

# Low-risk patient
(risk = predict(fit.final, 
                newdata = data.frame(treat = "3", extloc = "1", regimp = "0", noposnd = 1, age = 50, 
                                     time_dea = 365, status_dea = 0, log_noposnd = 0), type = "lp"))

  # 1 year survival rate of low risk
  predict(fit.final, 
          newdata = data.frame(treat = "3", extloc = "1", regimp = "0", noposnd = 1, age = 50, 
                              time_dea = 365, status_dea = 0, log_noposnd = 0), type = "survival")
  
  # 5 year survival rate of low risk
  predict(fit.final, 
          newdata = data.frame(treat = "3", extloc = "1", regimp = "0", noposnd = 1, age = 50, 
                               time_dea = 365*5, status_dea = 0, log_noposnd = 0), type = "survival")

# High-risk patient
(risk = predict(fit.final, 
                newdata = data.frame(treat = "1", extloc = "4", regimp = "1", noposnd = 33, age = 50, 
                                     time_dea = 365, status_dea = 0, log_noposnd = log(33)), type = "lp"))
  
  # 1 year survival rate of low risk
  predict(fit.final, 
          newdata = data.frame(treat = "1", extloc = "4", regimp = "1", noposnd = 33, age = 50, 
                               time_dea = 365, status_dea = 0, log_noposnd = log(33)), type = "survival")
  
  # 5 year survival rate of low risk
  predict(fit.final, 
          newdata = data.frame(treat = "1", extloc = "4", regimp = "1", noposnd = 33, age = 50, 
                               time_dea = 365*5, status_dea = 0, log_noposnd = log(33)), type = "survival")
  

  
#------------------------------------------#
#                                          #
#  Adjusted estimation of treatment effect #
#                                          #
#------------------------------------------#

  # Cox PH PMLE
  fit.final <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + extloc:age, data = data.vs, method = 'breslow')
  summary(fit.final)
  fit.final$loglik
  
  # Log rank tests #
  Lev.Obs.fit.dea <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + extloc:age, data = data.vs %>% filter(treat !=3), method = 'breslow')
  summary(Lev.Obs.fit.dea)
  (-0.752)**2
  
  FLe.Obs.fit.dea <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + extloc:age, data = data.vs %>% filter(treat !=2), method = 'breslow')
  summary(FLe.Obs.fit.dea)
  (-3.177)**2
  
  Lev.FLe.fit.dea <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + extloc:age, data = data.vs %>% filter(treat !=1), method = 'breslow')
  summary(Lev.FLe.fit.dea)
  (-2.348)**2

  

#----------------------------------------#
#                                        #
#            Model checking              #
#                                        #
#----------------------------------------#

# How many records are missing?
apply(X = data.vs, MARGIN = 2, function(x) sum(x %in% c(".",NA)))
# -- noposnd are missing in 18 patients

# Imputation of noposnd in data.vs as its median
data.vs$noposnd[is.na(data.vs$noposnd)] = median(data.vs$noposnd, na.rm = T)
data.vs <- data.vs %>% mutate(log_noposnd = log(pmax(noposnd,1)))

# Data export for SAS
write.table(x = data.vs, file = "final_data.txt",row.names = F,quote = F)


#----------------------------------------#
#                                        #
#     Functional Forms of covariates     #
#                                        #
#----------------------------------------#

# Age vs log(Age) vs Age^2
fit.age <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd, data = data.vs, method = 'breslow')

p1 <- ggplot(data = NULL, aes(data.vs$age, fit.age$residuals)) + 
  geom_point() +
  geom_smooth(method = "loess") +
  xlab("Age") + ylab("Martingale residual")

p2 <- ggplot(data = NULL, aes(log(data.vs$age), fit.age$residuals)) + 
  geom_point() +
  geom_smooth(method = "loess") +
  xlab("log(Age)") + ylab("Martingale residual")

p3 <- ggplot(data = NULL, aes((data.vs$age)^2, fit.age$residuals)) + 
  geom_point() +
  geom_smooth(method = "loess") +
  xlab("Age^2") + ylab("Martingale residual")

plot_grid(p1, p2, p3, nrow = 1)


# noposnd vs log(noposnd) vs noposnd^2
fit.noposnd <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + age, data = data.vs, method = 'breslow')

p1 <- ggplot(data = NULL, aes(data.vs$noposnd, fit.age$residuals)) + 
  geom_point() +
  geom_smooth(method = "loess") +
  xlab("Nodal involvement") + ylab("Martingale residual")

p2 <- ggplot(data = NULL, aes(log(data.vs$noposnd), fit.age$residuals)) + 
  geom_point() +
  geom_smooth(method = "loess") +
  xlab("log(Nodal involvement)") + ylab("Martingale residual")

p3 <- ggplot(data = NULL, aes((data.vs$noposnd)^2, fit.age$residuals)) + 
  geom_point() +
  geom_smooth(method = "loess") +
  xlab("Nodal involvement^2") + ylab("Martingale residual")

plot_grid(p1, p2, p3, nrow = 1)


#----------------------------------------------#
#                                              #
#   Martingale Residual plot for final model   #
#                                              #
#----------------------------------------------#

fit.final <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + extloc:age, data = data.vs, method = 'breslow')

# Treatment
ggplot(data = NULL, aes(as.factor(data.vs$treat), fit.final$residuals)) +
  geom_violin() +
  geom_jitter(height = 0, width = 0.1, size = 0.5) +
  geom_boxplot(width = 0.1, alpha = 0.5, col = "blue") +
  xlab("Treatment") + ylab("Martingale residual") +
  scale_x_discrete(labels=c("1" = "Observation", 
                            "2" = "Levamisole",
                            "3" = "5-FU + Lev"))

# extloc
ggplot(data = NULL, aes(as.factor(data.vs$extloc), fit.final$residuals)) +
  geom_violin() +
  geom_jitter(height = 0, width = 0.1, size = 0.5) +
  geom_boxplot(width = 0.1, alpha = 0.5, col = "blue") +
  xlab("Extent of local spread") + ylab("Martingale residual") +
  scale_x_discrete(labels=c("1" = "Submucosa", 
                            "2" = "Muscular",
                            "3" = "Serosa",
                            "4" = "Contiguos"))

# regimp
ggplot(data = NULL, aes(as.factor(data.vs$regimp), fit.final$residuals)) +
  geom_violin() +
  geom_jitter(height = 0, width = 0.1, size = 0.5) +
  geom_boxplot(width = 0.1, alpha = 0.5, col = "blue") +
  xlab("Regional implants") + ylab("Martingale residual") +
  scale_x_discrete(labels=c("0" = "No", 
                            "1" = "Exists"))

# log_noposnd
ggplot(data = NULL, aes(data.vs$log_noposnd, fit.final$residuals)) + 
  geom_point() +
  geom_smooth(method = "loess") +
  xlab("log(Nodal involvement)") + ylab("Martingale residual")

# age
ggplot(data = NULL, aes(data.vs$age, fit.final$residuals)) + 
  geom_point() +
  geom_smooth(method = "loess") +
  xlab("Age") + ylab("Martingale residual")



# Using package

ggcoxzph(cox.zph(fit.final))

ggcoxdiagnostics(fit.final, type = "deviance")

  

#----------------------------------------#
#                                        #
#    Model validation and prediction     #
#                                        #
#----------------------------------------#
  
data.val <- data.2 %>% dplyr::select(treat, extloc, regimp, noposnd, age, time_dea, status_dea)

# How many records are missing?
apply(X = data.val, MARGIN = 2, function(x) sum(x %in% c(".",NA)))
# -- noposnd are missing in 13 patients
# -- extloc and regimp are missing in all (243) patients

# Imputation of noposnd in data.val as its median
data.val$noposnd[is.na(data.val$noposnd)] = median(data.val$noposnd, na.rm = T)
data.val <- data.val %>% mutate(log_noposnd = log(pmax(noposnd,1)))

# Imputation of extloc in data.val using proportional odds model
library(MASS)
m <- polr(as.factor(extloc) ~ noposnd + age, data = data.vs)
summary(m)
pred_prob = predict(m, data.val, type = "probs")

data.val$extloc <- as.character(apply(pred_prob, MARGIN = 1, function(x) sample(1:4, 1, replace = F, prob = x)))

table(data.vs$extloc)
table(data.val$extloc)

detach("package:MASS",  unload=TRUE)

# Imputation of regimp in data.val using logistic regression
m <- glm(as.factor(regimp) ~ extloc + noposnd + age, data = data.vs, family = "binomial")
summary(m)
pred_prob = predict(m, data.val, type = "response")

data.val$regimp <- as.character(sapply(pred_prob, function(x) sample(0:1, 1, replace = F, prob = c(1-x,x))))

table(data.vs$regimp)
table(data.val$regimp)


# Partitioning data.val patients by risk group
fit.final <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + extloc:age, data = data.vs, method = 'breslow')
risk.val <- predict(fit.final, newdata = data.val, type = "lp")
low.risk <- which(risk.val < quantile(risk.val, probs = 1/3))
med.risk <- which(quantile(risk.val, probs = 1/3) <= risk.val & risk.val < quantile(risk.val, probs = 2/3))
high.risk <- which(risk.val >= quantile(risk.val, probs = 2/3))
data.val$risk = "1"
data.val$risk[med.risk] = "2"
data.val$risk[high.risk] = "3"


# prediction of KM curves # = average of the predicted survival curves within each group
predict(fit.final, newdata = data.val, type = "survival")

max.time = max(data.val$time_dea) # 4755

pred_surv = data.frame(t=1:max.time, low_risk = 0, med_risk = 0, high_risk = 0)

for(t in 1:max(data.val$time_dea)){
  pred <- predict(fit.final, newdata = data.val %>% mutate(time_dea = t), type = "survival")
  pred_surv$low_risk[t] = mean(pred[low.risk])
  pred_surv$med_risk[t] = mean(pred[med.risk])
  pred_surv$high_risk[t] = mean(pred[high.risk])
}

# KM plots #
fit.val <- survfit(Surv(time_dea, status_dea) ~ risk, data = data.val)
plot(fit.val, col = 1:3, xlab = "Days", ylab = "Survival Probability")
lines(pred_surv$t, pred_surv$low_risk, lty = 2, col=1)
lines(pred_surv$t, pred_surv$med_risk, lty = 2, col=2)
lines(pred_surv$t, pred_surv$high_risk,lty = 2, col=3)
legend(x = 0, y = .3, title="Observed", legend = c("Low risk", "Medium risk", "High risk"), lty = 1, bty = "n", border = F, col = 1:3)
legend(x = 1000, y = .3, title="Expected", legend = c("Low risk", "Medium risk", "High risk"), lty = 2, bty = "n", border = F, col = 1:3)


#------------------------------------------------------------#
#                                                            #
#   Adjusted estimation of treatment effect on total data    #
#                                                            #
#------------------------------------------------------------#

# Cox PH PMLE on Validation dataset
fit.final.val <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + extloc:age, 
                       data = data.val, method = 'breslow')
summary(fit.final.val)

# Cox PH PMLE on Total dataset
fit.final.total <- coxph(Surv(time_dea, status_dea) ~ treat + extloc + regimp + log_noposnd + age + extloc:age, 
                         data = rbind(data.vs, data.val %>% select(-risk)), method = 'breslow')
summary(fit.final.total)






