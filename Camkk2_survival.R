##Loading from Park SJ, Yoon BH, Kim SK*, Kim SY*. GENT2: an updated gene expression database for normal and tumor tissues.
##BMC Med Genomics. 2019 Jul 11;12(Suppl 5):101. doi: 10.1186/s12920-019-0514-7.

library(readr)
dat <- read_table2("CaMKK2_expression.subtype_GENT2.txt", col_names = T)


##Removing data that doesn't have subtype annotation
library(tidyverse)
dat <- dat %>% filter(Subtype != "None") %>% mutate(Subtype = as.factor(Subtype)) 

#fitting anova
av.fit <- dat %>%  
  aov(Gene_expression ~ Subtype, data = .)
summary(av.fit)

#plotting CaMKK2 v Tumor Grade with Anova and t.test of all v. Grade I
dat %>% ggplot(aes(Subtype,Gene_expression)) + 
  geom_boxplot() + geom_jitter(position = position_jitter(.2)) + 
  stat_compare_means(method = "anova") + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "I", label.y = 12, size = 12) +
  xlab("Brain Tumor Grade") + ylab("log2(CaMKK2 Expression)") + theme_classic2()

#loading packages for survival plotting
install.packages("survminer")
library(survival)
library(survminer)

##subsetting for survival anlaysis

gbm <- dat %>% filter(Subtype == "IV" & Prognosis != "None") %>% 
  mutate(CaMKK2 = ifelse(Gene_expression > median(Gene_expression), "High CaMKK2", "Low CaMKK2")) %>% 
  mutate(CaMKK2 = as.character(CaMKK2)) %>% mutate(Prognosis = ifelse(Prognosis == "Alive",1,2)) %>% 
  mutate(Prognosis = as.numeric(Prognosis)) %>%
  mutate(Prognosis_time = as.numeric(Prognosis_time))
table(gbm$CaMKK2) 

bt2 <- dat %>% filter(Subtype == "II" & Prognosis != "None") %>% 
  mutate(CaMKK2 = ifelse(Gene_expression > median(Gene_expression), "High CaMKK2", "Low CaMKK2")) %>% 
  mutate(CaMKK2 = as.character(CaMKK2)) %>% mutate(Prognosis = ifelse(Prognosis == "Alive",1,2)) %>% 
  mutate(Prognosis = as.numeric(Prognosis)) %>%
  mutate(Prognosis_time = as.numeric(Prognosis_time))
table(bt2$CaMKK2)

bt3 <- dat %>% filter(Subtype == "III" & Prognosis != "None") %>% 
  mutate(CaMKK2 = ifelse(Gene_expression > median(Gene_expression), "High CaMKK2", "Low CaMKK2")) %>% 
  mutate(CaMKK2 = as.character(CaMKK2)) %>% mutate(Prognosis = ifelse(Prognosis == "Alive",1,2)) %>% 
  mutate(Prognosis = as.numeric(Prognosis)) %>%
  mutate(Prognosis_time = as.numeric(Prognosis_time))
table(bt3$CaMKK2) 

##plotting grade IV (GBM)

install.packages("survival")

sfit <- survfit(Surv(time = Prognosis_time,event = Prognosis) ~ CaMKK2, data = gbm)
summary(sfit)
sfit
plot(sfit)

ggsurvplot(sfit, pval = T, linetype = "strata",
           legend.labs = c("> Median Log2(CaMKK2) n = 15", "< Median Log2(CaMKK2) n = 17"),
           title = "Grade IV Brain Tumor (GBM)")

surv_diff <- survdiff(Surv(Prognosis_time, Prognosis)~ CaMKK2, data = gbm)
surv_diff

## Plotting Grade II

sfit2 <- survfit(Surv(time = Prognosis_time,event = Prognosis) ~ CaMKK2, data = bt2)
summary(sfit2)
sfit2

ggsurvplot(sfit2, pval = T, linetype = "strata", 
           legend.labs = c("> Median Log2(CaMKK2) n = 2", "< Median Log2(CaMKK2) n = 3"),
           title = "Grade II Brain Tumor")

surv_diff2 <- survdiff(Surv(Prognosis_time, Prognosis)~ CaMKK2, data = bt2)
surv_diff2

## plotting Grade III

sfit3 <- survfit(Surv(time = Prognosis_time,event = Prognosis) ~ CaMKK2, data = bt3)
summary(sfit3)
sfit3

ggsurvplot(sfit3, pval = T, linetype = "strata", 
           legend.labs = c("> Median Log2(CaMKK2) n = 6", "< Median Log2(CaMKK2) n = 7"),
           title = "Grade III Brain Tumor")

surv_diff3 <- survdiff(Surv(Prognosis_time, Prognosis)~ CaMKK2, data = bt3)
surv_diff3 + ggtitle("Grade III Brain Tumor")
