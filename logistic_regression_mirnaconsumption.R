### logistic regressions with and without candidate miRNAs
# EVS 3/2023

library(tidyverse)
library(rstatix)
library(ggpubr)
library(MASS)


## get WFL-Z quantiles at or above 95th percentile
weights <- read.table(file = paste0(PATH, "growth.txt"), sep = "\t", header = TRUE)  %>% 
  filter(age == "mo12") %>% 
  #mutate(cwg85 = if_else(cwg.perc> 0.85, TRUE, FALSE)) %>% 
  mutate(wflQuant = round(pnorm(wflZ), 3)) %>% 
  mutate(is.95 = if_else(wflQuant >= 0.95, "High", "Normal")) %>% 
  dplyr::select(ID, age, wfl.cat, cwg, wflZ, wflQuant, is.95) %>% 
  rename(wfl.cat12 = wfl.cat,
         cwg12 = cwg,
         wfl12 = wflZ)

## join data for model
formod <- weights %>% 
  dplyr::select(ID, is.95) %>% 
  #dplyr::select(ID, wfl.cat12) %>% 
  # mom pre-preg BMI
  full_join(bmi %>% dplyr::select(-momBMI_SVN)) %>% 
  # sleep status at 6 months - remove for now (a lot of NAs)
  #full_join(sleep %>% dplyr::select(ID, wk24_waking)) %>% 
  # breastfeeding
  full_join(bf %>% dplyr::select(ID, bf_at6mo)) %>% 
  # demographics 
  full_join(demo %>% 
              dplyr::select(ID, mom_white, mom_ins, preg_gdm)) %>% 
  column_to_rownames(var = "ID") %>% 
  # binary-ize
  mutate(is.95 = if_else(is.95 == "High", 1, 0)) %>% 
  #  mutate(wfl.cat12 = if_else(wfl.cat12 %in% c("obese", "over"), 1, 0)) %>% 
  # DUMMY CODES
  mutate(momNorm = if_else(momBMIcat %in% c("under", "normal"), 1, 0),
         #momOver = if_else(momBMIcat %in% c("overweight", "obese"), 1, 0),
         bf_at6mo = if_else(bf_at6mo == "BF", 1, 0),
         ins1 = if_else(mom_ins == 1, 1, 0)) %>% 
         #ins2 = if_else(mom_ins == 2, 1, 0)) %>% 
         #ins34 = if_else(mom_ins %in% c(3, 4), 1, 0)) %>% 
  dplyr::select(-c(momBMIcat, mom_ins))

## make full model
mod <- glm(is.95 ~ ., data = formod, family = "binomial")
summary(mod) 

## predict probability of 95% percentile
probs <- predict(mod, type = "response")
pred.cl <- ifelse(probs > 0.5, 1, 0)
length(which(pred.cl == 1))  # none are predicted


##check influential values

plot(mod, which = 4, id.n = 3) 

# get model results
model.data <- augment(mod) %>% 
  mutate(index = 1:n())

# plot standardized residuals
ggplot(model.data, aes(index, .std.resid)) +
  geom_point(aes(color = is.95), alpha = 0.5)

# get high st res
model.data %>% 
  filter(abs(.std.resid) > 3) # none

## test multicollinearity

car::vif(mod) # low (none above 4 or under 0.25)


### ----get miRNA data -----

## get candidates
source("2023-Analyses/mirnaconsumption_linearmods.R")

## using miRNA consumption data; normalized ("normdf")
### no ages (this is 'total')

## -----  full model ----

# get data
full <- normdf %>% 
  mutate(ID = as.character(ID)) %>% 
  # subset to cnadidat mirRNAS
  dplyr::select(ID, all_of(cands)) %>% 
  left_join(formod %>% rownames_to_column(var = "ID")) %>% 
  column_to_rownames(var = "ID")

## make reduced model without demo
mod1.5 <- glm(is.95 ~ ., data = full %>% dplyr::select(starts_with("hsa"), is.95), family = "binomial")
summary(mod1.5) # hsa.miR.224.5p
## predict probability of 95% percentile
probs <- predict(mod1.5, type = "response")
pred.cl <- ifelse(probs > 0.5, 1, 0)
length(which(pred.cl == 1))  # 0 predicted
# make matrix of correct and incorrect
pred.df <- full %>% rownames_to_column(var = "ID") %>% dplyr::select(ID, is.95) %>% 
  full_join(data.frame(
    ID = names(pred.cl),
    preds = unname(pred.cl)
  ))
## summarize
length(which(pred.df$is.95 == 1 &  pred.df$preds == 1)) # 0
length(which(pred.df$is.95 == 1 &  pred.df$preds == 0)) # 15
length(which(pred.df$is.95 == 0 &  pred.df$preds == 1)) # 0


## test colinearity
car::vif(mod1.5) # a lot of the miRNAs are colinear
plot(hsa.miR.30a.5p ~ hsa.miR.141.3p, data = full)

##check influential values
plot(mod1.5, which = 4, id.n = 3)

## make full model with demo data
mod2 <- glm(is.95 ~ ., data = full, family = "binomial")
summary(mod2) # hsa.miR.224.5p, none of the demo variables

## predict probability of 95% percentile
probs <- predict(mod2, type = "response")
pred.cl <- ifelse(probs > 0.5, 1, 0)
length(which(pred.cl == 1))  # 1 predicted
# make matrix of correct and incorrect
pred.df <- full %>% rownames_to_column(var = "ID") %>% dplyr::select(ID, is.95) %>% 
  full_join(data.frame(
    ID = names(pred.cl),
    preds = unname(pred.cl)
  ))

## summarize
length(which(pred.df$is.95 == 1 &  pred.df$preds == 1)) # 1
length(which(pred.df$is.95 == 1 &  pred.df$preds == 0)) # 14
length(which(pred.df$is.95 == 0 &  pred.df$preds == 1)) # 0

## test colinearity
car::vif(mod2) # four of the miRNAs are colinear (>4)

##check influential values
plot(mod2, which = 4, id.n = 3) # none with high Cook sd

#### ---- use just miR-224-5p ----

mod3 <- glm(is.95 ~ hsa.miR.224.5p, data = full, family = "binomial")
summary(mod3)


## predict probability of 95% percentile
probs <- predict(mod3, type = "response")
pred.cl <- ifelse(probs > 0.5, 1, 0)
length(which(pred.cl == 1))  # 1 predicted
# make matrix of correct and incorrect
pred.df <- full %>% rownames_to_column(var = "ID") %>% dplyr::select(ID, is.95) %>% 
  full_join(data.frame(
    ID = names(pred.cl),
    preds = unname(pred.cl)
  ))

## summarize
length(which(pred.df$is.95 == 1 &  pred.df$preds == 1)) # 1
length(which(pred.df$is.95 == 1 &  pred.df$preds == 0)) # 14
length(which(pred.df$is.95 == 0 &  pred.df$preds == 1)) # 0

### ---- plot; miR-224-5p ----

library(cowplot)
# make pretty
forplot <- full %>% 
  mutate(fp95 = if_else(is.95 == 1, "> 95th", "< 95th"))

# boxplot
ggboxplot(forplot, x = "fp95", y = "hsa.miR.224.5p", xlab = "CWG percentile at 12 months", ylab = "Normalized miR-224-5p", size = 1)

# save
ggsave(filename = "2023-Analyses/prelim-figures/boxplot_miR2245p.png", dpi = 300)

# plot model predictions
ggplot(data = forplot, aes(x = hsa.miR.224.5p, y = fp95)) +
  geom_jitter(width = 0, height = 0.05) +
  
  geom_smooth(method="glm",  method.args = list(family="binomial")) +
  labs(x = "Normalized miR-224-5p", y = "CWG percentile at 12 months") +
  theme_cowplot() 


