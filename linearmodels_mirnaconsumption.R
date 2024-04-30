## linear regression to get candidate miRNAs
# from consumption data
# EVS 12/2023

library(tidyverse)
library(rstatix)
library(ggpubr)
library(car) # for box cox transformation

#get data
## CHANGE THIS PATHWAY TO THE PATH TO THE DATA ON YOUR COMPUTER
PATH <- "my-path"

# get miRNA data by consumption
dat <- read.csv("../mirnas.csv", sep = ",", header = TRUE)
mir <- dat %>% dplyr::select(ID, starts_with("hsa")) %>% 
  column_to_rownames(var = "ID")

## normalize for linear regression; Box cox?
normdf <- data.frame(ID = row.names(mir))
for(j in 1:ncol(mir)) {
  
  # get lambda
  pt <- powerTransform(mir[,j], family = "bcPower")
  lambda <- unname(pt$lambda)
  
  # transform and add to new dataframe
  normdf[,j+1] <- ((mir[,j] ^ lambda) - 1) / lambda
  colnames(normdf)[j+1] <- names(mir)[j]
}


## ----- linear regression for each ----

## get weight outcome data
metadata <- read.table(file = paste0(PATH, "growth.txt"), sep = "\t", header = TRUE)  %>% 
  #filter(age %in% c("birth", "mo1", "mo4", "mo12")) %>% 
  filter(age %in% "mo12")


# join
withmeta <- normdf %>% 
  mutate(ID = as.numeric(ID)) %>% 
  left_join(metadata)

# do linear regression in a loop
hsas <- names(withmeta %>% dplyr::select(starts_with("hsa")))
outdf <- data.frame()
for(i in 1:length(hsas)) {
  
  # 1. CWG at 12 months
  mod1 <- lm(withmeta[,i+1] ~ cwg.perc, data = withmeta)
  # get summary stats
  m1s <- data.frame(
    mirna = hsas[i],
    adjr2 = summary(mod1)$adj.r.squared,
    est = round(summary(mod1)$coefficients[2,1], 3),
    pval = round(summary(mod1)$coefficients[2,4], 3),
    mod = "CWG"
  )
  # 2. WFL Z-score
  mod2 <- lm(withmeta[,i+1] ~ wflZ, data = withmeta)
  m2s <- data.frame(
    mirna = hsas[i],
    adjr2 = summary(mod2)$adj.r.squared,
    est = round(summary(mod2)$coefficients[2,1], 3),
    pval = round(summary(mod2)$coefficients[2,4], 3),
    mod = "WFLZ"
  )
  
  # save output
  outdf <- rbind(outdf, m1s, m2s)
}

## adjust p values
cwg <- outdf %>% 
  filter(mod == "CWG") %>% 
  mutate(padj = p.adjust(pval, method = "fdr"))

wfl <- outdf %>% 
  filter(mod == "WFLZ") %>% 
  mutate(padj = p.adjust(pval, method = "fdr"))

# none are significant after multiple comparisons - get significant raw p values
sigs <- outdf %>% 
  filter(pval < 0.05)

# get values for both models 
outdf %>% filter(mirna %in% sigs$mirna) %>% dplyr::select(-est) %>% pivot_wider(names_from = mod, values_from = c(adjr2, pval)) %>% dplyr::select(mirna, adjr2_CWG, pval_CWG, adjr2_WFLZ, pval_WFLZ)

unique(sigs$mirna) # 8

# save this
cands <- unique(sigs$mirna)

