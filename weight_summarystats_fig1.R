### summary statistics and generating figure 1
# EVS 2/2023

library(tidyverse)
library(rstatix)
library(ggpubr)
library(vegan)
library(MASS)

# get Odds Ratio script
source("myOddsRatio.R")

# get growth data
dat <- read.table(file = "all-growth.txt", sep = "\t", header = TRUE)

## ---- summary stats: CWG 12 mo----

## how many were over the 85th percentile
dat %>% 
  filter(age == "mo12") %>% 
  filter(cwg.perc >= 0.85) %>% 
  count()

## ---- summary stats; WFL-Z  ----

# Z scores
dat %>% 
  group_by(age) %>% 
  get_summary_stats(wflZ, type = "common")

# WHO obesity categories
dat %>% 
  group_by(age) %>% 
  count(wfl.cat)

## ----- PANELED FIGURE; WFL-Z at 4 months predicts WFL-Z at 12 months and FIGURE; change in WFL-Z from birth to 12 months compared to WHO category----

# order the dataframe to plot properly and re-name
dat <- dat %>% 
  mutate(wfl.cat = recode_factor(wfl.cat,
                                 severe.under = "Severe Underweight",
                                 under = "Underweight",
                                 normal = "Normal",
                                 over = "Overweight",
                                 obese = "Obese")) %>% 
  mutate(wfl.cat = factor(wfl.cat, ordered = TRUE, levels = c("Severe Underweight", "Underweight",
                                                              "Normal", "Overweight", "Obese")))

## assign category at 12 months to all other ages 
forplot <- dat %>% 
  filter(age == "mo12") %>% 
  mutate(cat.12mo = wfl.cat) %>% 
  dplyr::select(ID, cat.12mo) %>% 
  full_join(dat) %>% 
  ## change age and order
  mutate(age = recode_factor(age,
                             birth = "Birth",
                             mo1 = "1",
                             mo4 = "4",
                             mo6 = "6",
                             mo12 = "12",
                             mo24 = "24")) %>% 
  mutate(age = factor(age, ordered = TRUE, levels = c("Birth", "1", "4", "6", "12", "24")))

### Change in WFL-Z over time
pa <- ggline(forplot %>% filter(!age == "24"), x = "age", y = "wflZ", group = "cat.12mo", color = "cat.12mo",
       # add error bars to the line plot
       add = "mean_se",
       # clean up axis and legend titles
       xlab = "Age (months)", ylab = "WFL", legend.title = "Category at 12 months",
       legend = "right",
       size = 1) +
  # add horizontal line at 0 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  theme(text = element_text(size = 18))

# calculate change WFlZ 12 months to birth
forplot <- dat %>% 
  ## get category at 12 months
  filter(age == "mo12") %>% 
  mutate(cat.12mo = wfl.cat) %>% 
  dplyr::select(ID, cat.12mo) %>% 
  full_join(dat)%>% 
  # calculate change from 12 mo
  dplyr::select(ID, age, wflZ, cat.12mo) %>% 
  pivot_wider(names_from = "age", values_from = "wflZ") %>% 
  mutate(delta12 = mo12 - birth)

# plot correlation
pb <- ggscatter(forplot, x = "mo12", y = "delta12", color = "cat.12mo",
          legend = "right", legend.title = "Category at 12 months",
          xlab = "WFL 12 months", ylab = expression(paste(Delta, "WFL")),
          size = 1.5) +
  # vertical and horizontal lines at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  theme(text = element_text(size = 18))

### combine
ggarrange(pa, pb, labels = c("A.", "B."), common.legend = TRUE, legend = "bottom")


# save
ggsave(file = "2023-Analyses/prelim-figures/paneled-wfl.png", dpi = 600,
       height = 7, width = 12, units = "in")
