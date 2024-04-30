## PCAs (Figure 2)
# EVS 3/2023

library(tidyverse)
library(phyloseq)
library(microViz)
library(RColorBrewer)
library(ggpubr)

## ---- make one phyloseq with all MBM for birth, 1 mo, 4 mo----

# using miRNA count data
mirna <- read.table(paste0(PATH, "/my-path.csv"), sep = ",", header = TRUE) %>% 
  rename(ID = Participant,
         age = Timepoint) %>% 
  # fix discrepancy in age categories
  mutate(age = case_when(
    str_detect(age, "1") ~ "mo1",
    str_detect(age, "4") ~ "mo4",
    str_detect(age, "Enroll") ~ "birth"
  ))  %>% 
  # make names based on timepoint
  mutate(IDage = paste(ID, age, sep = "_")) %>% 
  dplyr::select(-c(ID, age, Sample.ID)) %>% 
  column_to_rownames(var = "IDage") %>% 
  # remove empty miRNAs (all 0 count)
  dplyr::select_if(colSums(.) > 0)

## get weight outcome data
metadata <- read.table(file = paste0(PATH, "growth.txt"), sep = "\t", header = TRUE)  %>% 
  filter(age %in% c("birth", "mo1", "mo4", "mo12")) %>% 
  mutate(IDage = paste(ID, age, sep = "_")) %>% 
  column_to_rownames(var = "IDage")


## make phyloseq

ps <- phyloseq(otu_table(mirna, taxa_are_rows = FALSE))

# subset metadata
## add 12 month infant weight status
mo12 <- metadata %>% filter(age == "mo12") %>% 
  #mutate(cwg85 = if_else(cwg.perc> 0.85, TRUE, FALSE)) %>% 
  mutate(wflQuant = round(pnorm(wflZ), 3)) %>% 
  dplyr::select(ID, wfl.cat, cwg.perc, wflZ, wflQuant) %>% 
  rename(wfl.cat12 = wfl.cat,
         cwg12 = cwg.perc,
         wfl12 = wflZ)

sampdf <- metadata %>% 
  filter(rownames(metadata) %in% rownames(mirna)) %>% 
  rownames_to_column(var = "colID") %>% 
  left_join(mo12, by = "ID") %>% 
  column_to_rownames(var = "colID")

# make phyloseq
tax <- as.matrix(data.frame(mirna = names(mirna),
                            id = names(mirna)) %>% 
                   column_to_rownames(var = "id"))
ps <- phyloseq(otu_table(mirna, taxa_are_rows = FALSE),
               sample_data(sampdf),
               tax_table(tax))

## beautify 
ps <- ps %>% 
  ps_mutate(Age = case_when(
    age %in% "birth" ~ "Birth",
    age %in% "mo1" ~ "1",
    age %in% "mo4" ~ "4"
  ),
  mom = case_when(
    mom.cat %in% "under" ~ "Underweight",
    mom.cat %in% "normal" ~ "Normal",
    mom.cat %in% "overweight" ~ "Overweight",
    mom.cat %in% "obese" ~ "Obese")
  ) %>% 
  ps_mutate(Age = factor(Age, ordered = TRUE, levels = c("Birth", "1", "4")),
            mom = factor(mom, ordered = TRUE, levels = c(
              "Underweight", "Normal", "Overweight", "Obese"
            )))# %>% 
 # ps_mutate(cwg = if_else(cwg85 == TRUE, "High", "Normal"),
    #        cwg = factor(cwg, ordered = TRUE, levels = c("High", "Normal"))) 

### ---- mbm age ----

pa <- ps %>% 
  tax_transform("clr") %>% 
  #dist_calc("aitchison") %>% 
  ord_calc() %>% 
  ord_plot(color = "Age", shape = "Age",
           auto_caption = NA) +
  scale_color_brewer(palette = "Dark2") +
  stat_ellipse(aes(color = Age)) +
  theme_pubr(legend = "bottom",
             base_size = 16)

## comfirm with test
ps %>% 
  dist_calc("aitchison") %>% 
  dist_permanova(variables = "age",
                 n_processes = 3) # sig (r2 = 6.4%, p = 0.001)

### ---- mom obesity status (birth) ----

## 
pb <- ps %>% 
  ps_filter(age == "birth") %>% 
  tax_transform("clr") %>% 
  #dist_calc("aitchison") %>% 
  ord_calc() %>% 
  ord_plot(color = "mom", 
           auto_caption = NA) +
  scale_color_brewer(palette = "Paired", name = "Maternal BMI") +
  stat_ellipse(aes(color = mom)) +
  theme_pubr(legend = "bottom",
             base_size = 16) 

## comfirm with test
ps %>% ps_filter(age == "birth") %>% 
  dist_calc("aitchison") %>% 
  dist_permanova(variables = "mom.cat",
                 n_processes = 3) 

## ---- infant obesity status (all times) ----

### add 95% percentile WFL-Z
ps <- ps %>% 
  ps_mutate(is.95 = if_else(wflQuant >= 0.95, "High", "Normal"))

##
pc <- ps %>% 
  ps_filter(age == "birth") %>% 
  tax_transform("clr") %>% 
  #dist_calc("aitchison") %>% 
  ord_calc() %>% 
  ord_plot(color = "is.95", 
           auto_caption = NA) +
  scale_color_brewer(palette = "Set2", name = "12 mo Status") +
  stat_ellipse(aes(color = is.95)) +
  theme_pubr(legend = "bottom",
             base_size = 16) 

## comfirm with test
ps %>% 
  ps_filter(age == "birth") %>% 
  dist_calc("aitchison") %>% 
  dist_permanova(variables = "is.95",
                 n_processes = 3) # not sig

## 1 mo
pd <- ps %>% 
  ps_filter(age == "mo1") %>% 
  tax_transform("clr") %>% 
  #dist_calc("aitchison") %>% 
  ord_calc() %>% 
  ord_plot(color = "is.95", 
           auto_caption = NA) +
  scale_color_brewer(palette = "Set2", name = "12 mo Status") +
  stat_ellipse(aes(color = is.95)) +
  theme_pubr(legend = "bottom",
             base_size = 16) 

## comfirm with test
ps %>% 
  ps_filter(age == "mo1") %>% 
  dist_calc("aitchison") %>% 
  dist_permanova(variables = "is.95",
                 n_processes = 3)  # not sig

## mo 4
pe <- ps %>% 
  ps_filter(age == "mo4") %>% 
  tax_transform("clr") %>% 
  #dist_calc("aitchison") %>% 
  ord_calc() %>% 
  ord_plot(color = "is.95", 
           auto_caption = NA) +
  scale_color_brewer(palette = "Set2", name = "12 mo Status") +
  stat_ellipse(aes(color = is.95)) +
  theme_pubr(legend = "bottom",
             base_size = 16) 

## comfirm with test
ps %>% 
  ps_filter(age == "mo4") %>% 
  dist_calc("aitchison") %>% 
  dist_permanova(variables = "is.95",
                 n_processes = 3) # not sig

## add all 
ggarrange(pa, pb, pc, pd, pe,
          labels = c("A.", "B.", "C.", "D.", "E."),
          font.label = list(size = 18))

panela <- ggarrange(pa, pb, ncol = 2, labels = c("A.", "B."),
          font.label = list(size = 18))

panelb <- ggarrange(pc, pd, pe, ncol = 3, nrow = 1,
                    labels = c("C.", "D.", "E."),
                    font.label = list(size = 18))

ggarrange(panela, panelb, ncol = 1)

## save
ggsave(filename = "2023-Analyses/prelim-figures/pcas_biggercaptions.png", dpi = 300,
       height = 10, width = 14, units = "in")

