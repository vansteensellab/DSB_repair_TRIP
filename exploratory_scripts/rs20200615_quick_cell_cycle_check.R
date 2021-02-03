# Quick cell cycle check
library(ggplot2)
# library(ggpubr)
# library(reshape2)
# library(tibble)
# library(GenomicRanges)
# library(rtracklayer)
# library(corrr)
# library(Hmisc)
library(ggbeeswarm)
# library(RColorBrewer)
# library(data.table)
library(dplyr)
# library(plyr)
library(tidyr)
# library(stringr)
# library(plotly)
# library(ggpmisc)
# library(glmnet)
# library(cowplot)
# library(mda)
# library(earth)
# library(yaml)
# library(vip)
# library(caret)
# library(scales)
# library(glmnet)
# library(gridExtra)
# library(ggcorrplot)
# library(bnlearn)
# library(pheatmap)
# library(ppcor)
# library(parallel)
# library(stringr)

# load data
# in.dir.date = 20200410
# in.dir = paste0("/DATA/projects/DSBrepair/data/R/rs", in.dir.date, "/")
# load(paste0(in.dir, "RSTP2_Indel_Chromatin_2kb.RData"))

cell_cyc_tib <- trip_tib_mut_2000 %>% filter(replicate == "XVFUCCIHORep1", plasmid == "LBR2", siRNA == "oxia")


cell_cyc_tib %>% 
  ggplot(., aes(mutation, freq, color = color)) +
  geom_quasirandom() +
  theme_bw(base_size = 16) +
  ylim(0, 1) + facet_grid(ssODN ~ drug)


trip_tib_2000 %>% filter(replicate == "XVFUCCIHORep1", plasmid == "LBR2", siRNA == "oxia") %>% 
  ggplot(., aes(drug, MMEJ_MMEJNHEJ, color = color)) +
  geom_quasirandom() +
  theme_bw(base_size = 16) +
  ylim(0, 1) + facet_grid(~ ssODN)

