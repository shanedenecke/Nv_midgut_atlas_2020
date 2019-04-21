library(dplyr)
library(data.table)
library(xlsx)

setwd('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/Nv_P450id/expression')
nv.cyp=paste0('TRINITY_',readLines('./Nv_P450s.txt'))
trans=fread('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/raw_omics/transcription_summary_means.csv',
            select=c('Nv_gene','M1_mean','M2_mean','M3_mean','M4H_mean','C_mean'))

trans[Nv_gene %in% nv.cyp] %>% write.xlsx2('./Nv_cyp_expression.xlsx',row.names = F)


nv.protease=paste0('TRINITY_',readLines('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/protease_genes.txt'))
trans[Nv_gene %in% nv.protease]  %>% write.xlsx2('./Nv_protease_expression.xlsx',row.names = F)
