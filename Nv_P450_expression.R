library(dplyr)
library(data.table)
library(xlsx)

setwd('/data2/shane/Documents/Nezara_midgut_atlas/Nv_P450')
trans=fread('../raw_omics/transcription_summary_means.csv',
            select=c('Nv_gene','M1_mean','M2_mean','M3_mean','M4H_mean','C_mean'))

nv.cyp=paste0('TRINITY_',readLines('./Nv_P450s.txt'))

nv.slcs=fread('/data2/shane/Documents/SLC_id/final_SLC_dicts/NezVirFinal_SLC_table.csv')
slc.sub=nv.slcs[grepl('SLC_5_|SLC_6_|SLC_7_|SLC_15_|SLC_36_|SLC_2_|SLC_45_|SLC_50_',name)]
slc.sub$code=paste0('TRINITY_',slc.sub$code)
colnames(slc.sub)[1]='Nv_gene'
merged=merge(trans,slc.sub,by='Nv_gene') %>% separate(name,into=c('good1','good2',sep='_')) %>% unite(col='SLC',good1,good2,sep='_') 
fwrite(merged,'./expression/SLC_expression.csv',row.names = F)



trans[Nv_gene %in% nv.cyp] %>% write.xlsx2('./Nv_cyp_expression.xlsx',row.names = F)


nv.protease=paste0('TRINITY_',readLines('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/protease_genes.txt'))
trans[Nv_gene %in% nv.protease]  %>% write.xlsx2('./Nv_protease_expression.xlsx',row.names = F)
