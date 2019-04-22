library(data.table)
library(dplyr)
library(stringr)
setwd('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas')
#fread('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/raw_omics/kallisto.gene.TPM.not_cross_norm')
trans=fread('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/raw_omics/transcription_summary_means.csv')
clean.annot=fread('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/annotation/CLEAN_NV_TRINITY_UNIREF_ANNOTATION.csv')
dir.create('./Top300')

top300=function(x){
a=trans %>% arrange(desc(get(x))) %>% head(300) %>% data.table()
fwrite(a,paste0('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/Top300/',x,'_Top300.csv'),row.names = F)
writeLines(a$Nv_gene,paste0('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/Top300/',x,'_Top300.txt'))
}

top300("M1_mean")
top300("M2_mean")
top300("M3_mean")
top300('M4H_mean')
