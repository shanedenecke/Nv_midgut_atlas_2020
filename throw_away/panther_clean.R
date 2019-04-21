library(tidyr)

a=fread('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/annotation/names.tab',header=F,col.names = c('code_name','description'))

mag=a[grep("mag",code_name)]

final=separate(mag,col=code_name,into=c('code_name','junk1','junk2'),sep='\\.') %>% select(code_name,description)

fwrite(final,'./annotation/Panther_clean_key.csv',row.names = F)
