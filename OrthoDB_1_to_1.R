a=fread('/home/shanedenecke/Dropbox/wp4_drug_target_search/target_identify/nezara/input/ortho.table_raw.csv',select=c('Nv_prot','dm_FBgn'),col.names = c('Nv_prot','Dm_FBgn'))


### dm key filter by gene
dm.key=fread('/home/shanedenecke/Documents/omics_data/Drosophila_melanogaster/keys/Dm_master_key_FB_fasta.csv',select=c('Dm_FBgn','CG','name','len'))
dm2=arrange(dm.key,desc(len)) %>% arrange(name) %>% data.table()
dm3=dm2[!duplicated(Dm_FBgn)]
fwrite(dm3,'/home/shanedenecke/Documents/omics_data/Drosophila_melanogaster/keys/Dm_master_key_by_gene.csv',row.names = F)
dm.tf=fread('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/annotation/Transcription_factor_list.csv')

old.trans=fread('/home/shanedenecke/Documents/omics_data/Bayer_transcriptomes_clean/Nezara_SPKM.csv')

b=a[!grepl("_",Dm_FBgn) & grepl("FB",Dm_FBgn)]


dm.names=table(b$Dm_FBgn)[table(b$Dm_FBgn)==1] %>% names()
c=b[Dm_FBgn %in% dm.names]
merge(dm3,c,by='Dm_FBgn') %>% select(-len)

e=c[table(c$Dm_FBgn)==1] %>% merge(old.trans,by='Nv_prot') %>% merge(dm3,by='Dm_FBgn')

fwrite(e,'/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/annotation/Drosophila_1_to_1_orthologues.csv') 
tf.merge=merge(e,dm.tf,by='Dm_FBgn')
fwrite(tf.merge,'/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/annotation/TF_Drosophila_1_to_1_orthologues.csv') 


dm.names=names(table(c$Nv_prot)==1)

d=c[Nv_prot %in% dm.names]


b[!duplicated(b$dm_FBgn)]

b[!(duplicated(b) | duplicated(b, fromLast = TRUE)), ]

uniqueN(c(1,1,2,2,3))
