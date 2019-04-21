## import packages
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)



## set WD
setwd('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/')
dir.create('./DE_outputs/')
filtered.genes=readLines('./raw_omics/filtered_gene_list.txt')


#### CREATE ANNOTATION (TAKES A WHILE)
annot=fread('/home/shanedenecke/Documents/omics_data/Nezara_viridula/Nv_uniref_blast_annotation_trinity_april_15_2019.txt',
            select=c('V1','V11','V16'),col.names = c('Nv_gene','evalue','annot')) 
l=list()
for(i in unique(annot$Nv_gene)){
  l[[i]]=annot[Nv_gene==i] %>% arrange(evalue) %>% head(1)
}
top.annot=rbindlist(l)

clean.annot=top.annot %>% separate(col=annot,into=c('V2','V3'),sep=c('UniRef50_[:alnum:]+')) %>% separate(col=V3,into=c('annotation','V5'),sep='n=[0-9]+') %>%
  separate(V5,into=c('V6','V7'),sep='Tax=') %>%  separate(V7,into=c('taxa','V9'),sep=' ') %>% select(Nv_gene,annotation,evalue,taxa) %>% unite(col='annotation',annotation,taxa,evalue)
clean.annot$Nv_gene=paste('TRINITY',clean.annot$Nv_gene,sep='_')
fwrite(clean.annot,'./annotation/CLEAN_NV_TRINITY_UNIREF_ANNOTATION.csv')

### READ ANNOTATION
clean.annot=fread('./annotation/CLEAN_NV_TRINITY_UNIREF_ANNOTATION.csv')

## add proper row names to each file. First is blank by default for some reason
for(i in list.files('./DE_expression_tables/',full.names = T)[grepl(".csv",list.files('./DE_expression_tables/'))]){
  a=fread(i)
  colnames(a)[1]='Nv_gene'
  b=a %>% rowwise() %>% mutate(M1_mean=mean(c(M1_1,M1_2,M1_3,M1_4)),M1_sd=sd(c(M1_1,M1_2,M1_3,M1_4))) %>% 
    mutate(M2_mean=mean(c(M2_1,M2_2,M2_3,M2_4)),M2_sd=sd(c(M2_1,M2_2,M2_3,M2_4))) %>% 
    mutate(M3_mean=mean(c(M3_1,M3_2,M3_3,M3_4)),M3_sd=sd(c(M3_1,M3_2,M3_3,M3_4))) %>%
    mutate(M4H_mean=mean(c(M4H_1,M4H_2,M4H_3,M4H_4)),M4H_sd=sd(c(M4H_1,M4H_2,M4H_3,M4H_4))) %>%
    mutate(C_mean=mean(c(C1,C2,C3,C4)),C_sd=sd(c(C1,C2,C3,C4))) %>% 
    filter(Nv_gene %in% filtered.genes) %>% 
    data.table()
  fwrite(b,i)
}



#### Carcass vrs gut comparison

gut.sp=list.files('./DE_expression_tables/',full.names = T)[grepl("C__",list.files('./DE_expression_tables/')) & grepl("M[0-9]-UP",list.files('./DE_expression_tables/'))]
gut.sp.tl=list()
gut.sp.vl=list()

for(i in gut.sp){
  de=fread(i)
  gut.sp.vl[[i]]=de$Nv_gene
  gut.sp.tl[[i]]=de
}
gut.specific.list=Reduce(intersect,gut.sp.vl)
writeLines(gut.specific.list,'./DE_expression_tables/DE_outputs/Common_gut_specific_genes.txt')
gut.specific.table=lapply(gut.sp.tl,function(x) x[Nv_gene %in% gut.specific.list])[[1]] %>% 
  select(Nv_gene,M1_mean,M2_mean,M3_mean,M4H_mean,C_mean) %>% 
  merge(clean.annot,by='Nv_gene')
fwrite(gut.specific.table,'./DE_expression_tables/DE_outputs/Common_gut_specific_table.csv',row.names = F)


#### Specific to each compartment comparsion
compartment.list=c('M1','M2','M3','M4')
for(i in compartment.list){
  comp.up=paste0(i,'-UP.csv')
  temp.list=list()
  for(j in list.files('./DE_expression_tables/',full.names = T)[grep(comp.up,list.files('./DE_expression_tables/'))]){ ### loop over lists where each compartment is upregulated
    temp.list[[j]]=fread(j)$Nv_gene
  }
  comp.common=Reduce(intersect,temp.list)
  writeLines(comp.common,paste0('./DE_expression_tables/DE_outputs/',i,'_specific.txt'))
  comp.common.table=fread(j)[Nv_gene %in% comp.common] %>% select(Nv_gene,M1_mean,M2_mean,M3_mean,M4H_mean,C_mean) %>% merge(clean.annot,by='Nv_gene')
  fwrite(comp.common.table,paste0('./DE_expression_tables/DE_outputs/',i,'_specific.csv'),row.names = F)
}
    






###################### BIN


l=list()
for(i in list.files()){
  a=fread(paste0(getwd(),"/",i))
  colnames(a)[1]='Nv_gene'
  one=a[1,'sampleA'] %>% as.character()
  two=a[1,'sampleB'] %>% as.character()
  b=a %>% select(Nv_gene:FDR,matches(one),matches(two)) ## select relvant columns including raw values
  c=select(b,Nv_gene:FDR) ## remove raw values 
  up=gsub('.csv','',unlist(strsplit(i,split='__'))[3])
  c$up=up
  l[[i]]=c ## add to list
}

dir.create('~/Dropbox/WP2_omics/Helicoverpa_spatial_gut_atlas/clean_de_tables')
lapply(l,function(x) write.csv(x,file=paste('~/Dropbox/WP2_omics/Helicoverpa_spatial_gut_atlas/clean_de_tables/',x[1,'sampleA'],"__",x[1,'sampleB'],"__",x[1,'up'],"_UP.csv"),row.names = F))





###Bin


#m.gene=function(x){
#  merge(x,fread('~/Documents/omics_data/Helicoverpa_armigera/keys/updated_key_jan_2019_by_gene.csv',select=c('ha_geneid','annot','Geneid')),by='Geneid',all.x=T)
#}
