##########Perform Pfam enrichment
### Lacking GO terms


library(data.table)
library(dplyr)
library(stringr)

setwd('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas')
dir.create('Bulk_gene_analysis')
filtered.genes=readLines('./raw_omics/filtered_gene_list.txt')


### FUNCTIONS

#### CALCULATES FREQUENCEY
freq.table=function(x){
  total=x$Nv_gene %>% unique() %>% length()
  l=list()
  for(i in unique(x$code)){
    sub=x[code==i]
    num.genes=unique(sub$Nv_gene) %>% length()
    freq=num.genes/total
    l[[i]]=data.table(code_name=i,total=num.genes,frequency=freq)
  }
  return(rbindlist(l))
} 


### function which takes frequency tables for test gene set and genome wide set (obtained through freq.table function
vlepo.enrich=function(x,y){
  fischer_p=c()
  test.list=list()
  ref.list=list()
  for(i in x$names){
    test=x[x$names==i,c('N','total')] %>% as.numeric()
    ref=y[y$names==i,c('N','total')] %>% as.numeric()
    p.value=fisher.test(matrix(c(test[1],ref[1],test[2],ref[2]),nrow=2,ncol=2),alternative="greater")$p.value ## 3extract P value from fischer test
    fischer_p=c(fischer_p,p.value)
    
    names(test)=c('test_N','test_total')
    test.list[[i]]=as.data.frame(t(test))
    
    names(ref)=c('ref_N','ref_total')
    ref.list[[i]]=as.data.frame(t(ref))
    
  }
  
  test.add=rbindlist(test.list)
  ref.add=rbindlist(ref.list)
  x$fischP=fischer_p 
  x$fdr=p.adjust(x$fischP,method='fdr')
  combined=cbind(test.add,ref.add,x)
  enriched= combined %>% arrange(fdr) %>% filter(fdr<1e-05) %>% select(-freq,-N,-total) %>% mutate(test_freq=test_N/test_total, ref_freq=ref_N/ref_total) %>% 
    select(names,fischP,fdr,everything()) %>% data.table()
  return(enriched)
}


###subset gut specific genes from annotations
gene.subset=function(geneset,database){
  return(database[Nv_gene %in% geneset])
}


### IMPORT DATA
ipscan=fread('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/annotation/Trinity.fasta.tsv',header=F,fill=T,select=c(1,4,5,6),
             col.names = c('Nv_gene','database','code','description'))
ipscan$Nv_gene=gsub("_i[0-9].+$","",ipscan$Nv_gene)
ipscan.filter=ipscan[Nv_gene %in% filtered.genes]
pfam.key=fread('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/annotation/Pfam_key.tsv',header=F,col.names = c('names','annotation'))


#### FILTER DATA
panther=ipscan.filter[database=='PANTHER'] %>% unique.data.frame()
pfam=ipscan.filter[database=='Pfam'] %>% unique.data.frame()


###CALCULATE GENOME WIDE FREQUENCIES
panth.freq=panther %>% freq.table() %>% data.table()
pfam.freq=pfam %>% freq.table() %>% data.table()


#### PERFORM GROUP ANLAYSIS

gut.specific.genes='./DE_expression_tables/DE_outputs/Common_gut_specific_genes.txt'
comp.de.files=list.files('./DE_expression_tables/DE_outputs/',full.names = T)[grepl("specific.txt$",list.files('./DE_expression_tables/DE_outputs/'))]
comp.fuzz.files=list.files('./Mfuzz',full.names = T)[grepl("txt$",list.files('./Mfuzz'))]

files=c(gut.specific.genes,comp.de.files,comp.fuzz.files)

comparisons.list=list()
for(i in files){comparisons.list[[i]]=readLines(i)}

panth.specific=lapply(comparisons.list,gene.subset,panther)
pfam.specific=lapply(comparisons.list,gene.subset,pfam)



##make frequency tables
panth.specific.freq=lapply(panth.specific, freq.table())
pfam.specific.freq=lapply(pfam.specific, function(x) x$code %>% table() %>% freq.table() %>% data.table())


##caluclate enrichment 
panth.enrich=lapply(panth.specific.freq,vlepo.enrich,panth.freq)
pfam.enrich=lapply(pfam.specific.freq,vlepo.enrich,pfam.freq)



##anotate PFAM file
annotate_enrichment=function(enrich.table,key){
  l=list()
  for(i in names(enrich.table)){
    enrich=enrich.table[[i]]
    annot=merge(enrich,key,by='names',all.y=F,sort=F)
    l[[i]]=(annot[!duplicated(annot)])
  } 
  return(l)
}
pfam.annotate=annotate_enrichment(pfam.enrich,pfam.key)



### write.files

for(i in names(panth.enrich)){
  short=strsplit(i,'\\/')[[1]][length(strsplit(i,'\\/')[[1]])]
  fwrite(panth.enrich[[i]],paste0('./Bulk_gene_analysis/',short,'_PANTHER_enrichment.csv'))
  }
for(i in names(pfam.annotate)){
  short=strsplit(i,'\\/')[[1]][length(strsplit(i,'\\/')[[1]])]
  fwrite(pfam.annotate[[i]],paste0('./Bulk_gene_analysis/',short,'_PFAM_enrichment_annotate.csv'))
  }









############### BIN








###PERFORM COMPARTMENT SPECIFIC ANALYSIS

#read in data
comp.de.files=list.files('./DE_expression_tables/DE_outputs/',full.names = T)[grepl("specific.txt",list.files('./DE_expression_tables/DE_outputs/'))]
comp.de.list=list()
for(i in comp.de.files){comp.de.list[[i]]=readLines(i)}

panth.comp.specific=lapply(comp.de.list,comp.gene.subset,panther)
pfam.comp.specific=lapply(comp.de.list,comp.gene.subset,pfam)


##make frequency tables
panth.comp.specific.freq=lapply(panth.comp.specific, function(x) x$code %>% table() %>% freq.table() %>% data.table())
pfam.comp.specific.freq=lapply(pfam.comp.specific, function(x) x$code %>% table() %>% freq.table() %>% data.table())


##caluclate enrichment #### WRITE files
panth.comp.enrich=lapply(panth.comp.specific.freq,vlepo.enrich,panth.freq)
for(i in names(panth.comp.enrich)){fwrite(panth.comp.enrich[[i]],paste0(i,'_PANTHER_enrichment.csv'))}

pfam.comp.enrich=lapply(pfam.comp.specific.freq,vlepo.enrich,pfam.freq)

##anotate
annotate_enrichment=function(enrich.table,key){
  l=list()
  for(i in names(enrich.table)){
    enrich=enrich.table[[i]]
    annot=merge(enrich,key,by='names',all.y=F,sort=F)
    l[[i]]=(annot[!duplicated(annot)])
  } 
  return(l)
}

pfam.comp.annotate=annotate_enrichment(pfam.comp.enrich,pfam.key)
for(i in names(pfam.comp.annotate)){fwrite(pfam.comp.annotate[[i]],paste0(i,'_PFAM_enrichment_annotate.csv'))}






###Mfuzz

#read in data
comp.fuzz.files=list.files('./Mfuzz',full.names = T)[grepl("txt",list.files('./Mfuzz'))]
comp.fuzz.list=list()
for(i in comp.fuzz.files){comp.fuzz.list[[i]]=readLines(i)}

#subset gut specific genes from annotations
comp.gene.subset=function(geneset,database){
  return(database[Nv_gene %in% geneset])
}
panth.comp.specific=lapply(comp.fuzz.list,comp.gene.subset,panther)
pfam.comp.specific=lapply(comp.fuzz.list,comp.gene.subset,pfam)


##make frequency tables
panth.comp.specific.freq=lapply(panth.comp.specific, function(x) x$code %>% table() %>% freq.table() %>% data.table())
pfam.comp.specific.freq=lapply(pfam.comp.specific, function(x) x$code %>% table() %>% freq.table() %>% data.table())


##caluclate enrichment #### WRITE files
panth.comp.enrich=lapply(panth.comp.specific.freq,vlepo.enrich,panth.freq)
for(i in names(panth.comp.enrich)){fwrite(panth.comp.enrich[[i]],paste0(i,'_PANTHER_enrichment.csv'))}

pfam.comp.enrich=lapply(pfam.comp.specific.freq,vlepo.enrich,pfam.freq)

##anotate
annotate_enrichment=function(enrich.table,key){
  l=list()
  for(i in names(enrich.table)){
    enrich=enrich.table[[i]]
    annot=merge(enrich,key,by='names',all.y=F,sort=F)
    l[[i]]=(annot[!duplicated(annot)])
  } 
  return(l)
}

pfam.comp.annotate=annotate_enrichment(pfam.comp.enrich,pfam.key)
for(i in names(pfam.comp.annotate)){fwrite(pfam.comp.annotate[[i]],paste0(i,'_PFAM_enrichment_annotate.csv'))}





################ BIN



##PERFORM GUT SPECIFIC ANALYSIS

#read in gut specific data
gut.specific.genes=readLines('./DE_expression_tables/DE_outputs/Common_gut_specific_genes.txt')

#subset gut specific genes from annotations
panth.gut.specific=panther[Nv_gene %in% gut.specific.genes]
pfam.gut.specific=pfam[Nv_gene %in% gut.specific.genes]

## make test frequency tables
panth.gut.specific.freq=panth.gut.specific$code %>% table() %>% freq.table() %>% data.table()
pfam.gut.specific.freq=pfam.gut.specific$code %>% table() %>% freq.table() %>% data.table()


##caluclate enrichment
panth.gut.enrich=vlepo.enrich(panth.gut.specific.freq,panth.freq)
pfam.gut.enrich=vlepo.enrich(pfam.gut.specific.freq,pfam.freq)





