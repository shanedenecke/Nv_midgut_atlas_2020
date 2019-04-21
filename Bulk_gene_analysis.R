##########Perform Pfam enrichment
### Lacking GO terms and updated (filtered) dataset


library(data.table)
library(dplyr)
library(stringr)

setwd('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas')
filtered.genes=readLines('./raw_omics/filtered_gene_list.txt')


### FUNCTIONS

#### CALCULATES FREQUENCEY
freq.table=function(x){
  return(data.table(names=names(x),x) %>% select(-.) %>% mutate(total=sum(N),freq=N/total))
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
  enriched= combined %>% arrange(fdr) %>% select(-freq) %>% mutate(test_freq=test_N/test_total, ref_freq=ref_N/ref_total) %>% 
    select(names,fischP,fdr,everything()) %>% data.table()
  return(enriched)
}




### IMPORT DATA
ipscan=fread('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/annotation/Trinity.fasta.tsv',header=F,fill=T,select=c(1,3,4,5,6),
             col.names = c('Nv_gene','len','database','code','description'))
ipscan$Nv_gene=gsub("_i[0-9].+$","",ipscan$Nv_gene)
ipscan=ipscan[Nv_gene %in% filter
pfam.key=fread('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/annotation/Pfam_key.tsv',header=F,col.names = c('names','annotation'))


#### FILTER DATA
panther=ipscan[database=='PANTHER']
pfam=ipscan[database=='Pfam']


###CALCULATE GENOME WIDE FREQUENCIES
panth.freq=panther$code %>% table() %>% freq.table() %>% data.table()
pfam.freq=pfam$code %>% table() %>% freq.table() %>% data.table()


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



###PERFORM COMPARTMENT SPECIFIC ANALYSIS

#read in data
comp.de.files=list.files('./DE_expression_tables/DE_outputs/',full.names = T)[grepl("specific.txt",list.files('./DE_expression_tables/DE_outputs/'))]
comp.de.list=list()
for(i in comp.de.files){comp.de.list[[i]]=readLines(i)}

#subset gut specific genes from annotations
comp.gene.subset=function(geneset,database){
  return(database[Nv_gene %in% geneset])
}
panth.comp.specific=lapply(comp.de.list,comp.gene.subset,panther)
pfam.comp.specific=lapply(comp.de.list,comp.gene.subset,pfam)


##make frequency tables
panth.comp.specific.freq=lapply(panth.comp.specific, function(x) x$code %>% table() %>% freq.table() %>% data.table())
pfam.comp.specific.freq=lapply(pfam.comp.specific, function(x) x$code %>% table() %>% freq.table() %>% data.table())


##caluclate enrichment
panth.comp.enrich=lapply(panth.comp.specific.freq,vlepo.enrich,panth.freq)
pfam.comp.enrich=lapply(pfam.comp.specific.freq,vlepo.enrich,pfam.freq)

##anotate
pfam.comp.annotate=lapply(pfam.comp.enrich,merge,y=pfam.key,by='names',sort=F)


