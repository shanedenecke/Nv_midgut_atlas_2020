##########Perform Pfam enrichment
### Lacking GO terms


library(data.table)
library(dplyr)
library(stringr)

setwd('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas')
dir.create('Bulk_gene_analysis')
filtered.genes=readLines('./raw_omics/filtered_gene_list.txt')
clean.annot=fread('./annotation/CLEAN_NV_TRINITY_UNIREF_ANNOTATION.csv')


### FUNCTIONS

#### CALCULATES FREQUENCEY
freq.table=function(x){
  total=x$Nv_gene %>% unique() %>% length()
  l=list()
  for(i in unique(x$code)){
    sub=x[code==i]
    num.genes=unique(sub$Nv_gene) %>% length()
    freq=num.genes/total
    l[[i]]=data.table(code_name=i,code_count=num.genes,frequency=freq,total=total)
  }
  return(rbindlist(l))
} 


### function which takes frequency tables for test gene set and genome wide set (obtained through freq.table function

#trial=panth.specific.freq[[1]]
#omic=panther %>% freq.table() %>% data.table()

vlepo.enrich=function(trial,omic){
  fischer_p=c()
  test.list=list()
  ref.list=list()
  for(i in trial$code_name){
    test=trial[trial$code_name==i,c('code_count','total')] %>% as.numeric()
    ref=omic[omic$code_name==i,c('code_count','total')] %>% as.numeric()
    p.value=fisher.test(matrix(c(test[1],ref[1],test[2],ref[2]),nrow=2,ncol=2),alternative="greater")$p.value ## 3extract P value from fischer test
    fischer_p=c(fischer_p,p.value)
    
    names(test)=c('test_N','test_total')
    test.list[[i]]=as.data.frame(t(test))
    
    names(ref)=c('ref_N','ref_total')
    ref.list[[i]]=as.data.frame(t(ref))
    
  }
  
  test.add=rbindlist(test.list)
  ref.add=rbindlist(ref.list)
  trial$fischP=fischer_p 
  trial$fdr=p.adjust(trial$fischP,method='fdr')
  combined=cbind(test.add,ref.add,trial)
  enriched= combined %>% arrange(fdr) %>% filter(fdr<1e-02) %>% filter(test_N>10) %>% select(-frequency,-code_count) %>% 
    mutate(test_freq=test_N/test_total, ref_freq=ref_N/ref_total) %>% 
    select(code_name,fischP,fdr,everything()) %>% data.table()
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
pfam.key=fread('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/annotation/Pfam_key.tsv',header=F,col.names = c('code_name','annotation'))
panth.key=fread('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/annotation/Panther_clean_key.csv')

#### FILTER DATA
#panther=ipscan.filter[database=='PANTHER'] %>% unique.data.frame()
pfam=ipscan.filter[database=='Pfam'] %>% unique.data.frame()


###CALCULATE GENOME WIDE FREQUENCIES
#panth.freq=panther %>% freq.table() %>% data.table()
pfam.freq=pfam %>% freq.table() %>% data.table()



#### GROUP ANNOTATION

gut.specific.genes='./DE_expression_tables/DE_outputs/Common_gut_specific_genes.txt'
#comp.de.files=list.files('./DE_expression_tables/DE_outputs/',full.names = T)[grepl("specific.txt$",list.files('./DE_expression_tables/DE_outputs/'))]
comp.fuzz.files=list.files('./Mfuzz',full.names = T)[grepl("txt$",list.files('./Mfuzz'))]
top300.files=list.files('./Top300',full.names = T)[grepl("txt$",list.files('./Top300'))]

files=c(gut.specific.genes,comp.fuzz.files,top300.files)
short=c()
for(i in files){short=c(short,strsplit(i,'\\/')[[1]][length(strsplit(i,'\\/')[[1]])])}
comparisons.list=list()
for(i in files){comparisons.list[[i]]=readLines(i)}


### Get annotations
annotations=lapply(comparisons.list,function(x) clean.annot[Nv_gene %in% x])
for(i in short){fwrite(annotations[[grep(i,names(annotations))]],paste0('./Bulk_gene_analysis/',i,'_individual_gene_annotation.csv'),row.names = F)}

##make frequency tables
#panth.specific=lapply(comparisons.list,gene.subset,panther)
pfam.specific=lapply(comparisons.list,gene.subset,pfam)


#panth.specific.freq=list()
#for(i in names(panth.specific)){panth.specific.freq[[i]]=freq.table(panth.specific[[i]])}

pfam.specific.freq=list()
for(i in names(pfam.specific)){pfam.specific.freq[[i]]=freq.table(pfam.specific[[i]])}


##caluclate enrichment 
#panth.enrich=lapply(panth.specific.freq,vlepo.enrich,panth.freq)
pfam.enrich=lapply(pfam.specific.freq,vlepo.enrich,pfam.freq)



##anotate PFAM file
annotate_enrichment=function(enrich.table,key){
  l=list()
  for(i in names(enrich.table)){
    enrich=enrich.table[[i]]
    annot=merge(enrich,key,by='code_name',all.y=F,sort=F)
    l[[i]]=(annot[!duplicated(annot)])
  } 
  return(l)
}
pfam.annotate=annotate_enrichment(pfam.enrich,pfam.key)
#panth.annotate=annotate_enrichment(panth.enrich,panth.key)
setwd('./Bulk_gene_analysis/')
for(i in short){fwrite(pfam.annotate[[grep(i,names(pfam.annotate))]],paste0(i,'PFAM_enrichment.csv'))}


for(i in list.files()){
  if(grepl("Cluster_7",i)){
    a=gsub("7","M1_specific",i)
    file.rename(i,a)

