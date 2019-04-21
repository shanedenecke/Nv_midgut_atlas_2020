library(Mfuzz)
library(data.table)
library(dplyr)
library(stringr)

setwd('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas')
dir.create('./Mfuzz')

trans.raw=fread('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/raw_omics/kallisto.gene.TPM.not_cross_norm')
filtered.genes=readLines('./raw_omics/filtered_gene_list.txt')
pheno=read.table('./raw_omics/transcriptomic_phenotype.csv',row.names=1,header=T,sep=',')
colnames(trans.raw)[1]='Nv_gene'
trans.sum=fread('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/raw_omics/transcription_summary_means.csv')



#### Do not run. Takes a while
trans.sum= trans.raw %>% rowwise() %>% mutate(M1_mean=mean(c(M1_1,M1_2,M1_3,M1_4)),M1_sd=sd(c(M1_1,M1_2,M1_3,M1_4))) %>% 
  mutate(M2_mean=mean(c(M2_1,M2_2,M2_3,M2_4)),M2_sd=sd(c(M2_1,M2_2,M2_3,M2_4))) %>% 
  mutate(M3_mean=mean(c(M3_1,M3_2,M3_3,M3_4)),M3_sd=sd(c(M3_1,M3_2,M3_3,M3_4))) %>%
  mutate(M4H_mean=mean(c(M4H_1,M4H_2,M4H_3,M4H_4)),M4H_sd=sd(c(M4H_1,M4H_2,M4H_3,M4H_4))) %>%
  mutate(C_mean=mean(c(C1,C2,C3,C4)),C_sd=sd(c(C1,C2,C3,C4))) %>% mutate(length_sd=sd(M1_mean,M2_mean,M3_mean,M4_mean)) %>% 
  filter(Nv_gene %in% filtered.genes) %>%
  data.table()



### Make expression dataset
trans.matrix=as.matrix(select(trans.sum,Nv_gene,M1_mean,M2_mean,M3_mean,M4H_mean),rownames = 'Nv_gene')
gut.exp=ExpressionSet(assayData=trans.matrix)
gut.exp.f=filter.std(gut.exp,min.std=0,visu=F)
gut.exp.s <- standardise(gut.exp.f)

## Run Mfuzz
m1 <- mestimate(gut.exp.s)
c.8 <- mfuzz(gut.exp.s,c=8,m=m1)

## create graph
pdf('./Mfuzz/8_Clusetrs_fuzzy.pdf')
mfuzz.plot2(gut.exp.s,cl=c.8,mfrow=c(2,4),time.labels = c('M1','M2','M3','M4'),xlab='Gut Compartment',x11=F)
dev.off()


## divide up clusters
cluster=data.table(best.clust=c.8[[3]],Nv_gene=names(c.8[[3]]))
membership=c.8[[4]] %>% data.table()
for(i in 1:length(colnames(membership))){colnames(membership)[i]=paste0('A',colnames(membership)[i])}
membership=membership %>% rowwise() %>% mutate(maximum=max(A1,A2,A3,A4,A5,A6,A7,A8)) %>% data.table() 
membership=cbind(cluster,membership)
fin.membership=membership %>% filter(maximum>.6) %>% arrange(desc(best.clust)) %>% data.table() 
fin.membership$best.clust %>% table()


## write clusters to .csv
for(i in unique(fin.membership$best)){
  fil=fin.membership %>% filter(best.clust==i)
  #write.csv(fil,paste0('./Mfuzz/Cluster_',as.character(i),'_Gene_Table.csv'),row.names = F)
  writeLines(fil$Nv_gene,paste0('./Mfuzz/Cluster_',as.character(i),'_Gene_List.txt'))
}











############### BIN

pheno.update=new("AnnotatedDataFrame",data=pheno, varMetadata=meta)
#meta=data.frame(labelDescription='Self')
#gut.exp <- ExpressionSet(assayData=trans.matrix,phenoData=pheno.update)
gut.exp=minimalSet

load(url("http://duffel.rail.bio/recount/SRP049355/rse_gene.Rdata"))




l2=split(membership,membership$best)
lapply(l2,
}
l[[1]]


trans.sum=trans.sum %>% rowwise() %>% mutate(length_sd=sd(c(M1_mean,M2_mean,M3_mean,M4H_mean))) %>% data.table()
fwrite(trans.sum,'./raw_omics/transcription_summary_means.csv',row.names = F)

test = new('ExpressionSet', exprs=as.matrix(trans.sum))
pheno.update=new("AnnotatedDataFrame",data=pheno, varMetadata=meta)
c1=cselection(gut.exp.s,m1)

data(yeast)
yeast.r <- filter.NA(yeast, thres=0.25)
tmp <- filter.std(yeast.f,min.std=0)



l=list()
for(i in c(4,8,12,16)){
  cl <- mfuzz(gut.exp.s,c=i,m=m1)
  l[[grep(i,c(4,8,12,16))]]=cl
  pdf(paste0(as.character(i),'_Clusetrs_TEST_fuzzy.pdf'))
  mfuzz.plot2(gut.exp.s,cl=cl,mfrow=c(4,(i/4)),time.labels = c('M1','M2','M3','M4'),xlab='Gut Compartment',x11=F)
  dev.off()
} 

#membership=membership %>% rowwise() %>% mutate(maximum=max(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16)) %>% data.table() 



for(i in l){
  cluster=i[[3]] %>% data.table()
  membership=i[[4]] %>% data.table()
  colnames(membership)=c('M1','M2','M3','M4H')
  membership=membership %>% rowwise() %>% mutate(maximum=max(M1,M2,M3,M4H)) %>% data.table() 
  membership$best=cluster
  fin.membership=membership %>% filter(maximum>.8) %>% data.table() 
  fi
  
