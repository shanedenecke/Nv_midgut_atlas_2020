library(data.table)
library(dplyr)
library(stringr)
library(VennDiagram)

setwd('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas')
dir.create('./Venn_diagrams')
sum.trans=fread('/home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/raw_omics/transcription_summary_means.csv') %>% 
  select(Nv_gene,M1_mean,M2_mean,M3_mean,M4H_mean,C_mean) 




l=list()
for(i in c('M1_mean','M2_mean','M3_mean','M4H_mean','C_mean')){
  a=sum.trans[,c('Nv_gene',i),with=FALSE]
  colnames(a)=c('Nv_gene','expression')
  table=a %>% filter(expression>1) %>% data.table()
  l[[i]]=table$Nv_gene
}

intra.gut=l[1:4]
venn.diagram(intra.gut,filename='./Venn_diagrams/Intragut_comparison.tiff',main='Midgut Comparison',fill=c('coral','cadetblue','gold','lightgreen'))
  
pan.gut=Reduce(intersect,intra.gut)
gut.car=list(gut=pan.gut,carcass=l[[5]])

venn.diagram(gut.car,filename='./Venn_diagrams/gutcar_comparison.tiff',main='Gut Carcass Comparison',fill=c('blue','grey70'))


venn.diagram(l,filename='./Venn_diagrams/Messy_comparison.tiff',main='All Comparisons',fill=c('blue','red','green','yellow','grey70'),width=5000)


