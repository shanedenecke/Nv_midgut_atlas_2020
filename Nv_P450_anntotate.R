### annotate P450s
library(dplyr)
library(data.table)
library(tidyr)


setwd('/data2/shane/Documents/Nezara_midgut_atlas/Nv_P450/annotate')
blast.output=fread('./Nv_P450_blast_output.tsv',col.names = c('Nv_gene','annot','pident','evalue','querycov'))


blast.output$annot=substr(blast.output$annot,1,20)


l=list()
for(i in blast.output$Nv_gene){
  sub=subset(blast.output,Nv_gene==i)[1]
  if(sub$pident>55){
    l[[i]]=data.table(Nezara=sub$Nv_gene,annotation=sub$annot,pident=sub$pident,type='SUBFAMILY LEVEL')
  }else if(sub$pident>35){
    l[[i]]=data.table(Nezara=sub$Nv_gene,annotation=sub$annot,pident=sub$pident,type='FAMILY LEVEL')
  }else{
    l[[i]]=data.table(Nezara=sub$Nv_gene,annotation=sub$annot,pident=sub$pident,type='NEW FAMILY')
  }
}
database.annot=rbindlist(l)
fwrite(rbindlist(l),'Annotate_output.csv',row.names = F)

  
self=fread('./SELF_Nv_P450_blast_output.tsv',col.names = c('Nv_gene','annot','pident','evalue','querycov'))[pident!=100]
colnames(self)=c('Nv_gene','self_annot','Self_pident','self_evalue','self_querycov')
colnames(database.annot)[1]='Nv_gene'
a=merge(database.annot,self,by='Nv_gene')

fwrite(a,'Annotate_SELF_output.csv',row.names = F)
