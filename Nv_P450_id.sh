#### Nv_P450id
#cd /home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/Nv_P450id
cd /data2/shane/Documents/Nezara_midgut_atlas/Nv_P450

############### BLAST
#~/Applications/custom/recip_blast.sh ./ref/curated_P450s/Combined_P450_Manual.faa ../raw_omics/nviri_ogs.trinity.cl100.faa 1e-7 60
/data2/shane/Applications/custom/recip_blast.sh ./ref/curated_P450s/Combined_P450_Manual.faa ../raw_omics/nviri_ogs.trinity.cl100.faa 1e-7 60


########################################## Length and Pfam filters
mkdir filter
################ Lenght filter
#~/Applications/custom/fasta_seq_len_filter.py ./blast/recip_blast.faa 350 650 > ./filter/Nv_P450_length.faa
/data2/shane/Applications/custom/fasta_seq_len_filter.py ./blast/recip_blast.faa 350 650 > ./filter/Nv_P450_length.faa

echo "Numer of genes After length filter: ______________" $(grep ">" ./filter/Nv_P450_length.faa | wc -l)

###############Xref with IPscan

/home/pioannidis/Programs/interproscan-5.30-69.0/interproscan.sh -appl pfam -i ./filter/Nv_P450_length.faa -o ./filter/Nv_P450_ipscan.tsv -f TSV
grep "PF00067" ./filter/Nv_P450_ipscan.tsv | cut -f 1 | sort -u | /data2/shane/Applications/custom/unigene_fa_sub.sh ./filter/Nv_P450_length.faa - > ./filter/Nv_P450_length_ipscan.faa


####USING ALREADY PERFORMED IPSCAN
#grep "PF00067" ../annotation/nviri_ogs.trinity.cl100.faa.GOdesc.tsv | cut -f 1 | sort -u > ./filter/IPscan_P450s.txt
#echo "Number of genes with PFAM:___________" $(wc -l ./filter/IPscan_P450s.txt)
#rm ./filter/Nv_P450_length_ipscan.faa
#while read i
#do
#grep -A 1 $i ./filter/Nv_P450_length.faa >> ./filter/Nv_P450_length_ipscan.faa
#done < ./filter/IPscan_P450s.txt

echo "Number of genes with After PFAM filter:______________" $(grep ">" ./filter/Nv_P450_length_ipscan.faa | wc -l)




############## ALIGN AND TREE
### manual remove
mkdir align_tree
cp ./filter/Nv_P450_length_ipscan.faa ./align_tree/Nv_draft_P450.faa
a=./align_tree/Nv_draft_P450.faa 
mafft --threadtb 24 $a > $a'.aln'

#~/Applications/custom/fasta_remove.py $a ./ref/Manual_remove.txt > ./align_tree/Nv_P450s_manual.faa
/data2/shane/Applications/custom/fasta_remove.py $a ./ref/Manual_remove.txt > ./align_tree/Nv_P450s_manual.faa


cat ./align_tree/Nv_P450s_manual.faa ./ref/curated_P450s/DroMel_AcyPis_P450_Manual.faa > ./align_tree/Nv_reference_combined_P450.faa
### align and 
b=./align_tree/Nv_reference_combined_P450.faa
mafft --anysymbol --threadtb 24 $b > $b'.aln'
#~/Applications/trimal.v1.2rev59/trimAl/source/trimal -in $b'.aln' -out $b'.aln.trimm'
#~/Applications/custom/fasta_2_phylip.sh $b'.aln.trimm' > $b'.aln.trimm.phy'

/home/pioannidis/Programs/trimAl/source/trimal -in $b'.aln' -out $b'.aln.trimm'
/data2/shane/Applications/custom/fasta_2_phylip.sh $b'.aln.trimm' > $b'.aln.trimm.phy'


raxfile=/data2/shane/Documents/Nv_P450id/Nv_P450s_manual.faa.aln.trimm.phy
raxdir=/data2/shane/Documents/Nv_P450id/tree
/data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T 32 -m PROTGAMMAAUTO -s $raxfile -n 'Nv_P450.tre' -w $raxdir 








##### expression cross reference
mkdir expression
grep ">" ./align_tree/Nv_P450s_manual.faa | cut -d' ' -f 1 | sed 's/>//g' > ./expression/Nv_P450s.txt
Rscript ../scripts/Nv_P450_expression.R








##initial
mkdir blast
blastp -query ./ref/curated_P450s/Combined_P450_Manual.faa -db ../raw_omics/nviri_ogs.trinity.cl100.faa -outfmt "6 qseqid sseqid evalue qcovs" -evalue 1e-7 -qcov_hsp_perc 60 > ./blast/initial_blast.tsv
cat ./blast/initial_blast.tsv | cut -f 2 | sort -u > ./blast/initial_blast_names.txt
echo "Numer of genes in initial blast"
cat ./blast/initial_blast.tsv | cut -f 2 | sort -u | wc -l

## extract initial genes
rm ./blast/initial_blast.faa
while read i
do
grep -A 1 $i ../raw_omics/nviri_ogs.trinity.cl100.faa >> ./blast/initial_blast.faa
done < ./blast/initial_blast_names.txt

##recip blast
makeblastdb -in ./curated_P450s/Combined_P450_Manual.faa -parse_seqids -dbtype prot 
blastp -query ./blast/initial_blast.faa -db ./curated_P450s/Combined_P450_Manual.faa -outfmt "6 qseqid sseqid evalue qcovs" -evalue 1e-7 -qcov_hsp_perc 60 > ./blast/recip_blast.csv
cat ./blast/recip_blast.csv | grep "Cyp" | cut -f 1 | sort -u > ./blast/recip_blast_names.txt
echo "Numer of genes in Reciprocal blast"
cat ./blast/recip_blast.csv | grep "Cyp" | cut -f 1 | sort -u | wc -l

##extract recip blast
rm ./blast/recip_blast.faa
while read i
do
grep -A 1 $i ../raw_omics/nviri_ogs.trinity.cl100.faa >> ./blast/recip_blast.faa
done < ./blast/recip_blast_names.txt

########################################## Length and Pfam filters
mkdir filter
################ Lenght filter
~/Applications/custom/fasta_seq_len_filter.py ./blast/recip_blast.faa 300 600 > ./filter/Nv_P450_length.faa
echo "Numer of genes After length filter"
grep ">" ./filter/Nv_P450_length.faa | wc -l
###############Xref with IPscan
grep "PF00067" ../annotation/nviri_ogs.trinity.cl100.faa.GOdesc.tsv | cut -f 1 | sort -u > ./filter/IPscan_P450s.txt

echo "Number of genes with PFAM "
wc -l ./filter/IPscan_P450s.txt

rm ./filter/Nv_P450_length_ipscan.faa
while read i
do
grep -A 1 $i ./filter/Nv_P450_length.faa >> ./filter/Nv_P450_length_ipscan.faa
done < ./filter/IPscan_P450s.txt

echo "Number of genes with After PFAM filter"
grep ">" ./filter/Nv_P450_length_ipscan.faa | wc -l



############## ALIGN AND TREE

### manual remove

mkdir align_tree
cp ./filter/Nv_P450_length_ipscan.faa ./align_tree/Nv_draft_P450.faa
a=./align_tree/Nv_draft_P450.faa 
muscle -in $a -out $a'.aln'

~/Applications/custom/fasta_remove.py $a ./ref/Manual_remove.txt > ./align_tree/Nv_P450s_manual.faa
cat ./align_tree/Nv_P450s_manual.faa ./curated_P450s/AcyPis_P450.faa ./curated_P450s/DroMel_P450.faa > ./align_tree/Nv_reference_combined_P450.faa
### align and 
b=./align_tree/Nv_reference_combined_P450.faa
muscle -in $b -out $b'.aln'
~/Applications/trimal.v1.2rev59/trimAl/source/trimal -in $b'.aln' -out $b'.aln.trimm'
~/Applications/custom/fasta_2_phylip.sh $b'.aln.trimm' > $b'.aln.trimm.phy'


raxfile=/data2/shane/Documents/Nv_P450id/Nv_P450s_manual.faa.aln.trimm.phy
raxdir=/data2/shane/Documents/Nv_P450id/tree
/data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T 32 -m PROTGAMMAAUTO -s $raxfile -n 'Nv_P450.tre' -w $raxdir 


##### expression cross reference



mkdir expression
grep ">" ./align_tree/Nv_P450s_manual.faa | cut -d' ' -f 1 | sed 's/>//g' > ./expression/Nv_P450s.txt
Rscript Nv_P450_expression.R
