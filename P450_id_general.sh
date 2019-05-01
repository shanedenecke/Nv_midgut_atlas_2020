#### Nv_P450id
cd /data2/shane/Documents/Nezara_midgut_atlas/Nv_P450/validation
mkdir P450_search
mkdir draft_outputs
for i in $(ls ./proteomes/* | grep -E "*unigene.faa$")
do
  #### Linearize fasta
  /data2/shane/Applications/custom/linearize.fasta.sh $i > temp.faa
  mv temp.faa $i
  
  
  ### get basename and make directory
  a=$(basename $i)
  mkdir ./P450_search/$a
  cd ./P450_search/$a
  
  
  ############### BLAST
  /data2/shane/Applications/custom/recip_blast.sh /data2/shane/Documents/Nezara_midgut_atlas/Nv_P450/ref/curated_P450s/P450_homepage_combined.faa ../../proteomes/$a 1 10 4
  
  ########################################## Length and Pfam filters
  mkdir filter
  ################ Lenght filter
  /data2/shane/Applications/custom/fasta_seq_len_filter.py ./blast/recip_blast.faa 350 650 > ./filter/$a'_P450_length.faa'
  echo "Numer of genes After length filter: ______________" $(grep ">" ./filter/$a'_P450_length.faa' | wc -l) >> summary.txt
  
  ###############Xref with IPscan
  /home/pioannidis/Programs/interproscan-5.30-69.0/interproscan.sh -appl pfam -i ./filter/$a'_P450_length.faa' -o ./filter/$a'_ipscan.tsv' -f TSV
  grep "PF00067" ./filter/$a'_ipscan.tsv' | cut -f 1 | sort -u | /data2/shane/Applications/custom/unigene_fa_sub.sh ./filter/$a'_P450_length.faa' - > ./filter/$a'_length_ipscan.faa'
  echo "Number of genes with After PFAM filter:______________" $(grep ">" ./filter/$a'_length_ipscan.faa' | wc -l) >> summary.txt
  
  ############## ALIGN AND TREE
  ### manual remove
  mafft --threadtb 2 ./filter/$a'_length_ipscan.faa' > ./filter/$a'_length_ipscan.faa.aln'
  cp ./filter/$a'_length_ipscan.faa' /data2/shane/Documents/Nezara_midgut_atlas/Nv_P450/validation/draft_outputs/$a'_length_ipscan.faa'
  cp ./filter/$a'_length_ipscan.faa.aln' /data2/shane/Documents/Nezara_midgut_atlas/Nv_P450/validation/draft_outputs/$a'_length_ipscan.faa.aln'
 
  cd /data2/shane/Documents/Nezara_midgut_atlas/Nv_P450/validation

done

for i in ./draft_outputs/*.aln; do echo $i; grep ">" $i | wc -l; done