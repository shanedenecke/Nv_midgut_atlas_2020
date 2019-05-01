### Annotate P450s

cd /data2/shane/Documents/Nezara_midgut_atlas/Nv_P450

## makeblast db
#makeblastdb -in ./ref/All_insect_P450_Nelson/All_insect_P450.faa -parse_seqids -dbtype prot


## remove P450s with no conserved domain
/data2/shane/Applications/custom/fasta_remove.py ./validation/draft_outputs/NezVir_unigene.faa_length_ipscan.faa ./ref/Manual_remove.txt > ./annotate/Nv_final_P450.faa

### blast query
blastp -query ./annotate/Nv_final_P450.faa -db ./ref/All_insect_P450_Nelson/All_insect_P450.faa -outfmt "6 qseqid sseqid pident evalue qcovs" -evalue 1e-3 -max_target_seqs 6 -max_hsps 1 -num_threads 10 > ./annotate/Nv_P450_blast_output.tsv

### self blast
blastp -query ./annotate/Nv_final_P450.faa -db /data2/shane/Documents/Nezara_midgut_atlas/raw_omics/nviri_ogs.trinity.cl100.faa -outfmt "6 qseqid sseqid pident evalue qcovs" -evalue 1e-3 -max_target_seqs 2 -max_hsps 1 -num_threads 10 > ./annotate/SELF_Nv_P450_blast_output.tsv




#### align + tree


############## ALIGN AND TREE
b=/data2/shane/Documents/Nezara_midgut_atlas/Nv_P450/align_tree/Hemipteran_dros_P450.faa
mafft --anysymbol --threadtb 24 $b > $b'.aln'
#~/Applications/trimal.v1.2rev59/trimAl/source/trimal -in $b'.aln' -out $b'.aln.trimm'
#~/Applications/custom/fasta_2_phylip.sh $b'.aln.trimm' > $b'.aln.trimm.phy'

/home/pioannidis/Programs/trimAl/source/trimal -in $b'.aln' -out $b'.aln.trimm'
/data2/shane/Applications/custom/fasta_2_phylip.sh $b'.aln.trimm' > $b'.aln.trimm.phy'


raxfile=/data2/shane/Documents/Nezara_midgut_atlas/Nv_P450/align_tree/Hemipteran_dros_P450.faa.aln.trimm.phy
raxdir=/data2/shane/Documents/Nezara_midgut_atlas/Nv_P450/align_tree/
/data2/shane/Applications/raxml/raxmlHPC-PTHREADS-AVX -f a -x 12345 -p 12345 -N 500 -T 32 -m PROTGAMMAAUTO -s $raxfile -n 'Nv_P450.tre' -w $raxdir 



