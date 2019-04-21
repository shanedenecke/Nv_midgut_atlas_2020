
## copy relevant files to new directory

ls /home/shanedenecke/Documents/omics_data/Nezara_viridula/DE_expression_raw* | grep 'UP.subset' > ~/temp_files.txt
cd /home/shanedenecke/Documents/omics_data/Nezara_viridula/DE_expression_raw
while read i
do
cp $i /home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/DE_expression_tables
done < ~/temp_files.txt
rm ~/temp_files.txt


## rename files
cd /home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/
for i in /home/shanedenecke/Dropbox/WP2_omics/Nezara_midgut_atlas/DE_expression_tables/*
do 
#echo $i
newname=$(echo $i | perl -pe 's/kallisto.gene.counts.+([C|M].+)vs_([M|C].+).edgeR.+C2.(.+UP).+$/$1_$2__$3.csv/g')
mv $i $newname
done


