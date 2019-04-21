
## copy relevant files to new directory
rm -R ~/Dropbox/WP2_omics/Helicoverpa_spatial_gut_atlas/raw_de_tables
mkdir ~/Dropbox/WP2_omics/Helicoverpa_spatial_gut_atlas/raw_de_tables
ls /home/shanedenecke/Documents/omics_data/Helicoverpa_armigera/expression/02_analyze_DE/* | grep "L5" | grep 'UP.subset' | grep -Ev "L4|L3|L2|L1|Mp|Hp|Fp" > ~/temp_files.txt
while read i
do
cp $i ~/Dropbox/WP2_omics/Helicoverpa_spatial_gut_atlas/raw_de_tables
done < ~/temp_files.txt

rm ~/temp_files.txt


## rename files

for i in ~/Dropbox/WP2_omics/Helicoverpa_spatial_gut_atlas/raw_de_tables/*
do
new_name=$(echo $i | perl -pe 's/counts.+(L.+)_vs_(L.+).edgeR.+(L.+UP).+$/$1__$2__$3.csv/g' | perl -pe 's/L([0-9])_C/C$1/g' | sed 's/L5_/L5/g')
mv $i $new_name
done
