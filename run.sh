fname_full=$( cat "$1" | awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gbff"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}')
wget ${fname_full}.gz
gunzip -f $(basename "${fname_full}.gz")
fname=$(basename $fname_full .gbff)
GB2faa.pl $fname.gbff > $fname.faa
mkdir -p blast
mkdir -p names
mkdir -p blast_tab
mkdir -p blast_stat
blastp -query $fname.faa -db ami_BGC  -outfmt 6 -evalue 1e-5 -out blast/$fname.out
Rscript look_blast.R $fname.out
rm $fname.gbff 

