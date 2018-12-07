cd /Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data
# these database files are too large for github - add to gitignore
for filename in ./*_annotated_Oct2018.fasta; do
  echo "$filename"
  echo "${filename:2:${#filename}-26}"
  diamond blastp -d ./bigg_proteins_updated -q $filename -o "${filename:2:${#filename}-26}_BiGG.tsv"
done
mkdir ./diamond_output_BiGG
mv ./*_BiGG.tsv ./diamond_output_BiGG

for filename in ./*_annotated_Oct2018.fasta; do
echo "$filename"
echo "${filename:2:${#filename}-26}"
diamond blastp -d ./orthoMCL_proteins -q $filename -o "${filename:2:${#filename}-26}_orthoMCL.tsv"
done
mkdir ./diamond_output_orthoMCL
mv ./*_orthoMCL.tsv ./diamond_output_orthoMCL
mkdir ./genomes
mv ./*_annotated_Oct2018.fasta ./genomes
