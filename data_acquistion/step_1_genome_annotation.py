import os
import pandas as pd
import glob

data_path = "/home/mac9jc/paradigm/data"
os.chdir(data_path)

# process for visualization
# read in BiGG annotations
os.chdir(data_path+'/diamond_output_BiGG')
annotations_dict = dict()
columns = ['query_gene', 'BiGG_gene', 'pident', 'length', 'mismatch', 'gapopen','qstart', 'qend', 'sstart', 'send', 'evalue', 'score']
for filename in glob.glob(os.path.join(data_path+'/diamond_output_BiGG', '*_BiGG.tsv')):
    annotations_dict[filename.split('/')[len(filename.split('/'))-1]] = pd.read_table(filename, sep = '\t', names=columns)

# BiGG annotations
os.chdir(data_path+'/diamond_output_BiGG')
annotations_dict = dict()
columns = ['query_gene', 'BiGG_gene', 'pident', 'length', 'mismatch', 'gapopen','qstart', 'qend', 'sstart', 'send', 'evalue', 'score']
for filename in glob.glob(os.path.join(data_path+'/diamond_output_BiGG', '*_BiGG.tsv')):
    annotations_dict[filename.split('/')[len(filename.split('/'))-1]] = pd.read_table(filename, sep = '\t', names=columns)
annotations_df = pd.DataFrame(index = annotations_dict.keys(), columns=['species','BiGG_genes'])
for species, annotations in annotations_dict.items():
    annotations_df.species.loc[species] = species
    annotations_df['BiGG_genes'].loc[species] = set(annotations.BiGG_gene)
os.chdir(data_path)
annotations_df.to_csv('bigg_annotations_per_genome_feb2019.csv')

# count BiGG annotations
annot_count = dict()
annot_count_unique_genes = dict()
for key, value in annotations_dict.items():
    annotations_dict[key] = value[value.evalue < 0.01] # remove low confidence
    annot_count[key] = annotations_dict[key].shape[0] # count annotations found for each organism
    annot_count_unique_genes[key] = len(annotations_dict[key].query_gene.unique()) # count unique gene annotations found for each organism
pd.DataFrame.from_dict(annot_count, orient="index").to_csv("annotation_count_all_species_feb2019.csv")
pd.DataFrame.from_dict(annot_count_unique_genes, orient="index").to_csv("unique_gene_annotation_count_all_species_feb2019.csv")

# ORTHO_MCL annotations
os.chdir(data_path+'/diamond_output_orthoMCL')
annotations_dict_ortho = dict()
columns = ['query_gene', 'OrthoMCL_gene', 'pident', 'length', 'mismatch', 'gapopen','qstart', 'qend', 'sstart', 'send', 'evalue', 'score']
for filename in glob.glob(os.path.join(data_path+'/diamond_output_orthoMCL', '*_orthoMCL.tsv')):
    annotations_dict_ortho[filename.split('/')[len(filename.split('/'))-1]] = pd.read_table(filename, sep = '\t', names=columns)
annotations_ortho_df = pd.DataFrame(index = annotations_dict_ortho.keys(), columns=['species','OrthoMCL_gene'])
for species, annotations in annotations_dict_ortho.items():
    annotations_ortho_df.species.loc[species] = species
    annotations_ortho_df['OrthoMCL_gene'].loc[species] = set(annotations.OrthoMCL_gene)
os.chdir(data_path)
annotations_ortho_df.to_csv('ortho_annotations_per_genome_feb2019.csv')

# use resultant files for plotting
