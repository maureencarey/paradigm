import pandas as pd
import os

data_path = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data"
os.chdir(data_path)
mapping = pd.read_table("./orthology_mapping_jun22/pf3d7_plasmodb.txt")

# modify spcies column in mapping dataframe to map to filenames NECESSARY IF READING TXT MAPPING FILE

file_names = [SPECIES_ID]
mapping['[Organism]'] = [x.replace(". ", "") for x in mapping['[Organism]']]
mapping['[Organism]']  = [x.replace("strain", "") for x in mapping['[Organism]'] ]
mapping['[Organism]']  = [x.replace(" ", "") for x in mapping['[Organism]'] ]
mapping['[Organism]']  = [x.replace("Pfragilenilgiri", "PfragileNilgiri") for x in mapping['[Organism]'] ]
mapping['[Organism]']  = [x.replace("PknowlesiMalayanStrainPk1A", "PknowlesiMalayanPk1A") for x in mapping['[Organism]'] ]
mapping['[Organism]']  = [x.replace("PvivaxSal-1", "PvivaxSal1") for x in mapping['[Organism]'] ]

mapping.to_csv('plasmodium_gene_mapping.csv'.format(SPECIES_ID), header=True, index=True)
