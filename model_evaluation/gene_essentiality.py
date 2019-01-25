import cobra
import os
import pandas as pd
import numpy as np
import glob
from cobra.core import Gene, Metabolite, Reaction

## ##### get essential
os.chdir("/home/mac9jc/paradigm/models")
essentiality_screen_models = dict()
essentiality_screen_models['TgondiiGT1'] = cobra.io.load_json_model('gf_TgondiiGT1.json')
essentiality_screen_models['TgondiiME49'] = cobra.io.load_json_model('gf_TgondiiME49.json')
essentiality_screen_models['Pfalciparum3D7'] = cobra.io.load_json_model('gf_Pfalciparum3D7.json')
essentiality_screen_models['PbergheiANKA'] = cobra.io.load_json_model('gf_PbergheiANKA.json')
essentiality_screen_models['PvivaxSal1'] = cobra.io.load_json_model('gf_PvivaxSal1.json')
essentiality_screen_models['ChominisTU502_2012'] = cobra.io.load_json_model('gf_ChominisTU502_2012.json')
essentiality_screen_models['CparvumIowaII'] = cobra.io.load_json_model('gf_CparvumIowaII.json')
essentiality_screen_models['PknowlesiH'] = cobra.io.load_json_model('gf_PknowlesiH.json')
essentiality_screen_models['PcynomolgiB'] = cobra.io.load_json_model('gf_PcynomolgiB.json')

essentiality_screen_models['Pfalciparum3D7_without_ortho'] = cobra.io.load_json_model('gf_without_ortho_Pfalciparum3D7.json')
essentiality_screen_models['PbergheiANKA_without_ortho'] = cobra.io.load_json_model('gf_without_ortho_PbergheiANKA.json')
essentiality_screen_models['PvivaxSal1_without_ortho'] = cobra.io.load_json_model('gf_without_ortho_PvivaxSal1.json')
essentiality_screen_models['PknowlesiH_without_ortho'] = cobra.io.load_json_model('gf_without_ortho_PknowlesiH.json')
essentiality_screen_models['PcynomolgiB_without_ortho'] = cobra.io.load_json_model('gf_without_ortho_PcynomolgiB.json')

###### get essential genes
# use essentiality_screen_models

gene_essentiality_screen_results_raw= dict()
gene_essentiality_screen_results_interpreted = dict()

for species, model in essentiality_screen_models.items():
    raw_results = dict()
    interpreted_results = dict()

    print(species+', gene essenitality screen')
    # set objective
    model.objective = "generic_biomass"

    # don't accidentally use other biomass reaction
    if 'biomass' in [rxn.id for rxn in model.reactions]:
        model.reactions.get_by_id('generic_biomass').upper_bound = 1000.
        model.reactions.get_by_id('generic_biomass').lower_bound = 0.
        model.reactions.get_by_id('biomass').upper_bound = 0.
        model.reactions.get_by_id('biomass').lower_bound = 0.

    max_biomass = model.slim_optimize()

    # knockout and record growth
    for gene in model.genes:
        with model as cobra_model:
            cobra.manipulation.delete_model_genes(cobra_model, [gene.id], cumulative_deletions=True)
            f = cobra_model.slim_optimize()
            if f < 0.1*max_biomass:
                interpreted_results[gene.id] = 'lethal'
            else:
                interpreted_results[gene.id] = 'nonlethal'
            raw_results[gene.id] = f/max_biomass
            cobra.manipulation.undelete_model_genes(cobra_model)

    # save
    gene_essentiality_screen_results_raw[species+'_generic_biomass'] = raw_results
    gene_essentiality_screen_results_interpreted[species+'_generic_biomass'] = interpreted_results

    if species.startswith('P'):

        print(species+', gene essenitality screen')
        # set objective
        model.objective = "biomass"

        # don't accidentally use other biomass reaction
        model.reactions.get_by_id('biomass').upper_bound = 1000.
        model.reactions.get_by_id('biomass').lower_bound = 0.
        model.reactions.get_by_id('generic_biomass').upper_bound = 0.
        model.reactions.get_by_id('generic_biomass').lower_bound = 0.

        max_biomass = model.slim_optimize()

        # knockout and record growth
        for gene in model.genes:
            with model as cobra_model:
                cobra.manipulation.delete_model_genes(cobra_model, [gene.id], cumulative_deletions=True)
                f = cobra_model.slim_optimize()
                if f < 0.1*max_biomass:
                    interpreted_results[gene.id] = 'lethal'
                else:
                    interpreted_results[gene.id] = 'nonlethal'
                raw_results[gene.id] = f/max_biomass
                cobra.manipulation.undelete_model_genes(cobra_model)

        # save
        gene_essentiality_screen_results_raw[species+'_species_biomass'] = raw_results
        gene_essentiality_screen_results_interpreted[species+'_species_biomass'] = interpreted_results

iPfal18 = cobra.io.load_json_model('iPfal18.json')
iPfal18.reactions.get_by_id('biomass').upper_bound = 1000.
iPfal18.reactions.get_by_id('biomass').lower_bound = 0.
iPfal18.objective = "biomass"
raw_results = dict()
max_biomass = iPfal18.slim_optimize()
for gene in iPfal18.genes:
    with iPfal18 as cobra_model:
        cobra.manipulation.delete_model_genes(cobra_model, [gene.id], cumulative_deletions=True)
        f = cobra_model.slim_optimize()
        raw_results[gene.id] = f/max_biomass
        cobra.manipulation.undelete_model_genes(cobra_model)

cols = ['gene_id', 'normalized_growth']
pd.DataFrame.from_dict(raw_results, orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_iPfal18.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['TgondiiME49_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_TgondiiME49_generic_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['TgondiiGT1_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_TgondiiGT1_generic_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['Pfalciparum3D7_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_Pfalciparum3D7_generic_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PvivaxSal1_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PvivaxSal1_generic_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PbergheiANKA_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PbergheiANKA_generic_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['ChominisTU502_2012_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_ChominisTU502_2012_generic_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['CparvumIowaII_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_CparvumIowaII_generic_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PknowlesiH_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PknowlesiH_generic_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PcynomolgiB_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PcynomolgiB_generic_biomass.csv")

pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['Pfalciparum3D7_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_Pfalciparum3D7_species_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PvivaxSal1_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PvivaxSal1_species_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PbergheiANKA_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PbergheiANKA_species_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PknowlesiH_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PknowlesiH_species_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PcynomolgiB_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PcynomolgiB_species_biomass.csv")

pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['Pfalciparum3D7_without_ortho_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_Pfalciparum3D7_without_ortho_generic_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PvivaxSal1_without_ortho_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PvivaxSal1_without_ortho_generic_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PbergheiANKA_without_ortho_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PbergheiANKA_without_ortho_generic_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PknowlesiH_without_ortho_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PknowlesiH_without_ortho_generic_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PcynomolgiB_without_ortho_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PcynomolgiB_without_ortho_generic_biomass.csv")

pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['Pfalciparum3D7_without_ortho_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_Pfalciparum3D7_without_ortho_species_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PvivaxSal1_without_ortho_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PvivaxSal1_without_ortho_species_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PbergheiANKA_without_ortho_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PbergheiANKA_without_ortho_species_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PknowlesiH_without_ortho_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PknowlesiH_without_ortho_species_biomass.csv")
pd.DataFrame.from_dict(gene_essentiality_screen_results_raw['PcynomolgiB_without_ortho_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/gene_essentiality_matrix_PcynomolgiB_without_ortho_species_biomass.csv")


