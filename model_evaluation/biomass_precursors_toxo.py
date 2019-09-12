import cobra
import os
import pandas as pd
import numpy as np
import glob
from cobra.core import Gene, Metabolite, Reaction
from datetime import datetime

day = datetime.now().strftime('%d_%m_%Y')

## ##### get essential reactions
os.chdir("/home/mac9jc/paradigm/models")
model_dict = dict()
model_dict['T. gondii ARI'] = cobra.io.read_sbml_model('gf_no_ortho_TgondiiARI.xml')
model_dict['T. gondii FOU'] = cobra.io.read_sbml_model('gf_no_ortho_TgondiiFOU.xml')
model_dict['T. gondii GAB2-2007-GAL-DOM2'] = cobra.io.read_sbml_model('gf_no_ortho_TgondiiGAB2-2007-GAL-DOM2.xml')
model_dict['T. gondii GT1'] = cobra.io.read_sbml_model('gf_no_ortho_TgondiiGT1.xml')
model_dict['T. gondii MAS'] = cobra.io.read_sbml_model('gf_no_ortho_TgondiiMAS.xml')
model_dict['T. gondii ME49'] = cobra.io.read_sbml_model('gf_no_ortho_TgondiiME49.xml')
model_dict['T. gondii p89'] = cobra.io.read_sbml_model('gf_no_ortho_Tgondiip89.xml')
model_dict['T. gondii RUB'] = cobra.io.read_sbml_model('gf_no_ortho_TgondiiRUB.xml')
model_dict['T. gondii TgCatPRC2'] = cobra.io.read_sbml_model('gf_no_ortho_TgondiiTgCatPRC2.xml')
model_dict['T. gondii VAND'] = cobra.io.read_sbml_model('gf_no_ortho_TgondiiVAND.xml')
model_dict['T. gondii VEG'] = cobra.io.read_sbml_model('gf_no_ortho_TgondiiVEG.xml')

print('models loaded')

met_list = ['spmd_c', 'amet_c',
            'thm_c','thmpp_c',
            'pydx5p_c','pydx_c', 
            'r5p_c',
            'thf_c','5mthf_c',
            'coa_c', 
            'ribflv_c', # cofactors
            'fad_c','nad_c','nadp_c',
            'tag_c','dag_c', 'bm_lipid_c', # generic lipids
            'atp_c','gtp_c','ctp_c','ttp_c', # nucleotides
            'datp_c','dgtp_c','dctp_c','dttp_c', # nucleotides
            'bm_protein_c', # all trna ligated amino acids at once
            'ala__L_c', 'arg__L_c', 'asn__L_c', 'asp__L_c', 'cys__L_c', 'gln__L_c', 
            'glu__L_c', 'gly_c', 'his__L_c', 'ile__L_c',  'leu__L_c', 'lys__L_c', 
            'met__L_c', 'phe__L_c', 'pro__L_c', 'ser__L_c', 'thr__L_c', 'trp__L_c', 
            'tyr__L_c', 'val__L_c'] # every amino acid

# list of all potential essential reaction ids and check if model is broadly functional?
list_o_reactions = list()
for species, model in model_dict.items():
    model.objective = 'generic_biomass'
    print({species+' produces X units of biomass':model.slim_optimize()})
   
    # add sinks in case there are backups generating one precursor
    for met in model.reactions.generic_biomass.reactants:
        if 'DM_'+met.id not in [r.id for r in model.reactions]:
            model.add_boundary(met, type = "demand")
    # remove biomass reactions from consideration
    model.remove_reactions([model.reactions.get_by_id('generic_biomass')])
    if species in ['P. falciparum 3D7','P. berghei ANKA','P. vivax Sal1']:
       model.remove_reactions([model.reactions.get_by_id('biomass')])
    list_o_reactions.append([x.id for x in model.reactions])
    model_dict[species] = model
list_o_reactions = [val for sublist in list_o_reactions for val in sublist]
list_o_reactions = list(set(list_o_reactions))

print('made list of reactions')

# perform knockouts
for met_id in met_list:
    print(met_id)
    matrix_of_results_for_a_met = pd.DataFrame(index = list_o_reactions,columns=model_dict.keys())
    
    for species, model in model_dict.items():
        results_for_species = dict()

        # sequential knockouts and record growth
        cobra_model = model.copy()
        if met_id not in [m.id for m in cobra_model.metabolites]:
            print('ERROR: '+met_id+' is not in '+species)

        else:
            # set objective to demand for metabolite
            if 'DM_'+met_id not in [r.id for r in cobra_model.reactions]:
                cobra_model.add_boundary(cobra_model.metabolites.get_by_id(met_id), type = 'demand')
            cobra_model.objective = 'DM_'+met_id

            # set threshold for interpreting knockouts
            max_biomass = cobra_model.slim_optimize()
            if (max_biomass <= 0.0000000000001):
                print('ERROR - {} doesnt produce: {}'.format(species, met_id))
                continue
            grow_no_grow_threshold = 0.5*max_biomass
            
            # for each demand as objective, do essentiality screen
            cobra_model2 = cobra_model.copy()
            for rxn in cobra_model2.reactions:
                cobra_model.reactions.get_by_id(rxn.id).knock_out()
                f = cobra_model2.slim_optimize()
                if f < grow_no_grow_threshold:
                    results_for_species[rxn.id] = 'lethal'
                else:
                    results_for_species[rxn.id] = 'nonlethal'

        # convert results dictionary into a matrix for PCoA
        for rxn_id in list_o_reactions:
            if rxn_id in [r.id for r in model_dict[species].reactions] and rxn_id in results_for_species.keys():
                matrix_of_results_for_a_met.loc[rxn_id,species] = results_for_species[rxn_id]
            else:
                matrix_of_results_for_a_met.loc[rxn_id,species] = 'not in model'

    # save
    matrix_of_results_for_a_met.to_csv("/home/mac9jc/paradigm/data/results/biomass_precursors/biomass_precursors_essentiality_matrix_toxo_{}_{}.csv".format(met_id,day))
    print('made file for ' +met_id)



