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
model_dict['TgondiiGT1'] = cobra.io.read_sbml_model('gf_no_ortho_TgondiiGT1.xml')
model_dict['TgondiiME49'] = cobra.io.read_sbml_model('gf_no_ortho_TgondiiME49.xml')
model_dict['Pfalciparum3D7'] = cobra.io.read_sbml_model('gf_Pfalciparum3D7.xml')
model_dict['PbergheiANKA'] = cobra.io.read_sbml_model('gf_PbergheiANKA.xml')
model_dict['PvivaxSal1'] = cobra.io.read_sbml_model('gf_PvivaxSal1.xml')
model_dict['ChominisTU502_2012'] = cobra.io.read_sbml_model('gf_no_ortho_ChominisTU502_2012.xml')
model_dict['CparvumIowaII'] = cobra.io.read_sbml_model('gf_no_ortho_CparvumIowaII.xml')
model_dict['EhistolyticaHM1IMSS-A'] = cobra.io.read_sbml_model('gf_no_ortho_EhistolyticaHM1IMSS-A.xml')
model_dict['LmajorSD75'] = cobra.io.read_sbml_model('gf_no_ortho_LmajorSD75.xml')
model_dict['TcruziCLBrener'] = cobra.io.read_sbml_model('gf_no_ortho_TcruziCLBrener.xml')
model_dict['TbruceiLister427'] = cobra.io.read_sbml_model('gf_no_ortho_TbruceiLister427.xml')
model_dict['VbrassicaformisCCMP3155'] = cobra.io.read_sbml_model('gf_no_ortho_VbrassicaformisCCMP3155.xml')
model_dict['BbovisT2Bo'] = cobra.io.read_sbml_model('gf_no_ortho_BbovisT2Bo.xml')
model_dict['gf_no_ortho_NfowleriATCC30863'] = cobra.io.read_sbml_model('gf_no_ortho_NfowleriATCC30863.xml')
print('models loaded')

met_list = ['spmd_c','thmpp_c','pydx5p_c', 'r5p_c','amet_c', 'thf_c','coa_c', 'ribflv_c', # cofactors
            'fad_c','nad_c','nadp_c',
            'tag_c','dag_c', 'bm_lipid_c', # generic lipids
            'atp_c','gtp_c','ctp_c','ttp_c', # nucleotides
            #'datp_c','dgtp_c','dctp_c','dttp_c', # nucleotides
            'bm_protein_c']#, # all trna ligated amino acids at once
#'ala__L_c', 'arg__L_c', 'asn__L_c', 'asp__L_c', 'cys__L_c', 'gln__L_c', 'glu__L_c', 'gly_c', 'his__L_c', 'ile__L_c',  'leu__L_c', 'lys__L_c', 'met__L_c', 'phe__L_c', 'pro__L_c', 'ser__L_c', 'thr__L_c', 'trp__L_c', 'tyr__L_c', 'val__L_c']] # every amino acid

# list of all potential essential reaction ids and check if model is broadly functional?
list_o_reactions = list()
for species, model in model_dict.items():
    model.objective = 'generic_biomass'
    if model.slim_optimize() < 0.01: print('ERROR - model doesnt grow:' + species)
    # remove biomass reactions from consideration
    model_dict[species] = model.remove_reactions[model.reactions.get_by_id('generic_biomass')]
    if species in ['Pfalciparum3D7','PbergheiANKA','PvivaxSal1']:
        model_dict[species] = model_dict[species].remove_reactions[model_dict[species].reactions.get_by_id('biomass')]
    list_o_reactions.append([x.id for x in model_dict[species].reactions])
list_o_reactions = [val for sublist in list_o_reactions for val in sublist]
list_o_reactions = list(set(list_o_reactions))

print('made list of reactions')

# perform knockouts
for met_id in met_list:
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
            max_biomass = model.slim_optimize()
            if max_biomass < 0.01:
                print('ERROR - {} doesnt produce: {}'.format(species, met_id))
                continue
            grow_no_grow_threshold = 0.5*max_biomass
            
            # for each demand as objective, do essentiality screen
            with cobra_model:
                for rxn in cobra_model.reactions:
                    cobra_model.reactions.get_by_id(rxn.id).knock_out()
                    f = cobra_model.slim_optimize()
                    if f < grow_no_grow_threshold:
                        results_for_species[rxn.id] = 'lethal'
                    else:
                        results_for_species[rxn.id] = 'nonlethal'
                            
        # convert results dictionary into a matrix for PCoA
        for rxn_id in list_o_reactions:
            if rxn_id in [r.id for r in model.reactions]:
                matrix_of_results_for_a_met.loc[rxn_id,species] = results_for_species[rxn_id]
            else:
                matrix_of_results_for_a_met.loc[rxn_id,species] = 'NA'

    # save
    matrix_of_results_for_a_met.to_csv("/home/mac9jc/paradigm/data/results/biomass_precursors/biomass_precursors_essentiality_matrix_{}_{}.csv".format(met,day))
    print('made file for ' +met_id)



