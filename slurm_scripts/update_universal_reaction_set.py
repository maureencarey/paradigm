import cobra
import os
import pandas as pd
from cobra.core import Gene, Metabolite, Reaction
import requests

#data_path = "/home/mac9jc/paradigm/data/"
#model_path = "/home/mac9jc/paradigm/models/"
data_path = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data"
model_path = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/models"

def met_ids_without_comp(met_id):
    # only one id listed
    # print list of metabolites without the compartment associated
    # this needs to be updated if you have different compartments than those listed below
    m = met_id
    if m.endswith('_c') or m.endswith('_e') or m.endswith('_f') or \
    m.endswith('_g') or m.endswith('_h') or m.endswith('_i') or \
    m.endswith('_l') or m.endswith('_m') or m.endswith('_n') or \
    m.endswith('_p') or m.endswith('_r') or m.endswith('_s') or \
    m.endswith('_u') or m.endswith('_v') or m.endswith('_x'):
        id_withou_c = m[:-2]
    elif m.endswith('_cx') or m.endswith('_um') or m.endswith('_im') \
    or m.endswith('_ap') or m.endswith('_fv') or m.endswith('_cm'):
        id_withou_c = m[:-3]
    else:
        print('unknown compartment')
        print(m)
        id_withou_c = ''
    
    return(id_withou_c)

def add_sbo_terms(model):
    # Add SBO terms to objects

    for met in model.metabolites:
        met.annotation['sbo'] = ['SBO:0000247']
    for gene in model.genes:
        gene.annotation['sbo'] = ['SBO:0000243']
    for rxn in model.reactions:
        if 'Biomass' in rxn.id or 'biomass' in rxn.id:
            rxn.annotation['sbo'] = ['SBO:0000629']
        elif rxn.id.startswith('EX_'):
            rxn.annotation['sbo'] = ['SBO:0000627']
        elif rxn.id.startswith('DM_'):
            rxn.annotation['sbo'] = ['SBO:0000628']
        elif rxn.id.startswith('SK_'):
            rxn.annotation['sbo'] = ['SBO:0000632']
        elif [met_ids_without_comp(met.id) for met in rxn.reactants] == [met_ids_without_comp(met.id) for met in rxn.products]:
            rxn.annotation['sbo'] = ['SBO:0000185']
        else: rxn.annotation['sbo'] = ['SBO:0000176']

    return(model)

def add_full_met_info(model, met, met_id):
    # add full met info from BiGG api
    
    id_map = {'EC Number':'ec-code','RHEA':'rhea', 'KEGG Reaction':'kegg.reaction',
            'KEGG Compound':'kegg.compound', 'SEED Reaction':'seed.reaction',
            'SEED Compound':'seed.compound', 'MetaNetX (MNX) Equation':'metanetx.reaction',
            'MetaNetX (MNX) Chemical':'metanetx.chemical', 'PubChem':'pubchem.compound',
            'BioCyc':'biocyc','Reactome':'reactome','Brenda':'brenda','LipidMaps':'lipidmaps',
            'Human Metabolome Database':'hmdb','CHEBI':'chebi','InChI':'inchikey'}

    m = requests.get('http://bigg.ucsd.edu/api/v2/universal/metabolites/{}'.format(met_id))
    x = m.json()
    
    if met.name == '': met.name = x['name']
        
    # fix current reaction.annotation field that is in list form
    # this fix will allow universal to be written as as a xml file, not just json
    # however save the info just in case it is not duplicated elsewhere

    # move inappropriately formated annotations object into notes field
    annot_list = met.annotation
    if isinstance(annot_list, list):
        if len(annot_list) >1:
            for sub_list in annot_list:
                if sub_list[0] not in met.notes.keys():
                    met.notes[sub_list[0]] = sub_list[1]
                else: met.notes[sub_list[0]] = [met.notes[sub_list[0]]].append(sub_list[1])
    
    list_o_problem_types = list() # database links is not a list
    list_o_problem_types_2 = list() # database links doesn't have IDs
    list_o_problem_mets = list() # no database links

    # get annotation for requests object
    if 'database_links' in x.keys():
        temp_annotation = x['database_links']
    else: 
        temp_annotation = dict()
        list_o_problem_mets.append(met.id)

    #if met.id in universal_model.metabolites:
    met.annotation = dict()
    met.annotation['bigg.metabolite'] = [met.id]

    for key in id_map.keys():
        if key in temp_annotation.keys():
            if isinstance(temp_annotation[key],list):
                list_o_ids = list()
                for item_in_list in temp_annotation[key]:
                    if isinstance(item_in_list,dict):
                        if 'id' in item_in_list.keys():
                            list_o_ids.append(item_in_list['id'])
                    else: list_o_problem_types_2.append({met.id:key})
                met.annotation[id_map[key]] = list(set(list_o_ids))
            else: list_o_problem_types.append({met.id:key})

    if x['formulae'] != [] and len(x['formulae'])>0:
        met.formula = x['formulae'][0]
    if x['charges'] != [] and len(x['charges'])>0:
        met.charge = x['charges'][0]
    
    if len(list_o_problem_types)>0:
        print(met.id, ' database links is not a list')
    if len(list_o_problem_types_2)>0:
       	print(met.id, ' database links are formated incorrectly (no id)')
    if len(list_o_problem_mets)>0:
       	print(met.id, ' has no database links')

#data_path = "/home/mac9jc/paradigm/data"
    data_path = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data"
    os.chdir(data_path)
    df = pd.read_table('metanetx_chem_prop.tsv', sep='\t', comment='#')
    if 'metanetx.chemical' in met.annotation.keys():
        id_string = met.annotation['metanetx.chemical'][0]
        if id_string in df['MNX_ID']:
            met.annotation['inchi'] = [str(df.loc[df['MNX_ID'] == id_string]['InChI'].values[0])]
            met.annotation['inchikey'] = [str(df.loc[df['MNX_ID'] == id_string]['InChIKey'].values[0])]

    return(model)

def add_full_rxn_info(model,rxn, rxn_id):
    # add full reaction information from BiGG api
    
    id_map = {'EC Number':'ec-code','RHEA':'rhea', 'KEGG Reaction':'kegg.reaction',
            'KEGG Compound':'kegg.compound','SEED Reaction':'seed.reaction',
            'SEED Compound':'seed.compound','MetaNetX (MNX) Equation':'metanetx.reaction',
            'MetaNetX (MNX) Chemical':'metanetx.chemical','PubChem':'pubchem.compound',
            'BioCyc':'biocyc','Reactome':'reactome','Brenda':'brenda','LipidMaps':'lipidmaps',
            'Human Metabolome Database':'hmdb','CHEBI':'chebi','InChI':'inchikey'}
    
    m = requests.get('http://bigg.ucsd.edu/api/v2/universal/reactions/{}'.format(rxn_id))
    x = m.json()
    
    if rxn.name == '': rxn.name = x['name']
    
    list_o_problem_types = list() # database links is not a list
    list_o_problem_types_2 = list() # database links doesn't have IDs
    list_o_problem_rxns = list() # no database links

    # fix current reaction.annotation field that is in list form
    # this fix will allow universal to be written as as a xml file, not just json
    # however save the info just in case it is not duplicated elsewhere
    # rxn.notes are currently {'original_bigg_id':[id_string]}
    annot_list = rxn.annotation
    if isinstance(annot_list, list):
        if len(annot_list) >1:
            for sub_list in annot_list:
                if sub_list[0] not in rxn.notes.keys():
                    rxn.notes[sub_list[0]] = sub_list[1]
                else:
                    rxn.notes[sub_list[0]] = [rxn.notes[sub_list[0]]].append(sub_list[1])

    if rxn.reaction == '':
        rxn.reaction = x['reaction_string']

    rxn.annotation = dict()

    rxn.annotation['bigg.reaction'] = [rxn.id]

    if 'database_links' in x.keys():
        temp_annotation = x['database_links']
    else: 
        temp_annotation = dict()
        list_o_problem_rxns.append(rxn.id)

    for key in id_map.keys():
        if key in temp_annotation.keys():
            if isinstance(temp_annotation[key],list):
                list_o_ids = list()
                for item_in_list in temp_annotation[key]:
                    if isinstance(item_in_list,dict):
                        if 'id' in item_in_list.keys():
                            list_o_ids.append(item_in_list['id'])
                    else: list_o_problem_types_2.append({rxn.id:key})
                rxn.annotation[id_map[key]] = list(set(list_o_ids))
            else: list_o_problem_types.append({met.id:key})
    # also have info on x['pseudoreaction']

    if len(list_o_problem_types)>0:
        print(rxn.id, ' database links is not a list')
    if len(list_o_problem_types_2)>0:
        print(rxn.id, ' database links are formated incorrectly (no id)')
    if len(list_o_problem_rxns)>0:
        print(rxn.id, ' has no database links')

    return(model)

def fix_charge_or_formula(model):
    # fix object type for metabolite charges and formulas
    
    for met in model.metabolites:
        if isinstance(met.charge, list):
            if len(met.charge) == 0:
                met.charge = int(0)
            else:
                met.charge = met.charge[0]
        if isinstance(met.formula, list):
            if len(met.formula) == 0:
                met.formula = ''
            else:
                met.formula = met.formula[0]

    return(model)

# load universal model
os.chdir(model_path)
#universal = cobra.io.load_json_model('universal_model_oct26_2018.json')
#cobra.manipulation.modify.escape_ID(universal)
#
## remove Biomass reactions
#rxn_list_to_delete = [r.id for r in universal.reactions if r.id.startswith('BIOMASS_')]
#universal.remove_reactions(rxn_list_to_delete)
#
## add full met and rxn info
#for met in universal.metabolites:
#    universal = add_full_met_info(universal, met, met_ids_without_comp(met.id))
#for rxn in universal.reactions:
#    universal = add_full_rxn_info(universal, rxn, rxn.id)
#
## fix charges and formulas
#universal = fix_charge_or_formula(universal)
#
## manual curation
#universal.reactions.PGM.lower_bound = -1000. # make reversible
##universal.reactions.ATPM # pseudoreaction for ATP maintenance and sampe as NTP1
##universal.reactions.RPEc # duplicate with RPEc
#universal.remove_reactions([universal.reactions.RPEc,universal.reactions.ATPM])
#
## Add SBO terms
#universal = add_sbo_terms(universal)
#
## SAVE
#os.chdir(model_path)
#cobra.io.save_json_model(universal, 'universal_model_updated.json')
#cobra.io.write_sbml_model(universal, 'universal_model_updated.xml')

universal = cobra.io.load_json_model('universal_model_updated.json')
rxn_list = [r.id for r in universal.reactions]
met_list = [m.id for m in universal.metabolites]

# same for Pf model
os.chdir(model_path)
model = cobra.io.load_json_model('iPfal19.json')
cobra.manipulation.modify.escape_ID(model)

# add full met info
for met in model.metabolites:
    met_id = met_ids_without_comp(met.id)
    if met_id+'_c' in met_list:
        model = add_full_met_info(model, met, met_id)

# add full reaction info
for rxn in model.reactions:
    if rxn.id in rxn_list: model = add_full_rxn_info(model, rxn, rxn.id)

# change id so that Memote doesn't think its the biomass reaction
model.reactions.get_by_id('DM_biomass_c').name = 'unblock biomass'
model.reactions.get_by_id('DM_biomass_c').id = 'DM_bm'

# fix charges and formulas
model = fix_charge_or_formula(model)

# fix food vacucole compartement
for met in model.metabolites:
    if met.compartment == 'food vacuole':
        met.compartment = 'food_vacuole'

# TO DO: ask Memote to recognized EuPathDB gene IDs
for gene in model.genes:
    gene.annotation['EuPathDB.genes'] = [gene.id]

# Add SBO terms
model = add_sbo_terms(model)

model.metabolites.get_by_id('5mti_c').charge = int(model.metabolites.get_by_id('5mti_c').charge)

os.chdir(model_path)
model.name = 'iPfal19, curated P. falciparum 3D7'
cobra.io.save_json_model(model,  'iPfal19_updated.json')
cobra.io.write_sbml_model(model,  'iPfal19_updated.xml')

