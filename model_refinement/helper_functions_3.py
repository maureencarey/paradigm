import cobra
import os
import pandas as pd
from cobra.core import Gene, Metabolite, Reaction
import requests
import time

def met_ids_without_comp(model,met_id):
    
    # print list of metabolites without the compartment associated
    # this needs to be updated if you have different compartments than those listed below
    if met_id in [met.id for met in model.metabolites]:
        for m in [model.metabolites.get_by_id(met_id)]:
            if m.id.endswith('_c') or m.id.endswith('_e') or m.id.endswith('_f') or \
            m.id.endswith('_g') or m.id.endswith('_h') or m.id.endswith('_i') or \
            m.id.endswith('_l') or m.id.endswith('_m') or m.id.endswith('_n') or \
            m.id.endswith('_p') or m.id.endswith('_r') or m.id.endswith('_s') or \
            m.id.endswith('_u') or m.id.endswith('_v') or m.id.endswith('_x'):
                id_withou_c = m.id[:-2]
            elif m.id.endswith('_cx') or m.id.endswith('_um') or m.id.endswith('_im') \
            or m.id.endswith('_ap') or m.id.endswith('_fv') or m.id.endswith('_cm'):
                id_withou_c = m.id[:-3]
            else:
                print('unknown compartment')
                print(m.id)
                id_withou_c = ''
    else:
        id_withou_c = ''
    return(id_withou_c)

def update_universal_model(model):

	# open/add exchanges and remove unnecessary biomass functions from universal
	for rxn in model.reactions:
		if rxn.id.startswith('EX_'):
			rxn.lower_bound = -1000.
			rxn.upper_bound = 1000.
	for met in model.metabolites:
		if met.id.endswith('_e'):
			if 'EX_'+met.id not in model.reactions:
				model.add_boundary(met, type = 'exchange')
	for rxn in [r for r in model.reactions if r.id.lower().startswith('biomass')]:
		rxn.remove_from_model()
		
	return(model)

def add_sbo_terms(model):
    # Add SBO terms to objects
    # all other annotation values are list, but memote will not recognize SBO terms embedded in lists

    for met in model.metabolites:
        met.annotation['sbo'] = 'SBO:0000247'
    for gene in model.genes:
        gene.annotation['sbo'] = 'SBO:0000243'
    for rxn in model.reactions:
        annotations = []
        if 'Biomass' in rxn.id or 'biomass' in rxn.id:
            annotations.append('SBO:0000629')
        elif rxn.id.startswith('EX_'):
            annotations.append('SBO:0000627')
        elif rxn.id.startswith('DM_'):
            annotations.append('SBO:0000628')
        elif rxn.id.startswith('SK_'):
            annotations.append('SBO:0000632')
        elif [met_ids_without_comp(model,met.id) for met in rxn.reactants] == [met_ids_without_comp(model,met.id) for met in rxn.products]:
            annotations.append('SBO:0000185')
        else: annotations.append('SBO:0000176')

        if len(annotations) > 1:
            rxn.annotation['sbo'] = annotations
            print(reaction.id + ' has more than one SBO annotation. Memote will not like this.')
        else:
            rxn.annotation['sbo'] = annotations[0]

    return(model)


def add_partial_met_info(model, met, met_id):
    # add met info from notes field

    id_map = {'EC Number':'ec-code','RHEA':'rhea', 'KEGG Reaction':'kegg.reaction',
            'KEGG Compound':'kegg.compound', 'SEED Reaction':'seed.reaction',
            'SEED Compound':'seed.compound', 'MetaNetX (MNX) Equation':'metanetx.reaction',
            'MetaNetX (MNX) Chemical':'metanetx.chemical', 'PubChem':'pubchem.compound',
            'BioCyc':'biocyc','Reactome':'reactome','Brenda':'brenda','LipidMaps':'lipidmaps',
            'Human Metabolome Database':'hmdb','CHEBI':'chebi','InChI':'inchikey'}

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

    met.annotation = dict()
    met.annotation['bigg.metabolite'] = [met.id]

    if isinstance(met.notes, dict):
        for key, value in met.notes.items():
            list_o_ids = list()
            if isinstance(value,list) and len(value) > 0:
                for list_item in value:
                    if isinstance(list_item,dict) and 'id' in list_item.keys():
                        list_o_ids.append(list_item['id'])
                    else: list_o_ids = [value]

            if key in id_map.keys():
                met.annotation[id_map[key]] = list_o_ids
            else:
                met.annotation[key] = list_o_ids

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


def add_full_rxn_info(model,rxn, rxn_id):
    # add full reaction information from BiGG api

    id_map = {'EC Number':'ec-code','RHEA':'rhea', 'KEGG Reaction':'kegg.reaction',
            'KEGG Compound':'kegg.compound','SEED Reaction':'seed.reaction',
            'SEED Compound':'seed.compound','MetaNetX (MNX) Equation':'metanetx.reaction',
            'MetaNetX (MNX) Chemical':'metanetx.chemical','PubChem':'pubchem.compound',
            'BioCyc':'biocyc','Reactome':'reactome','Brenda':'brenda','LipidMaps':'lipidmaps',
            'Human Metabolome Database':'hmdb','CHEBI':'chebi','InChI':'inchikey'}
    x = dict()
    m = ''
    while m == '':
        try:
            m = requests.get('http://bigg.ucsd.edu/api/v2/universal/reactions/{}'.format(rxn_id))
            x = m.json()
            break
        except:
            time.sleep(1)
            continue

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


def add_partial_rxn_info(model, rxn, rxn_id):
    # add rxn info from notes field

    id_map = {'EC Number':'ec-code','RHEA':'rhea', 'KEGG Reaction':'kegg.reaction',
            'KEGG Compound':'kegg.compound', 'SEED Reaction':'seed.reaction',
            'SEED Compound':'seed.compound', 'MetaNetX (MNX) Equation':'metanetx.reaction',
            'MetaNetX (MNX) Chemical':'metanetx.chemical', 'PubChem':'pubchem.compound',
            'BioCyc':'biocyc','Reactome':'reactome','Brenda':'brenda','LipidMaps':'lipidmaps',
            'Human Metabolome Database':'hmdb','CHEBI':'chebi','InChI':'inchikey'}

    # fix current reaction.annotation field that is in list form
    # this fix will allow universal to be written as as a xml file, not just json
    # however save the info just in case it is not duplicated elsewhere

    # move inappropriately formated annotations object into notes field
    annot_list = rxn.annotation
    if isinstance(annot_list, list):
        if len(annot_list) >1:
            for sub_list in annot_list:
                if sub_list[0] not in rxn.notes.keys():
                    rxn.notes[sub_list[0]] = sub_list[1]
                else: rxn.notes[sub_list[0]] = [rxb.notes[sub_list[0]]].append(sub_list[1])

    rxn.annotation = dict()
    rxn.annotation['bigg.metabolite'] = [rxn.id]

    if isinstance(rxn.notes, dict):
        for key, value in rxn.notes.items():
            list_o_ids = list()
            if isinstance(value,list) and len(value) > 0:
                for list_item in value:
                    if isinstance(list_item,dict) and 'id' in list_item.keys():
                        list_o_ids.append(list_item['id'])
            else: list_o_ids = [value]

            if key in id_map.keys():
                rxn.annotation[id_map[key]] = list_o_ids
            else:
                rxn.annotation[key] = list_o_ids

    return(model)

def add_full_met_info(model, met, met_id):
    # add full met info from BiGG api

    id_map = {'EC Number':'ec-code','RHEA':'rhea', 'KEGG Reaction':'kegg.reaction',
            'KEGG Compound':'kegg.compound', 'SEED Reaction':'seed.reaction',
            'SEED Compound':'seed.compound', 'MetaNetX (MNX) Equation':'metanetx.reaction',
            'MetaNetX (MNX) Chemical':'metanetx.chemical', 'PubChem':'pubchem.compound',
            'BioCyc':'biocyc','Reactome':'reactome','Brenda':'brenda','LipidMaps':'lipidmaps',
            'Human Metabolome Database':'hmdb','CHEBI':'chebi','InChI':'inchikey'}
    x = dict()
    m = ''
    while m == '':
        try:
            m = requests.get('http://bigg.ucsd.edu/api/v2/universal/metabolites/{}'.format(met_id))
            x = m.json()
            break
        except:
            time.sleep(1)
            continue

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

    data_path = "/home/mac9jc/paradigm/data"
    #data_path = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data"
    os.chdir(data_path)
    df = pd.read_table('metanetx_chem_prop.tsv', sep='\t', comment='#')
    if 'metanetx.chemical' in met.annotation.keys():
        id_string = met.annotation['metanetx.chemical'][0]
        if id_string in df['MNX_ID'].tolist():
            met.annotation['inchi'] = [str(df.loc[df['MNX_ID'] == id_string]['InChI'].values[0])]
            met.annotation['inchikey'] = [str(df.loc[df['MNX_ID'] == id_string]['InChIKey'].values[0])]

    return(model)
