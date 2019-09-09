import pandas as pd
import cobra
import requests
import time
from cobra.core import Gene, Metabolite, Reaction
import os


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

def get_comp(model,met_id):
    
    # get compartment associated with a metabolite(s)
    if met_id in [met.id for met in model.metabolites]:
        for m in [model.metabolites.get_by_id(met_id)]:
            if m.id.endswith('_c') or m.id.endswith('_e') or m.id.endswith('_f') or \
            m.id.endswith('_g') or m.id.endswith('_h') or m.id.endswith('_i') or \
            m.id.endswith('_l') or m.id.endswith('_m') or m.id.endswith('_n') or \
            m.id.endswith('_p') or m.id.endswith('_r') or m.id.endswith('_s') or \
            m.id.endswith('_u') or m.id.endswith('_v') or m.id.endswith('_x'):
                id_withou_c = m.id[-2:]
            elif m.id.endswith('_cx') or m.id.endswith('_um') or m.id.endswith('_im') or \
             m.id.endswith('_ap') or m.id.endswith('_fv') or m.id.endswith('_cm'):
                id_withou_c = m.id[-3:]
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

def check_biomass(model):
    print(model.objective.expression)

def merge_subunits(genes): # From CarveMe
    """ Merge list of protein subunit genes into complex
        Args: genes (pandas.Series): list of genes
        Returns: str: boolean rule
        """
    genes = genes.dropna()
    
    if len(genes) == 0:
        return None
    else:
        protein = ' and '.join(sorted(genes))
        if len(genes) > 1:
            return '(' + protein + ')'
        else:
            return protein

def merge_subunit_scores(scores): # From CarveMe
    """ Merge scores of all genes in a protein complex.
        Calculates the mean score among all subunits.
        Args: scores: individual gene scores
        Returns: float: merged score
        """
    return scores.fillna(0).mean()

def merge_proteins(proteins): # From CarveMe
    """ Merge all isozymes that catalyze a given reaction.
        Automatically removes all isozymes with missing score.
        Args: proteins (pandas.Series): list of proteins
        Returns: str: boolean rule
        """
    proteins = set(proteins.dropna())
    if not proteins:
        return None
    gpr_str = ' or '.join(sorted(proteins))
    if len(proteins) > 1:
        return '(' + gpr_str + ')'
    else:
        return gpr_str

def merge_protein_scores(scores): # From CarveMe
    """ Merge scores of all isozymes that catalyze a given reaction.
        Calculates the maximum score among all isozymes.
        Args: scores (pandas.Series): protein scores
        Returns: float: merged score
        """
    return scores.max(skipna=True)

def reaction_scoring(annotation, gprs, spontaneous_score=0.0, debug_output=None): # From CarveMe
    """ Calculate reaction scores using new eggnog output.
        Args: annotation (pandas.DataFrame): gene annotation results
        gprs (pandas.DataFrame): BiGG GPR rules
        spontaneous_score (float): score to give to spontaneous reactions (default: 0.0)
        Returns: pandas.DataFrame: reaction scores
        """
    
    # filter best match for each gene
    gene2gene = annotation.sort_values(by='score', ascending=False) \
        .groupby('BiGG_gene', as_index=False).apply(lambda x: x.iloc[0])
    # merge with gpr table
    gprs['BiGG_gene'] = gprs.apply(lambda row: '{}.{}'.format(row['model'], row['gene'][2:]), axis=1)
    gene_scores = pd.merge(gene2gene, gprs, how='right')
    # add default scores for spontaneous genes
    spontaneous = {'G_s0001', 'G_S0001', 'G_s_0001', 'G_S_0001', 'G_KPN_SPONT'}
    gene_scores.loc[gene_scores.gene.isin(spontaneous), 'score'] = spontaneous_score
    gene_scores.loc[gene_scores.gene.isin(spontaneous), 'query_gene'] = 'spontaneous'
    # from gene to protein scores
    protein_scores = gene_scores.groupby(['protein', 'reaction', 'model'], as_index=False) \
        .agg({'query_gene': merge_subunits, 'score': merge_subunit_scores})
    protein_scores.rename(columns={'query_gene': 'GPR'}, inplace=True)
    # from protein to reaction scores
    reaction_scores = protein_scores.groupby(['reaction'], as_index=False) \
        .agg({'GPR': merge_proteins, 'score': merge_protein_scores}).dropna()
    return(reaction_scores)


def prune_unused_metabolites2(cobra_model):
    """ USE THIS UNTIL AUG 31 UPDATES ARE INTEGRATED INTO MASTER COBRAPY BRANCH
    Remove metabolites that are not involved in any reactions and
    returns pruned model
    Parameters
    ----------
    cobra_model: class:`~cobra.core.Model.Model` object
        the model to remove unused metabolites from
    Returns
    -------
    output_model: class:`~cobra.core.Model.Model` object
        input model with unused metabolites removed
    inactive_metabolites: list of class:`~cobra.core.reaction.Reaction`
        list of metabolites that were removed
    """
    output_model = cobra_model.copy()
    inactive_metabolites = [m for m in output_model.metabolites if len(m.reactions) == 0]
    output_model.remove_metabolites(inactive_metabolites)
    return output_model, inactive_metabolites

def get_comp(model,met_id):

    # get compartment associated with a metabolite(s)
    if met_id in [met.id for met in model.metabolites]:
        for m in [model.metabolites.get_by_id(met_id)]:
            if m.id.endswith('_c') or m.id.endswith('_e') or m.id.endswith('_f') or \
            m.id.endswith('_g') or m.id.endswith('_h') or m.id.endswith('_i') or \
            m.id.endswith('_l') or m.id.endswith('_m') or m.id.endswith('_n') or \
            m.id.endswith('_p') or m.id.endswith('_r') or m.id.endswith('_s') or \
            m.id.endswith('_u') or m.id.endswith('_v') or m.id.endswith('_x'):
                id_withou_c = m.id[-2:]
            elif m.id.endswith('_cx') or m.id.endswith('_um') or m.id.endswith('_im') or \
                m.id.endswith('_ap') or m.id.endswith('_fv') or m.id.endswith('_cm'):
                    id_withou_c = m.id[-3:]
            else:
                print('unknown compartment')
                print(m.id)
                id_withou_c = ''
    else:
        id_withou_c = ''
    return(id_withou_c)

def flatten_mixed_list(list_of_interest):
    new_list = list()
    for x in list_of_interest:
        if isinstance(x,list):
            new_list.extend(x)
        else:
            new_list.append(x)
    return(new_list)

def unaccept_comp_intersection(rxn_list_compartments, acceptable_compartments):
    temp = set(acceptable_compartments)
    unacceptable_comp = [value for value in rxn_list_compartments if value not in temp]
    return(unacceptable_comp)


def update_universal_model(model):

    # open/add exchanges and remove unnecessary biomass functions from universal
    for rxn in model.reactions:
        if rxn.id.startswith('EX_'):
            rxn.lower_bound = -1000.
            rxn.upper_bound = 1000.
    for met in model.metabolites:
        if met.id.endswith('_e'):
            if 'EX_'+met.id not in model.reactions: model.add_boundary(met, type = 'exchange')
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

    old_annot = met.annotation
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
    if isinstance(old_annot,dict):
        for key in old_annot.keys():
            if key not in met.annotation.keys():
                met.annotation[key] = old_annot[key]
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
    
    old_annot = rxn.annotation
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
    
    if isinstance(old_annot,dict):
        for key in old_annot.keys():
            if key not in rxn.annotation.keys():
                rxn.annotation[key] = old_annot[key]
          
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
    old_annot = rxn.annotation
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
    if isinstance(old_annot,dict):
       	for key	in old_annot.keys():
       	    if key not in rxn.annotation.keys():
       	       	rxn.annotation[key] = old_annot[key]

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
    old_annot = met.annotation
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

    if isinstance(old_annot,dict):
       	for key	in old_annot.keys():
       	    if key not in met.annotation.keys():
       	       	met.annotation[key] = old_annot[key]

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
