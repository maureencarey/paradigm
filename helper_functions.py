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


def id_bad_compartment_rxns(cobra_model,compartment_list, compartment_options_list):
    """ Output list of reactions that are in relevant compartments or need to be moved to a relevant comp
    Parameters
    ----------
    cobra_model: class:`~cobra.core.Model.Model` object
        the model to evaluate reactions and their compartments
    Returns
    -------
    bad_rxn_list: list of class:`~cobra.core.reaction.Reaction.id`
        list of reaction ids that need to be moved to relevant compartments
    good_rxn_list: list of class:`~cobra.core.reaction.Reaction.id`
        list of reaction ids that are ok
    """

    total_compartments = ["_c","_e","_m","_ap","_fv","_k","_glc","_pm"]
    # cytosol, extracellular, mitochondrdia, apicoplast, food vacuole, kinetoplast, glycosome, pseudomitoc$

    not_compartments = list(set(compartment_options_list) - set(compartment_list))

    good_rxn_list = list()
    bad_rxn_list = list()

    for rxn_object in cobra_model.reactions: # if a reaction does not contain any bad compartments
        rxn_bad_counter = 0
        for x in not_compartments:
            if x in rxn_object.reaction or rxn_object.reaction.endswith(x):
                rxn_bad_counter = rxn_bad_counter + 1
        if rxn_bad_counter == 0:
            good_rxn_list.append(rxn_object.id)
        else:
            bad_rxn_list.append(rxn_object.id)
    return good_rxn_list, bad_rxn_list

def move_bad_rxns(cobra_model,bad_rxn_list,alternative_rxns, compartment_list):
    """ Output list of reactions that are in relevant compartments or need to be moved to a relevant comp
    Parameters
    ----------
    cobra_model: class:`~cobra.core.Model.Model` object
        the model to evaluate reactions and their compartments
    bad_rxn_list: list of class:`~cobra.core.reaction.Reaction.id`
        list of reaction ids that need to be moved to relevant compartments
    alternative_rxns:
    Returns
    -------
    remove_rxn_ids_list: list of class:`~cobra.core.reaction.Reaction.id`
        list of reaction ids that need to be removed, and replaced with an existing reaction in a relevant compartment
    add_reaction_ids_list: list of class:`~cobra.core.reaction.Reaction.id`
        list of reaction ids that need to be added, replacing an existing reaction that was located in an irrelevant compartment
    bad_rxns_keep_rewrite_list: list of class:`~cobra.core.reaction.Reaction.id`
        list of reaction ids that need to be removed, and replaced with a novel reaction in a relevant compartment
    """

    bad_rxns_keep_rewrite_list = list()
    add_reaction_list = list()
    remove_rxn_list = list()
    for rxn_id in bad_rxn_list:
        if len(alternative_rxns[rxn_id]['alternative_reactions']) == 0:
            bad_rxns_keep_rewrite_list.append(rxn_id) # no alternative, keep reaction - will have to change via strings
        else:
            alt_rxns = alternative_rxns[rxn_id]['alternative_reactions']
            keep_og = 0
            for alt_rxn_1, locations in alt_rxns.items():
                keep_alt = 0
                for loc in locations:
                    if loc in compartment_list: keep_alt = keep_alt
                    else: keep_alt = keep_alt + 1
                if keep_alt == 0:
                    keep_og = 1
                    add_reaction_list.append(alt_rxn_1)
                else:
                    keep_og = keep_og
            if keep_og == 0:
                bad_rxns_keep_rewrite_list.append(rxn_id) # no usable alternative - will have to change via strings
            else:
                remove_rxn_list.append(rxn_id)
    add_reaction_list = list(set(add_reaction_list))
    return remove_rxn_list, add_reaction_list, bad_rxns_keep_rewrite_list

def fixing_reaction_compartment(fix_these_reactions_list,cobra_model):

    # need to double check that 'hf.' is necessary for hf.met_ids_without_comp once this fxn is moved to hf
    """ move reactions to cytosol or other appropriate compartment
    Parameters
    ----------
    cobra_model: class:`~cobra.core.Model.Model` object
        the model to modify
    fix_these_reactions: list of class:`~cobra.core.reaction.Reaction`
        list of reaction that need to be moved to relevant compartments
    Returns
    -------
    cobra_model: class:`~cobra.core.Model.Model` object
        the modified model
    error_dict: class: `dict` object
        dictionary of any errors to output into logger (key = rxn.id, value = error)
    reactions_added_list: list of class:`~cobra.core.reaction.Reaction.id`
    transport_for_inappropariate_compartment_list: list of class:`~cobra.core.reaction.Reaction.id`
    """

    reactions_added_list = list()
    transport_for_inappropariate_compartment_list = list()
    error_dict = dict()
    for rxn in fix_these_reactions_list:

        if [met_ids_without_comp(cobra_model,x.id) for x in rxn.reactants] == [met_ids_without_comp(cobra_model,x.id) for x in rxn.products]:
            # remove things like x_p + y_p => x_e + y_e
            transport_for_inappropariate_compartment_list.append(rxn.id)
            new_rxn = list()
        else:

            new_rxn = Reaction()
            met_dict = dict()
            for met in rxn.metabolites:

                if get_comp(cobra_model,met.id) == '_p': # move periplasmic metabolites to extracellular instead of cytosol
                    if met_ids_without_comp(cobra_model,met.id)+'_e' not in [x.id for x in cobra_model.metabolites]:
                        met2 = met.copy()
                        met_id_without_comp = met_ids_without_comp(cobra_model,met.id)
                        met2.id = met_id_without_comp+'_e'
                        met2.compartment = 'extracellular'
                        met_dict[met2] = rxn.metabolites[met]
                        cobra_model.add_metabolites(met2) # []
                    else:
                        met2 = cobra_model.metabolites.get_by_id(met_ids_without_comp(cobra_model,met.id)+'_e')
                        met_dict[met2] = rxn.metabolites[met]
                else: # non periplasmic metabolite
                    if met_ids_without_comp(cobra_model,met.id)+'_c' not in [x.id for x in cobra_model.metabolites]:
                        met2 = met.copy()
                        met_id_without_comp     = met_ids_without_comp(cobra_model,met.id)
                        met2.id = met_id_without_comp+'_c'
                        met2.compartment = 'cytoplasm'
                        met_dict[met2] = rxn.metabolites[met]
                        cobra_model.add_metabolites(met2) # []
                    else:
                        met2 = cobra_model.metabolites.get_by_id(met_ids_without_comp(cobra_model,met.id)+'_c')
                        met_dict[met2] = rxn.metabolites[met]

        # fix reaction variables
        if new_rxn:
            new_rxn.add_metabolites(met_dict)
            new_rxn.name = rxn.name
            new_rxn.id = rxn.id+'c'
            new_rxn.lower_bound = rxn.lower_bound
            new_rxn.upper_bound = rxn.upper_bound
            new_rxn.gene_reaction_rule = rxn.gene_reaction_rule
            new_rxn.notes = rxn.notes
            new_rxn.notes['created for paradigm'] = 'true'
            new_rxn.annotation = rxn.annotation
            cobra_model.add_reactions([new_rxn])
            reactions_added_list.append(new_rxn.id)
        l = len(cobra_model.reactions)
        cobra_model.remove_reactions([rxn])
        cobra_model.repair()
        if len(cobra_model.reactions)>l:
            error_dict[rxn.id] = 'failed to remove a reaction'
    return cobra_model, error_dict, reactions_added_list, transport_for_inappropariate_compartment_list


def flatten(list_o_list):
    # convert list of lists to flat list
    return([item for sublist in list_o_list for item in sublist])

def transcript_to_gene_id(transcript_id):

    x = transcript_id

    if '-t' in x:
        y = x.split('-t')[0]
    else: y = x

    if y.endswith('-RA') or y.endswith('-T1'):
        y2 = y[:-3]
    else: y2 = y

    if y2.endswith('.mRNA') or y2.endswith(':mRNA'):
        y3 = y2[:-5]
    else: y3 = y2

    if y3.startswith('rna_'):
        y4 = y3[4:]
    else: y4 = y3

    if y4.endswith('.1') or y4.endswith('-1') or y4.endswith('.2') or y4.endswith('-2') or y4.endswith('.3'):
        y5 = y4[-2]
    else: y5 = y4

    if '.t' in y5:
        y6 = y5.split('.t')[0]
    else: y6 = y5

    if y6.endswith(':pseudogenic_transcript'):
        y7 = y6.split(':pseudogenic_transcript')[0]
    else: y7 = y6

    if y7.startswith('mRNA1_') or y7.startswith('mRNA2_'):
        y8 = y7[6:]
    else: y8 = y7

    if y8.endswith('-mRNA-1-add'):
        y9 = y8.split('-mRNA-1-add')[0]
    else: y9 = y8

    return(y9)

def get_KEGG_id(gene_id, KEGG_DB):
    # get KEGG ids associated with a gene
    KEGG_id_list = KEGG_DB.loc[KEGG_DB['Gene ID'] == gene_id][['Reaction']] # not yet a list
    return(flatten([x.split() for x in flatten(KEGG_id_list.values.tolist())]))

def transform_universal_to_KEGG(universal_model):
    universal_KEGG_dict = dict()
    for rxn in universal_model.reactions:
        annot = rxn.annotation
        k_list = list()
        if 'kegg.reaction' in annot.keys():
            # sometimes list
            if isinstance(annot['kegg.reaction'], list):
                for x in annot['kegg.reaction']:
                    k_list.append(x) 
            else:
                k_list.append(annot['kegg.reaction'])
        universal_KEGG_dict[rxn.id] = k_list
    return(universal_KEGG_dict)

def get_rxn_from_KEGG(KEGG_string_input, universal_model_dict_input):
    # get reactions in universal model with matching KEGG reaction ID
    # returns a list of reaction IDs
    # universal_model_dict is output of fxn transform_universal_to_KEGG
    
    r_list_list = []
    for rxn_name_keys,KEGG_ids_values in universal_model_dict_input.items():
        for x in KEGG_ids_values:
            if x == KEGG_string_input:
                r_list_list.append(rxn_name_keys)
    return(r_list_list)
            
def prune_protein_to_gene_id(protein, prune_sequence):
    if prune_sequence in protein: #'-t36_1-p1'
        gene_id = protein[:-9]
    else:
        gene_id = protein
    return(gene_id)
