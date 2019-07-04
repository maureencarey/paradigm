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
            or m.id.endswith('_ap') or m.id.endswith('_fv'):
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
             m.id.endswith('_ap') or m.id.endswith('_fv'):
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
