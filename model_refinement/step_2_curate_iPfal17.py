import cobra
import os
from cobra import Model, Reaction, Metabolite
import pandas as pd
import requests
import logging
from datetime import datetime

day = datetime.now().strftime('%d_%m_%Y')
logging.basicConfig(filename='step2_{}.log'.format(day), level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)

og_path = "/home/mac9jc/paradigm"
data_path = "/home/mac9jc/paradigm/data"
model_path = "/home/mac9jc/paradigm/models"

os.chdir(model_path)
pf_model = cobra.io.read_sbml_model("iPfal17.xml")
logging.info('finished loading model')

## adding notes from previous curation
os.chdir(data_path)
edits = pd.read_csv("iPfal17_edits_BMCGenomics_Carey.csv")

for x in edits.index:
    
    rxn_id_string = edits['Unnamed: 0'][x]
    
    if rxn_id_string in pf_model.reactions:
        if isinstance(edits.Subsystem[x],str) or str(float(edits.Subsystem[x])).lower() != 'nan':
            pf_model.reactions.get_by_id(rxn_id_string).notes['SUBSYSTEM'] = edits.Subsystem[x]
        if isinstance(edits['Confidence Score'][x],str) or str(float(edits['Confidence Score'][x])).lower() != 'nan':
            pf_model.reactions.get_by_id(rxn_id_string).notes['CONFIDENCE'] = edits['Confidence Score'][x]
       	if isinstance(edits['EC Number'][x],str) or str(float(edits['EC Number'][x])).lower() != 'nan':            
            pf_model.reactions.get_by_id(rxn_id_string).notes['EC_NUMBER'] = edits['EC Number'][x]
        if isinstance(edits.References[x],str) or str(float(edits.References[x])).lower() != 'nan':
            if isinstance(edits.Notes[x],str) or str(float(edits.Notes[x])).lower() != 'nan':
                pf_model.reactions.get_by_id(rxn_id_string).notes['iPfal17_notes'] = {'REFERENCE': edits.References[x], 'NOTES': edits.Notes[x]}
            else:
                pf_model.reactions.get_by_id(rxn_id_string).notes['iPfal17_notes'] = {'REFERENCE': edits.References[x]}
        elif isinstance(edits.Notes[x],str) or str(float(edits.Notes[x])).lower() != 'nan':
            pf_model.reactions.get_by_id(rxn_id_string).notes['iPfal17_notes'] = {'NOTES': edits.Notes[x]}
    elif '[' in rxn_id_string:
        rxn_id_string_brackets = rxn_id_string.replace('[','_LSQBKT_').replace(']','_RSQBKT_')
        if rxn_id_string_brackets in pf_model.reactions:
            if isinstance(edits.Subsystem[x],str) or str(float(edits.Subsystem[x])).lower() != 'nan':
                pf_model.reactions.get_by_id(rxn_id_string_brackets).notes['SUBSYSTEM'] = edits.Subsystem[x]
            if isinstance(edits['Confidence Score'][x],str) or str(float(edits['Confidence Score'][x])).lower() != 'nan':
                pf_model.reactions.get_by_id(rxn_id_string_brackets).notes['CONFIDENCE'] = edits['Confidence Score'][x]
            if isinstance(edits['EC Number'][x],str) or str(float(edits['EC Number'][x])).lower() != 'nan':
                pf_model.reactions.get_by_id(rxn_id_string_brackets).notes['EC_NUMBER'] = edits['EC Number'][x]
            if isinstance(edits.References[x],str) or str(float(edits.References[x])).lower() != 'nan':
                if isinstance(edits.Notes[x],str) or str(float(edits.Notes[x])).lower() != 'nan':
                    pf_model.reactions.get_by_id(rxn_id_string_brackets).notes['iPfal17_notes'] = {'REFERENCE': edits.References[x], 'NOTES': edits.Notes[x]}
                else:
                    pf_model.reactions.get_by_id(rxn_id_string_brackets).notes['iPfal17_notes'] = {'REFERENCE': edits.References[x]}
            elif isinstance(edits.Notes[x],str) or str(float(edits.Notes[x])).lower() != 'nan':
                pf_model.reactions.get_by_id(rxn_id_string_brackets).notes['iPfal17_notes'] = {'NOTES': edits.Notes[x]}
    elif '.' in rxn_id_string:
        rxn_id_string_period = rxn_id_string.replace('.','_PERIOD_')
        if rxn_id_string_period in pf_model.reactions: 
            if isinstance(edits.Subsystem[x],str) or str(float(edits.Subsystem[x])).lower() != 'nan':
                pf_model.reactions.get_by_id(rxn_id_string_period).notes['SUBSYSTEM'] = edits.Subsystem[x]
            if isinstance(edits['Confidence Score'][x],str) or str(float(edits['Confidence Score'][x])).lower() != 'nan':
                pf_model.reactions.get_by_id(rxn_id_string_period).notes['CONFIDENCE'] = edits['Confidence Score'][x]    
            pf_model.reactions.get_by_id(rxn_id_string_period).notes['EC_NUMBER'] = rxn_id_string
            if isinstance(edits.References[x],str) or str(float(edits.References[x])).lower() != 'nan':
                if isinstance(edits.Notes[x],str) or str(float(edits.Notes[x])).lower() != 'nan':
                    pf_model.reactions.get_by_id(rxn_id_string_period).notes['iPfal17_notes'] = {'REFERENCE': edits.References[x], 'NOTES': edits.Notes[x]}
                else:
                    pf_model.reactions.get_by_id(rxn_id_string_period).notes['iPfal17_notes'] = {'REFERENCE': edits.References[x]}
            elif isinstance(edits.Notes[x],str) or str(float(edits.Notes[x])).lower() != 'nan':
                pf_model.reactions.get_by_id(rxn_id_string_period).notes['iPfal17_notes'] = {'NOTES': edits.Notes[x]}
    
    else:
        logging.info(edits['Unnamed: 0'][x]+'is in the edits file, but not in the model???')

# update gene identifiers to latest EuPathDB/PlasmoDB version
os.chdir(data_path)
aliases = pd.read_csv('Pfalciparum3D7_GeneAliases.csv')
rename_dictionary1 = dict()
for index, row in aliases.iterrows():
    for x in ['name2','name3','name4','Unnamed: 4','Unnamed: 5','Unnamed: 6','Unnamed: 7']:
        if row[x] != 'nan': rename_dictionary1[row[x]] = row['name1']
    
rename_dictionary = dict()
for key, value in rename_dictionary1.items():
    if isinstance(key, str):
        if key in [x.id for x in pf_model.genes]:
            rename_dictionary[key] = value
        elif key.startswith('MAL'):
            s = key.replace('.','_')
            if s in [x.id for x in pf_model.genes]:
                rename_dictionary[s] = value
                
for key, value in rename_dictionary.items():
    if '.' in value:
        rename_dictionary[key] = value.split('.')[0]

def rename_genes_updated(cobra_model, rename_dict):

    from six import iteritems
    from ast import NodeTransformer
    from cobra.core import Gene, Metabolite, Reaction
    from cobra.core.gene import ast2str
    from cobra.manipulation.delete import get_compiled_gene_reaction_rules

    recompute_reactions = set()  # need to recomptue related genes
    remove_genes = []
    for old_name, new_name in iteritems(rename_dict):
        # undefined if there a value matches a different key because dict is unordered
        try:
            gene_index = cobra_model.genes.index(old_name)
        except ValueError:
            gene_index = None
        old_gene_present = gene_index is not None
        new_gene_present = new_name in cobra_model.genes
        if old_gene_present and new_gene_present:
            old_gene = cobra_model.genes.get_by_id(old_name)
            if old_gene != cobra_model.genes.get_by_id(new_name): 
                 # Added in case not renaming some
                remove_genes.append(old_gene)
                recompute_reactions.update(old_gene._reaction)
        elif old_gene_present and not new_gene_present: # rename old gene to new gene
            gene = cobra_model.genes[gene_index]
            # trick DictList into updating index
            cobra_model.genes._dict.pop(gene.id)  # ugh
            gene.id = new_name
            cobra_model.genes[gene_index] = gene
        elif not old_gene_present and new_gene_present: pass # already fixed
        else:  pass # not old gene_present and not new_gene_present
            # the new gene's _model will be set by repair
            # cobra_model.genes.append(Gene(new_name)) 
            # Removed, otherwise, adds genes that are unassigned to reactions
    cobra_model.repair()

    class Renamer(NodeTransformer):
        def visit_Name(self, node):
            node.id = rename_dict.get(node.id, node.id)
            return node

    gene_renamer = Renamer()
    for rxn, rule in iteritems(get_compiled_gene_reaction_rules(cobra_model)):
        if rule is not None:
            rxn._gene_reaction_rule = ast2str(gene_renamer.visit(rule))

    for rxn in recompute_reactions:
        rxn.gene_reaction_rule = rxn._gene_reaction_rule
    for i in remove_genes:
        cobra_model.genes.remove(i)
    
rename_genes_updated(pf_model,rename_dictionary)

logging.info('finished renaming curation')

# this set of curation steps comes from BioRxiv Untaroiu, Carey, Guler and Papin (2018)
pf_model.reactions.get_by_id('HMGLB').add_metabolites(\
{pf_model.metabolites.get_by_id('h2o2_c'):1.,
pf_model.metabolites.get_by_id('h_c'):-2.})
pf_model.reactions.get_by_id('HMGLB').notes['REFERENCES'] = \
'doi: 10.1073/pnas.0601876103; DOI: 10.1111/j.1365-2141.1975.tb00540.x'

pheme_c = pf_model.metabolites.get_by_id('pheme_c')
gthrd_c = pf_model.metabolites.get_by_id('gthrd_c')
gthox_c = pf_model.metabolites.get_by_id('gthox_c')
pheme_fv = pf_model.metabolites.get_by_id('pheme_fv')
h2o2_c = pf_model.metabolites.get_by_id('h2o2_c')
heme_degraded_c = Metabolite('heme_degraded_c', formula='',
    name='degraded heme', compartment='c')
heme_degraded_fv = Metabolite('heme_degraded_fv',formula='',
    name='degraded heme',compartment='fv')

new_rxn = Reaction()
new_rxn.name = 'gthrd_heme'
new_rxn.id = 'gthrd_heme'
new_rxn.add_metabolites({pheme_c : -1,gthrd_c : -1,
    gthox_c : +1,    heme_degraded_c : +1 })
new_rxn.lower_bound = 0.
new_rxn.upper_bound = 1000.
new_rxn.notes['REFERENCES'] = 'doi: 10.1074/jbc.270.42.24876'
pf_model.add_reactions([new_rxn])

new_rxn = Reaction()
new_rxn.name = 'perox_heme'
new_rxn.id = 'perox_heme'
new_rxn.add_metabolites({pheme_fv : -1, h2o2_c : -1, heme_degraded_fv : +1 })
new_rxn.lower_bound = 0.
new_rxn.upper_bound = 1000.
new_rxn.notes['REFERENCES'] = 'doi: 10.1042/bj1740893'
pf_model.add_reactions([new_rxn])

pf_model.add_boundary(heme_degraded_c, type="sink", reaction_id="SK_heme_degraded_c",lb=0, ub=1000.0)
pf_model.add_boundary(heme_degraded_fv, type="sink", reaction_id="SK_heme_degraded_fv",
                     lb=0, ub=1000.0)
logging.info('finished Anas curation')

# Make dictionary to make all metabolite IDs compatible with bigg
# IN FUTURE, EXPAND TO ALL MODELS # for model in [pf_curated, chominis, leish]:
os.chdir(model_path)
universal_model = cobra.io.load_json_model('universal_model_updated.json')

# switch _D_ to __D_ to be BiGG compatible
for met in pf_model.metabolites:
    if '_D_' in met.id: 
        met.id = met.id.replace('_D_','__D_')
    elif '_L_' in met.id: 
        met.id = met.id.replace('_L_','__L_')
       
met_dict = {'3oodcoa_c':'3ohodcoa_c', 'Asn_X_Ser_FSLASH_Thr_c':'Asn_X_Ser_Thr_c',
'citrul_c':'citr__L_c','Lcystin_c':'cysi__L_c',
'Lcystin_e':'cysi__L_e','doldp_L_c':'doldp_c',
'dolmanp_L_c':'dolmanp_c','dhor_S_m':'dhor__S_m',
'dhor_S_c':'dhor__S_c','g3m8mpdol_L_c':'g3m8mpdol_c',
'glc__D_e_c':'glc__D_c', # CHECK THAT ITS C
'hcys_l_e':'hcys__L_e','lgt_S_c':'lgt__S_c',
'm5mpdol_L_c':'m5mpdol_c','m6mpdol_L_c':'m6mpdol_c',
'm7mpdol_L_c':'m7mpdol_c','saccrp_L_c':'saccrp__L_c',
'sertrna_sec__c':'sertrna_sec_c','sphmyln_host_c':'sphmyln_hs_c',
'Asn_X_Ser_FSLASH_Thr_e':'Asn_X_Ser_Thr_e',
'sperm_c':'sprm_c','pyrdat_c':'4pyrdx_c',
'pnto_R_c':'pnto__R_c','pnto_R_e':'pnto__R_e',
'proteinSS_c':'protdt_c','proteinSHSH_c':'protds_c',
'proteinSS_ap':'protdt_ap','proteinSHSH_ap':'protds_ap',
'proteinSS_m':'protdt_m','proteinSHSH_m':'protds_m'}

for x in pf_model.metabolites:
    if x.id in met_dict.keys():
        if x.id != met_dict[x.id]:
            x.id = met_dict[x.id]
           
logging.info('finished met renaming')

# a few manual fixes
pf_model.reactions.get_by_id('hcys_ex').remove_from_model() # duplicate with EX_hcys___L_e
pf_model.metabolites.get_by_id('hcys_e').remove_from_model() # duplicate with hcys__L_e
pf_model.add_boundary(pf_model.metabolites.get_by_id('protein_t_e'), type="exchange")
pf_model.reactions.get_by_id('PNTH').remove_from_model() # no GPR, no longer in BiGG
pf_model.reactions.get_by_id('PNTK2').remove_from_model() # incorrect duplicate of PNTK
pf_model.reactions.get_by_id('superox_ex').remove_from_model() # replace with BiGG rxns below
rxn1 = universal_model.reactions.get_by_id('O2St').copy()
rxn2 = universal_model.reactions.get_by_id('O2Stm').copy()
pf_model.add_reactions([rxn1, rxn2])
pf_model.reactions.get_by_id('adpr_ex').remove_from_model() # replace with BiGG rxn below
pf_model.reactions.get_by_id('adpr_t').remove_from_model() # replace with BiGG rxn below
rxn1 = universal_model.reactions.get_by_id('EX_adprib_e').copy()
rxn2 = universal_model.reactions.get_by_id('ADPRIBt').copy()
pf_model.add_reactions([rxn1, rxn2])

r = pf_model.reactions.get_by_id('CYOR_u6m_mt')
r.id = 'CYOR_q8_m'
r.reaction = universal_model.reactions.get_by_id('CYOR_q8_m').reaction

# unconnected rxns with no GPRs below
pf_model.reactions.get_by_id('EX_ura_LPAREN_e_RPAREN_').remove_from_model() # not needed
pf_model.reactions.get_by_id('EX_uri_LPAREN_e_RPAREN_').remove_from_model() # not needed
pf_model.reactions.get_by_id('EX_xmp_LPAREN_e_RPAREN_').remove_from_model() # not needed
pf_model.reactions.get_by_id('EX_xtsn_LPAREN_e_RPAREN_').remove_from_model() # not needed
pf_model.reactions.get_by_id('NOt').remove_from_model() # not needed
pf_model.reactions.get_by_id('EX_f6p_LPAREN_e_RPAREN_').remove_from_model() # not needed
pf_model.reactions.get_by_id('EX_duri_LPAREN_e_RPAREN_').remove_from_model() # not needed
pf_model.reactions.get_by_id('EX_g6p_LPAREN_e_RPAREN_').remove_from_model() # not needed
pf_model.reactions.get_by_id('EX_cytd_LPAREN_e_RPAREN_').remove_from_model() # not needed
pf_model.reactions.get_by_id('EX_dcyt_LPAREN_e_RPAREN_').remove_from_model() # not needed

# should have small bounds because of inefficiency -
# add info from notes from curation for iPfal17 -
# they did not save in the model for some reason
#pf_model.reactions.MLTHFtap
#pf_model.reactions.MLTHFte3
#pf_model.reactions.MLTHFtmt

# replace some bad practice curation from iPfal17
pf_model.reactions.get_by_id('mthgxl_s').id = 'MGSA'
pf_model.reactions.get_by_id('2_PERIOD_1_PERIOD_1_PERIOD_12').id = 'METMT'
pf_model.reactions.get_by_id('METMT').name = 'Methionine methyltransferase'
pf_model.reactions.get_by_id('uri_gf').id = 'ATPUP'
pf_model.reactions.get_by_id('ATPUP').name = 'ATP:uridine 5-phosphotransferase'
pf_model.reactions.get_by_id('ATPUP').notes['KEGG'] = 'R00964'
pf_model.reactions.get_by_id('ATPUP').notes['EC NUMBER'] = '2.7.1.48 or 2.7.1.213'
pf_model.reactions.get_by_id('ATPUP').notes['RHEA'] = '16828'

pf_model.reactions.get_by_id('1_7_1_1_mt').id = 'NITRm'
pf_model.reactions.get_by_id('NITRm').name = 'Nitrate reductase (NADH), mitochondrial'
pf_model.reactions.get_by_id('1_7_1_3_mt').id = 'NTRIRym'
pf_model.reactions.get_by_id('NTRIRym').name = 'Nitrate reductase (NADPH), mitochondrial'


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

# get rid of metabolites that are never used
pf_model, unused = prune_unused_metabolites2(pf_model)
    
# # pf_model.reactions.PUNP8 # CURATED SOMETHING WRONG
# pf_model.reactions.UP4UH1
# # pf_model.reactions.get_by_id('PYRDAT') # SOMETHING WRONG, genes don't make sense


# for rxn in pf_model.reactions:
#     if rxn.id not in [r.id for r in universal_model.reactions]:
#         logging.info(rxn.id)
#         logging.info(rxn.reaction)
# 
# for met in pf_model.metabolites:
#     if met.id.split('_')[0] not in [m.id.split('_')[0] for m in universal_model.metabolites]:
#         logging.info(met.id)
#         logging.info(met.name)
        
        
os.chdir(data_path)
met_document = pd.read_table('bigg_metabolites.txt')

cobra.manipulation.modify.escape_ID(pf_model)
logging.info('finished escape')

# add metabolite info 
need_info = list()
for met in pf_model.metabolites:
    if met.id in [m.id for m in universal_model.metabolites]:
        if met.id.startswith('protein') or met.id.startswith('lipid'):
            test_string  = 's'
        else:
            met_row = met_document[met_document['bigg_id'] == met.id]
            met_row['universal_bigg_id']
            m = requests.get('http://bigg.ucsd.edu/api/v2/universal/metabolites/{}'.\
            format(met_row['universal_bigg_id'][met_row['universal_bigg_id'].index[0]]))
            x = m.json()
            met.name = x['name']
            met.notes = x['database_links']
            met.formula = x['formulae']
            met.charge = x['charges']
    elif met.id.endswith('_ap') or met.id.endswith('_fv'):
        if met.id[:-3] in [met.id[:-2] for met in universal_model.metabolites]:
            m = requests.get('http://bigg.ucsd.edu/api/v2/universal/metabolites/{}'.\
            format(met.id[:-3]))
            x = m.json()
            met.name = x['name']
            met.notes = x['database_links']
            met.formula = x['formulae']
            met.charge = x['charges']
        else:
            need_info.append(met.id)
    elif met.id[:-2] in [met.id[:-2] for met in universal_model.metabolites]:
        m = requests.get('http://bigg.ucsd.edu/api/v2/universal/metabolites/{}'.\
        format(met.id[:-2]))
        x = m.json()
        met.name = x['name']
        met.notes = x['database_links']
        met.formula = x['formulae']
        met.charge = x['charges']
    else:
        need_info.append(met.id)

# add some new metabolites
change_mets = {'hb_e':{'name':'host hemoglobin','formula':[], 'charge':[]}, 
'hb_c':{'name':'host hemoglobin','formula':[], 'charge':[]},
'hemozoin_fv':{'name':'hemozoin','formula':[], 'charge':[]},
'hemozoin_e':{'name':'hemozoin','formula':[], 'charge':[]}, 
'heme_degraded_fv':{'name':'degraded heme','formula':[], 'charge':[]},
'heme_degraded_c':{'name':'degraded heme','formula':[], 'charge':[]},
'5mti_c':{'name':'5-methyl thioinosine','formula':[], 'charge':[]}, 
'pc_e':{'name': 'Phosphatidylcholine','formula':[], 'charge':[]},
'pc_c':{'name': 'Phosphatidylcholine','formula':[], 'charge':[]},
'pe_e':{'name': 'Phosphatidylethanolamine','formula':[], 'charge':[]}, 
'pe_c':{'name': 'Phosphatidylethanolamine','formula':[], 'charge':[]}, 
'all_pe_c':{'name': 'Phosphatidylethanolamine (total)','formula':[], 'charge':[]},
'all_ps_c':{'name': 'Phosphatidylserine (total)','formula':[], 'charge':[]},
'all_pg_c':{'name': 'Phosphatidylglycerol (total)','formula':[], 'charge':[]},
'all_pc_c':{'name': 'Phosphatidylcholine (total)','formula':[], 'charge':[]},
'all_pi_c':{'name': 'Phosphatidylinositol  (total)','formula':[], 'charge':[]},
'all_apg_c':{'name': 'acyl - phosphatidylglycerol (total)','formula':[], 'charge':[]},
'all_dgl_c':{'name': 'diacyl - phosphatidylglycerol (total)','formula':[], 'charge':[]},
'acgpail_c':{'name': 'N-Acetyl-D-glucosaminylphosphatidylinositol','formula':[], 'charge':[]},
'pail_e':{'name': 'phosphatidylinositol','formula':['C9H16O9PRCO2R2CO2'], 'charge':[-1]}, 
'pail_c':{'name': 'phosphatidylinositol','formula':['C9H16O9PRCO2R2CO2'], 'charge':[-1]},
'gpail_c':{'name': 'D-glucosaminylphosphatidylinositol','formula':['C15H28NO13PRCO2R2CO2'], 'charge':[0]},
'gacpail_c':{'name': 'Glucosaminyl-acylphosphatidylinositol','formula':[], 'charge':[0]},
'mgacpail_c':{'name': 'Mannosyl-glucosaminyl-acylphosphatidylinositiol','formula':[], 'charge':[0]},
'crm_c':{'name': 'ceramide','formula':['C18H36NO2RCO'], 'charge':[0]},
'sphmyln_c':{'name': 'Sphingomyelin (generic)','formula':['C23H48N2O5PRCO'], 'charge':[0]},
'sphmyln_e':{'name': 'Sphingomyelin (generic)','formula':['C23H48N2O5PRCO'], 'charge':[0]},
'lipid':{'name': 'Lipid (total)','formula':[], 'charge':[]},
'dhcrm_c':{'name': 'Dihydroceramide','formula':['C18H38NO2RCO'], 'charge':[0]},
'tag_c':{'name': 'triacylglycerols (total)','formula':[], 'charge':[]},
'dag_c':{'name': 'diacylglycerols (total)','formula':[], 'charge':[]},
'dag_e':{'name': 'diacylglycerols (total)','formula':[], 'charge':[]},
'cdpdag_e':{'name': 'CDP-Diacylglycerol','formula':[], 'charge':[]},
'cdpdag_c':{'name': 'CDP-Diacylglycerol','formula':[], 'charge':[]},
'ptd1ino_e':{'name': 'Phosphatidyl 1D myo inositol','formula':[], 'charge':[]},
'ptd3ino_c':{'name': 'Phosphatidyl 3D myo inositol','formula':[], 'charge':[]},
'ptd1ino_c':{'name': 'Phosphatidyl 1D myo inositol','formula':[], 'charge':[]},
'ptd4ino_c':{'name': 'Phosphatidyl 4D myo inositol','formula':[], 'charge':[]}, 
'ptd145bp_c':{'name': '1 Phosphatidyl D myo inositol 4 5 bisphosphate ','formula':[], 'charge':[]},
'protein_t_e':{'name': 'glutathione-protein conjugate (extracellular) ','formula':[], 'charge':[]},
'protein_t_c':{'name': 'glutathione-protein conjugate (intracellular) ','formula':[], 'charge':[]},
'lpchol_c':{'name': 'Lysophosphatidylcholine','formula':['C8H19NO5PRCO2'], 'charge':[0]}, 
'xolest2_e':{'name': 'Cholesterol ester','formula':['C24H46NO7RCO'], 'charge':[0]},
'xolest2_c':{'name': 'Cholesterol ester','formula':['C24H46NO7RCO'], 'charge':[0]},
'gluside_c':{'name':'D-glucosyl-N-acylsphingosine','formula':[], 'charge':[0]},
'gluside_e':{'name':'D-glucosyl-N-acylsphingosine','formula':[], 'charge':[0]},
'ROH_c':{'name':'hydroxyl group on protein','formula':[], 'charge':[]}, 
'ROH_ap':{'name':'hydroxyl group on protein','formula':[], 'charge':[]}, 
'ROH_m':{'name':'hydroxyl group on protein','formula':[], 'charge':[]}, 
'ROOH_ap':{'name':'hydroperoxy group on protein','formula':[], 'charge':[]},
'ROOH_m':{'name':'hydroperoxy group on protein','formula':[], 'charge':[]},
'ROOH_c':{'name':'hydroperoxy group on protein','formula':[], 'charge':[]}, 
'gthox_protein_e':{'name':'oxidized glutathione and protein','formula':[], 'charge':[]}, 
'gthox_protein_c':{'name':'oxidized glutathione and protein','formula':[], 'charge':[]}, 
'up4u_c':{'name':'P1,P4-Bis(5-uridyl) tetraphosphate','formula':['C18H26N4O23P4'], 'charge':[]},
'psertrna_sec_c':{'name':'O-Phosphoseryl-tRNA(Sec)','formula':[], 'charge':[]}} 

# add compartments to all metabolites
for met in pf_model.metabolites:
    if met.id in change_mets.keys():
        met.name = change_mets[met.id]['name']
        met.formula = change_mets[met.id]['formula']
        met.charge = change_mets[met.id]['charge']
        if met.id in need_info: need_info.remove(met.id)
    if met.id.endswith('_ap'):
        met.compartment = 'apicoplast'
    elif met.id.endswith('_c'):
        met.compartment = 'cytosol'
    elif met.id.endswith('_fv'):
        met.compartment = 'food vacuole'
    elif met.id.endswith('_m'):
        met.compartment = 'mitochondria'
    elif met.id.endswith('_e'):
        met.compartment = 'extracellular'
    else:
        met.compartment = 'other'
logging.info('these mets arent in universal')
logging.info(need_info)
logging.info('-----------------------')

pf_model.reactions.get_by_id('ACCOAL').id = 'ACCOAL2_temp'
pf_model.reactions.get_by_id('ACCOAL2').id = 'ACCOAL'
pf_model.reactions.get_by_id('ACCOAL2_temp').id = 'ACCOAL2'

rxn_dict = {'ATPM':'NTP1',
'GTHOr_c':'GTHOr',
'EX_4HBA':'DM_4hba_c',
'ALCD19_D':'ALCD19y',
'GLUASPtmt':'ASPGLU2m',
'g3pi_t':'G3PIt',
'HDCEAtr':'HDCEAt',
'INSTt2r':'INSTt2',
'OCDCEAtr':'OCDCEAt',
'TTDCAt':'TTDCAtr',
'ACP1':'ACP1_FMN',
'PHEMEtmt':'PHEMEtm',
'LPLIPAL1E180pp':'LPLPS1AGPE180',
'PLIPA2E140pp':'PLA2PE_MYRS_MYRS_c',
'PLIPA2E161pp':'PLA2PE_HDE_HDE_c',
'PLIPA2E180pp':'PLA2PE_STC_STC_c',
'LPLIPAL1E140pp':'LPEH140_sn1_c',
'LPLIPAL1E160pp':'LPEH160_sn1_c',
'LPLIPAL1E181pp':'LPEH1819Z_sn1_c',
'progly_t':'PROGLyt',
'MLTHFt':'MLTHFte3',
'cys_br':'HMR_3996',
'KAT8':'ACACT8r',
'ACONTac':'ACONTa',
'G5SADr':'G5SADs',
'SPODMc':'SPODM',
'XPRT':'XPPT',
'EX_HMFURN':'DM_hmfurn_c',
'ALALtmt':'ALAtmi',
'EX_ribflv1':'EX_ribflv_e',
'HDCAtr':'HDCAt',
'NACUP':'NACt',
'ACCOALm':'ACS2',
'NO3t2rpp':'NO3t2',
'ORNtiDF':'ORNt',
'EX_cholesterol':'EX_chsterol_e',
'CYSGLYexR':'CYSGLYex',
'MI145P6K':'AMITP',
'trdrd_t':'r1441',
'PLIPA2E160pp':'PLA2PE_PALM_PALM_c',
'biomass_s':'DM_biomass_c',
'RPE':'RPEc',
'HEX8':'HEX4',
'DUMPK':'URIDK2r',
'EX_folate4':'EX_4abz_e',
'GLCt1r':'GLCt1',
'EX_nicotinamide1':'EX_ncam_e',
'NO2t2rpp':'NO2t2r',
'PPAt':'PPAtr',
'AKG_MALtmt':'AKGMALtm',
'ETHAt2pp':'ETHAt6',
'acp_exp':'SK_apoACP_c',
'FEROpp':'FEROc',
'LPLIPAL1E161pp':'LPEH1619Z_sn1_c',
'ACONTbc':'ACONTb',
'MDHc':'MDH',
'PGMT_2':'PPM',
'ATPH1e':'ATPH1',
'DHAPtmt':'DHAPtm',
'g3pc_t':'G3PCt',
'GTHRDti':'GTHRDt2',
'ATPH2e':'NDP1',
'OCDCAtr':'OCDCAt',
'SPMS2':'SPRMS',
'EX_nicotinamide2':'NCAMUP',
'EX_folate3':'EX_fol_e',
'AKG_OAAtmt':'OAAAKGtm',
'NADH5':'NADHOR',
'ACCOAtm':'ACCOAtm_1',
'EX_homocysteine':'EX_hcys__L_e',
'ILEt2r':'ILELAT1tc',
'FOLATEt':'r0963',
'citrul_ex':'SK_citr__L_c',
'SUCD2_u6m_mt':'SUCDH_q8_m',
'ACONTbm_mt':'ACN_b_m',
'ACONTam_mt':'ACN_a_m',
'PIOHex_mt':'PIt5m',
'ICDHyr_mt': 'ICDHym',
'ATPADPex_mt':'ATPtm',
'SUCOAS1m_mt':'SUCOAS1m',
'CYOOm_mt':'CYOOm',
'anth_ex':'DM_anth_c',
'EX_glc_e':'EX_glc__D_e',
'sperm_ex': 'DM_sprm_c',
'AKGDHm':'AKGDm'}

# rename reactions if not using BiGG canonical ID
for rxn in pf_model.reactions:
    if rxn.id in rxn_dict.keys():
        rxn.id = rxn_dict[rxn.id]
        
# replace improperly formated exchange and transport reactions
for rxn in pf_model.reactions:
    if '_LPAREN_e_RPAREN_' in rxn.id: 
        rxn.id = rxn.id.replace('_LPAREN_e_RPAREN_','_e')
    if '_LSQBKT_e_RSQBKT_' in rxn.id: 
        rxn.id = rxn.id.replace('_LSQBKT_e_RSQBKT_','_e')
    if '_LPAREN_ap_RPAREN_' in rxn.id: 
        rxn.id = rxn.id.replace('_LPAREN_ap_RPAREN_','_ap')
    if '_LSQBKT_ap_RSQBKT_' in rxn.id: 
        rxn.id = rxn.id.replace('_LSQBKT_ap_RSQBKT_','_ap')
    if '_D_' in rxn.id: 
        rxn.id = rxn.id.replace('_D_','__D_')
    if '_L_' in rxn.id: 
        rxn.id = rxn.id.replace('_L_','__L_')
    if rxn.id.endswith('_ex'): 
        rxn.id = 'EX_'+rxn.reactants[0].id
    if rxn.id.endswith('_mt'): 
        rxn.id = rxn.id.replace('_mt','m')

# rename reactions if not using BiGG canonical ID
for rxn in pf_model.reactions:
    if rxn.id in rxn_dict.keys():
        rxn.id = rxn_dict[rxn.id]
        
for met in pf_model.metabolites:
    if '___' in met.id:
        met.id = met.id.replace('___','__')

for rxn in pf_model.reactions:
    if '___' in rxn.id:
        rxn.id = rxn.id.replace('___','__')

# add exchange for all extracellular mets
for met in pf_model.metabolites:
    if met.id.endswith('_e'):
        if 'EX_'+met.id not in [r.id for r in pf_model.reactions]:
            pf_model.add_boundary(met, type="exchange")
        
logging.info('these mets arent used in more than one reaction and have no GPR:')
for met in pf_model.metabolites:
    if len(met.reactions) == 1:
        for rxn in met.reactions:
            if len(rxn.genes) < 1:
                logging.info(met.id)

logging.info('------------------------------')
        
# these reactions are not going to be in the universal reaction bag, but we wouldnt expect
# them to be anyways (dont print them)
need_info_rxn = dict()
for rxn in pf_model.reactions:
    if rxn.id not in [r.id for r in universal_model.reactions]:
        if rxn.id.startswith('EX_') or rxn.id.endswith('_ap') or rxn.id.endswith('tap'):
            s = rxn.id 
        else:
            need_info_rxn[rxn.id] = rxn.reaction
logging.info('\n')
logging.info('not in universal')
logging.info(need_info_rxn)

# no production of protein_t_e or protein_t_c

pf_model.metabolites.get_by_id('5mti_c').charge = 0.
pf_model.metabolites.get_by_id('5mti_c').formula = 'C11H14N4O4S'
pf_model.metabolites.get_by_id('5mti_c').notes = {'KEGG id' : 'C19787'}

# PRINT DUPLICATE REACTIONS
duplicates = dict()
temp_dict = dict() # get products and reactants for every reaction
for rxn in pf_model.reactions:
    rxn_dict = dict()
    check_rxn_products = rxn.products
    check_rxn_reactants = rxn.reactants
    rxn_dict['reactants'] = [x.id for x in check_rxn_reactants]
    rxn_dict['products'] = [x.id for x in check_rxn_products]
    temp_dict[rxn.id] = rxn_dict

for rxn in universal_model.reactions:
    for key in temp_dict.keys():
        if key != rxn.id:
            if [x.id for x in rxn.reactants] == temp_dict[key]['reactants'] and \
            [x.id for x in rxn.products] == temp_dict[key]['products']:
                if rxn.id not in duplicates.keys():
                    duplicates[rxn.id] = key
                elif duplicates[rxn.id] == key or key in duplicates[rxn.id]:
                    continue
                else:
                     duplicates[rxn.id] = duplicates[rxn.id]+', '+key
logging.info('\n')
logging.info('duplicates reactions in universal')
logging.info(duplicates)
#
os.chdir(model_path)
cobra.io.save_json_model(pf_model, "iPfal19.json")
cobra.io.write_sbml_model(pf_model, "iPfal19.xml")
