import cobra
import os
from cobra import Model, Reaction, Metabolite
import pandas as pd
import requests
import logging
from datetime import datetime
#import helper_functions_3 as hf3
import sys

sys.path.append(os.path.abspath("/home/mac9jc/paradigm/"))
import helper_functions as hf

log_path = "/home/mac9jc/paradigm/model_generation_logs/"

os.chdir(log_path)
day = datetime.now().strftime('%d_%m_%Y')
logging.basicConfig(filename='step2_{}.log'.format(day), level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)

og_path = "/home/mac9jc/paradigm"
data_path = "/home/mac9jc/paradigm/data"
model_path = "/home/mac9jc/paradigm/models"

os.chdir(model_path)
pf_model = cobra.io.read_sbml_model("iPfal17.xml")
universal_model = cobra.io.read_sbml_model('universal_model_updated.xml')
logging.info('finished loading model')

## adding notes from previous curation
os.chdir(data_path)
edits = pd.read_csv("iPfal17_edits_BMCGenomics_Carey.csv")

for x in edits.index:
    
    rxn_id_string = edits['Unnamed: 0'][x]    
    if rxn_id_string in pf_model.reactions:
        rxn_id_string_use = rxn_id_string
    elif '.' in rxn_id_string and rxn_id_string.replace('.','_PERIOD_') in pf_model.reactions:
        rxn_id_string_use = rxn_id_string.replace('.','_PERIOD_')
    elif '[' in rxn_id_string and rxn_id_string.replace('[','_LSQBKT_').replace(']','_RSQBKT_') in pf_model.reactions:
        rxn_id_string_use = rxn_id_string.replace('[','_LSQBKT_').replace(']','_RSQBKT_')
    else:
        logging.info(edits['Unnamed: 0'][x]+'is in the edits file, but not in the model???')
        continue

    if isinstance(edits.Subsystem[x],str) or str(float(edits.Subsystem[x])).lower() != 'nan':
        pf_model.reactions.get_by_id(rxn_id_string_use).subsystem = edits.Subsystem[x]
    if isinstance(edits['Confidence Score'][x],str) or str(float(edits['Confidence Score'][x])).lower() != 'nan':
        pf_model.reactions.get_by_id(rxn_id_string_use).annotation['CONFIDENCE'] = edits['Confidence Score'][x]
    if isinstance(edits['EC Number'][x],str) or str(float(edits['EC Number'][x])).lower() != 'nan':            
        pf_model.reactions.get_by_id(rxn_id_string_use).annotation['ec-code'] = edits['EC Number'][x]
    if isinstance(edits.References[x],str) or str(float(edits.References[x])).lower() != 'nan':
        pf_model.reactions.get_by_id(rxn_id_string_use).annotation['REFERENCE'] = edits.References[x]
        if isinstance(edits.Notes[x],str) or str(float(edits.Notes[x])).lower() != 'nan':
            pf_model.reactions.get_by_id(rxn_id_string_use).notes = {'REFERENCE': edits.References[x], 'NOTES': edits.Notes[x]}
        else:
            pf_model.reactions.get_by_id(rxn_id_string_use).notes = {'REFERENCE': edits.References[x]}
    elif isinstance(edits.Notes[x],str) or str(float(edits.Notes[x])).lower() != 'nan':
        pf_model.reactions.get_by_id(rxn_id_string_use).notes = {'NOTES':edits.Notes[x]}
    pf_model.reactions.get_by_id(rxn_id_string_use).annotation['CURATION'] = 'DOI: 10.1186/s12864-017-3905-1'

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
pf_model.reactions.get_by_id('HMGLB').annotation['REFERENCE'] = \
'doi: 10.1073/pnas.0601876103; DOI: 10.1111/j.1365-2141.1975.tb00540.x'
pf_model.reactions.get_by_id('HMGLB').annotation['CURATION'] = 'doi: 10.1186/s12859-019-2756-y'

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
    gthox_c : +1, heme_degraded_c : +1 })
new_rxn.lower_bound = 0.
new_rxn.upper_bound = 1000.
new_rxn.annotation['REFERENCE'] = 'doi: 10.1074/jbc.270.42.24876'
new_rxn.annotation['CURATION'] = 'doi: 10.1186/s12859-019-2756-y'
pf_model.add_reactions([new_rxn])

new_rxn = Reaction()
new_rxn.name = 'perox_heme'
new_rxn.id = 'perox_heme'
new_rxn.add_metabolites({pheme_fv : -1, h2o2_c : -1, heme_degraded_fv : +1 })
new_rxn.lower_bound = 0.
new_rxn.upper_bound = 1000.
new_rxn.annotation['REFERENCE'] = 'doi: 10.1042/bj1740893'
new_rxn.annotation['CURATION'] = 'doi: 10.1186/s12859-019-2756-y'
pf_model.add_reactions([new_rxn])

pf_model.add_boundary(heme_degraded_c, type="sink", reaction_id="SK_heme_degraded_c",lb=0, ub=1000.0)
pf_model.add_boundary(heme_degraded_fv, type="sink", reaction_id="SK_heme_degraded_fv",
                     lb=0, ub=1000.0)
logging.info('finished Anas curation')
pf_model.repair()

# Make dictionary to make all metabolite IDs compatible with bigg
# IN FUTURE, EXPAND TO ALL MODELS # for model in [pf_curated, chominis, leish]:
os.chdir(model_path)
universal_model = cobra.io.load_json_model('universal_model.json')
rxn_list = [r.id for r in universal_model.reactions]
met_list = [m.id for m in universal_model.metabolites]

# switch _D_ to __D_ to be BiGG compatible
for met in pf_model.metabolites:
    if '_D_' in met.id:
        met.id = met.id.replace('_D_','__D_')
    elif '_L_' in met.id:
        met.id = met.id.replace('_L_','__L_')

# these reactions are missing pi, h, etc. replace formula with BiGG version
notes = pf_model.reactions.GLUDxi.notes
gpr = pf_model.reactions.GLUDxi.gene_reaction_rule
pf_model.remove_reactions([pf_model.reactions.GLUDxi])
pf_model.add_reactions([universal_model.reactions.GLUDxi.copy()])
pf_model.reactions.GLUDxi.gene_reaction_rule = gpr
pf_model.reactions.GLUDxi.notes = notes
pf_model.repair()

gpr = pf_model.reactions.trdrd_exp.gene_reaction_rule
notes = pf_model.reactions.trdrd_exp.notes
pf_model.remove_reactions([pf_model.reactions.trdrd_exp])
pf_model.add_boundary(pf_model.metabolites.trdrd_c, type = "sink")
pf_model.reactions.SK_trdrd_c.notes = notes
pf_model.reactions.SK_trdrd_c.gene_reaction_rule = gpr
pf_model.reactions.SK_trdrd_c.name = 'thioredoxin expression'
pf_model.repair()

gpr = pf_model.reactions.fldox_exp.gene_reaction_rule
notes = pf_model.reactions.fldox_exp.notes
pf_model.remove_reactions([pf_model.reactions.fldox_exp])
pf_model.add_boundary(pf_model.metabolites.fldox_ap, type = "sink")
pf_model.reactions.SK_fldox_ap.notes = notes
pf_model.reactions.SK_fldox_ap.gene_reaction_rule = gpr
pf_model.reactions.SK_fldox_ap.name = 'flavodoxin expression'
pf_model.repair()

#notes = pf_model.reactions.PPPGOm.notes
#gpr = pf_model.reactions.PPPGOm.gene_reaction_rule
#pf_model.remove_reactions([pf_model.reactions.PPPGOm])
#pf_model.add_reactions([universal_model.reactions.PPPGOm.copy()])
#pf_model.reactions.PPPGOm.gene_reaction_rule = gpr
#pf_model.reactions.PPPGOm.notes = notes

notes = pf_model.reactions.PDX5POi.notes
gpr = pf_model.reactions.PDX5POi.gene_reaction_rule
pf_model.remove_reactions([pf_model.reactions.PDX5POi])
pf_model.add_reactions([universal_model.reactions.PDX5POi.copy()])
pf_model.reactions.PDX5POi.gene_reaction_rule = gpr
pf_model.reactions.PDX5POi.notes = notes
pf_model.repair()

notes = pf_model.reactions.PPA.notes
gpr = pf_model.reactions.PPA.gene_reaction_rule
pf_model.remove_reactions([pf_model.reactions.PPA])
pf_model.add_reactions([universal_model.reactions.PPA])
pf_model.reactions.PPA.gene_protein_rule = gpr
pf_model.reactions.PPA.notes = notes
pf_model.repair()

notes = pf_model.reactions.DOLPMT.notes
gpr = pf_model.reactions.DOLPMT.gene_reaction_rule
pf_model.remove_reactions([pf_model.reactions.DOLPMT])
pf_model.add_reactions([universal_model.reactions.DOLPMT.copy()])
pf_model.reactions.DOLPMT.gene_reaction_rule = gpr
pf_model.reactions.DOLPMT.notes = notes
pf_model.repair()

notes = pf_model.reactions.PIt2r.notes
gpr = pf_model.reactions.PIt2r.gene_reaction_rule
pf_model.remove_reactions([pf_model.reactions.PIt2r])
pf_model.add_reactions([universal_model.reactions.PIt2r.copy()])
pf_model.reactions.PIt2r.gene_reaction_rule = gpr
pf_model.reactions.PIt2r.notes = notes

# remove these reactions but add a related BiGG reaction
pf_model.add_reactions([universal_model.reactions.THD2.copy()])
pf_model.reactions.THD2.gene_reaction_rule = pf_model.reactions.THD2pp.gene_reaction_rule
pf_model.reactions.THD2.notes = pf_model.reactions.THD2pp.notes
pf_model.remove_reactions([pf_model.reactions.THD2pp])
pf_model.repair()

pf_model.add_reactions([universal_model.reactions.DOLPMT1_c.copy()])
pf_model.reactions.DOLPMT1_c.gene_reaction_rule = pf_model.reactions.DOLPMT1.gene_reaction_rule
pf_model.reactions.DOLPMT1_c.notes = pf_model.reactions.DOLPMT1.notes
pf_model.remove_reactions([pf_model.reactions.DOLPMT1])

pf_model.add_reactions([universal_model.reactions.G12MT3_c.copy()])
pf_model.reactions.G12MT3_c.gene_reaction_rule = pf_model.reactions.DOLPMT2.gene_reaction_rule
pf_model.reactions.G12MT3_c.notes = pf_model.reactions.DOLPMT2.notes
pf_model.remove_reactions([pf_model.reactions.DOLPMT2])

pf_model.add_reactions([universal_model.reactions.GLCNACPT_c.copy()])
pf_model.reactions.GLCNACPT_c.gene_reaction_rule = pf_model.reactions.GLCNACPT.gene_reaction_rule
pf_model.reactions.GLCNACPT_c.notes = pf_model.reactions.GLCNACPT.notes
pf_model.remove_reactions([pf_model.reactions.GLCNACPT])

pf_model.remove_reactions([pf_model.reactions.G_Protein_Ex,pf_model.reactions.HMBZ_out])
pf_model.add_boundary(pf_model.metabolites.gthox_protein_e, type = "exchange")
pf_model.reactions.EX_gthox_protein_e.annotation['REFERENCE'] = ['DOI: 10.3390/molecules200610511']
pf_model.add_boundary(pf_model.metabolites.hemozoin_e, type = "exchange")
pf_model.repair()

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
pf_model.repair()

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
pf_model.repair()

notes = pf_model.reactions.MGSA.notes
gpr = pf_model.reactions.MGSA.gene_reaction_rule
pf_model.remove_reactions([pf_model.reactions.MGSA])
pf_model.add_reactions([universal_model.reactions.MGSA.copy()])
pf_model.reactions.MGSA.gene_reaction_rule = gpr
pf_model.reactions.MGSA.notes = notes
pf_model.repair()

pf_model.reactions.get_by_id('2_PERIOD_1_PERIOD_1_PERIOD_12').id = 'METMT'
pf_model.reactions.get_by_id('METMT').name = 'Methionine methyltransferase'
pf_model.reactions.get_by_id('uri_gf').id = 'ATPUP'
pf_model.reactions.get_by_id('ATPUP').name = 'ATP:uridine 5-phosphotransferase'
pf_model.reactions.get_by_id('ATPUP').annotation['kegg.reaction'] = 'R00964'
pf_model.reactions.get_by_id('ATPUP').annotation['ec-code'] = ['2.7.1.48','2.7.1.213']
pf_model.reactions.get_by_id('ATPUP').annotation['rhea'] = '16828'

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
        
os.chdir(data_path)
#met_document = pd.read_table('bigg_metabolites.txt')

cobra.manipulation.modify.escape_ID(pf_model)
logging.info('finished escape')

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
pf_model.repair()

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
'AKGDHm':'AKGDm',
'EX_folate2':'6HMHPTtr'}

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

pf_model.repair()

# add exchange for all extracellular mets
for met in pf_model.metabolites:
    if met.id.endswith('_e'):
        if 'EX_'+met.id not in [r.id for r in pf_model.reactions]:
            pf_model.add_boundary(met, type="exchange")
for rxn in pf_model.reactions:
    if rxn.id.startswith('EX_'):
        pf_model.reactions.get_by_id(rxn.id).lower_bound = -1000.
        pf_model.reactions.get_by_id(rxn.id).upper_bound = 1000.

pf_model.reactions.get_by_id('PGK').lower_bound = -1000.
pf_model.reactions.get_by_id('PGK').upper_bound = 1000.

#logging.info('these mets arent used in more than one reaction and have no GPR:')
#for met in pf_model.metabolites:
#    if len(met.reactions) == 1:
#        for rxn in met.reactions:
#            logger.info(rxn.id)
#            if len(rxn.genes) < 1:
#                logging.info(met.id)

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
pf_model.metabolites.get_by_id('5mti_c').annotation['kegg.reaction'] = 'C19787'

# convert biomass to same aggregate format as for the rest of the models
pf_model.reactions.Protein.id = 'protein_bm'
pf_model.reactions.protein_bm.name = 'protein aggregate reaction for biomass'
pf_model.reactions.Lipid_prod.id = 'lipid_bm'
pf_model.reactions.lipid_bm.name   = 'lipid aggregate reaction for biomass'

# bm_rna = universal_model.metabolites.get_by_id('bm_rna_c').copy()
# bm_dna = universal_model.metabolites.get_by_id('bm_dna_c').copy()
# bm_lipid = universal_model.metabolites.get_by_id('bm_lipid_c').copy()
pf_model.metabolites.protein_c.id = 'bm_protein_c'
pf_model.metabolites.lipid_c.id = 'bm_lipid_c'

# add growth associated ATP demand
# atpm = universal_model.reactions.ATPM
# pf_model.add_reactions([atpm])
# pf_model.reactions.ATPM.lower_bound = 0.001
# pf_model.reactions.ATPM.upper_bound = 1000.

# some incorrectly formated exchange reactions
for rxn_id in ['EX_dag','EX_hb','EX_4ahmmp','EX_folate1','EX_inositol','EX_phosphatidyl1','EX_phosphatidyl2']:
    pf_model.remove_reactions([pf_model.reactions.get_by_id(rxn_id)])
pf_model.repair()

# some incorrectly formatteed exchange reactions
for rxn_id in [r.id for r in pf_model.reactions if r.id.endswith('_e_t')]:
    pf_model.reactions.get_by_id(rxn_id).id = rxn_id[:-4]+'t'

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
#logging.info(duplicates)

pf_model.repair()

# add full met info
met_counter = 0
for met in pf_model.metabolites:
    if met_counter % 10 == 0:
        logger.info(met_counter)
    met_counter = met_counter +1
    met_id = hf.met_ids_without_comp(pf_model,met.id)
    if met_id+'_c' in met_list: pf_model = hf.add_full_met_info(pf_model, met, met_id)
    else: pf_model = hf.add_partial_met_info(pf_model,met,met_id)
    if 'inchi' in met.annotation.keys():
        if met.annotation['inchi'][0] == 'nan': met.annotation['inchi'] = ['']

# add full reaction info
rxn_counter = 0
for rxn in pf_model.reactions:
    if rxn_counter % 10 == 0:
        logger.info(rxn_counter)
    rxn_counter = rxn_counter +1
    if rxn.id in rxn_list: pf_model = hf.add_full_rxn_info(pf_model, rxn, rxn.id)
    else: pf_model = hf.add_partial_rxn_info(pf_model, rxn, rxn.id)

# change id so that Memote doesn't think its the biomass reaction
pf_model.reactions.get_by_id('DM_biomass_c').name = 'unblock biomass'
pf_model.reactions.get_by_id('DM_biomass_c').id = 'DM_bm'

# fix some formatting issues that memote highlights
pf_model.reactions.ATPtm.annotation['kegg.reaction'] = 'R00124'
pf_model.metabolites.get_by_id('5mti_c').charge = int(pf_model.metabolites.get_by_id('5mti_c').charge)
pf_model.reactions.ATPtm.annotation['kegg.reaction'] = 'R00124' # incorrect as 'R00124#2'
pf_model.reactions.CHSTEROLt.annotation['rhea'] = ['39051', '39052', '39054', '39053'] # ['39051#1', '39052#1', '39054#1', '39053#1']
pf_model.reactions.OIVD1m.annotation['ec-code'] = ['1.2.1.25'] # ['1.2.1', '1.2.1.25']
pf_model.reactions.OIVD3m.annotation['ec-code'] = ['1.2.1.25'] #['1.2.1', '1.2.1.25']
pf_model.reactions.PPM.annotation['ec-code'] = ['5.4.2.7', '5.4.2.2'] # ['5.4.2', '5.4.2.7', '5.4.2.2']
pf_model.reactions.UDCPDPS.annotation['ec-code'] = ['2.5.1.31'] # ['2.5.1.M1', '2.5.1.31', '2.5.1']
pf_model.reactions.GLUTRS.annotation['ec-code'] = ['6.1.1.17', '6.1.1.24'] # ['6.1.1.17', '6.1.1', '6.1.1.24']
met_list_temp = ["adp_ap","adp_c","adp_m","cmp_ap","cmp_c","gdp_c","gdp_m","gdpfuc_c","gdpmann_c","malt_c","malt_e","uacgam_c","udp_c","udpg_c","gdp_ap"]
for met_id in met_list_temp:
    new_list_o_kegg_ids = list()
    if isinstance(pf_model.metabolites.get_by_id(met_id).annotation, dict):
        if 'kegg.compound' in pf_model.metabolites.get_by_id(met_id).annotation.keys():
            for option in pf_model.metabolites.get_by_id(met_id).annotation['kegg.compound']:
                if option.startswith('C'): # else starts with G -> glycan id
                    new_list_o_kegg_ids.append(option)
            pf_model.metabolites.get_by_id(met_id).annotation['kegg.compound'] = new_list_o_kegg_ids
        else:
            pf_model.metabolites.get_by_id(met_id).annotation['kegg.compound'] = []
pf_model.repair()

for rxn in pf_model.reactions:
    if 'bigg.reaction' in rxn.annotation.keys():
        temp_id = 'skip'
    elif 'bigg.reaction' in rxn.annotation.keys() and rxn.annotation['bigg.reaction'] != ['']:
        temp_id = 'skip'
    else:
        if rxn.id.endswith('ap'):
            temp_id = rxn.id[:-2] 
       	elif rxn.id.endswith('_ap'):
            temp_id = rxn.id[:-3] 
       	elif rxn.id.endswith('m'):
            temp_id = rxn.id[:-1]
       	elif rxn.id.endswith('_m'):
            temp_id = rxn.id[:-2]
        elif rxn.id.endswith('tmt'):
            temp_id = rxn.id[:-2]
        elif rxn.id.endswith('tap'):
            temp_id = rxn.id[:-2]
        elif rxn.id.endswith('t'):
            temp_id = rxn.id[:-1]
        else: temp_id = 'skip'

    if temp_id != 'skip':
        if temp_id in [r.id for r in universal_model.reactions]: 
            pf_model = hf.add_full_rxn_info(pf_model, rxn, temp_id) # in this case there is a reaction in a different compartmet in the universal
        elif temp_id+'tipp' in [r.id for r in universal_model.reactions]:
            pf_model = hf.add_full_rxn_info(pf_model, rxn, temp_id+'tipp') # this case is if there is an analogous periplasmic transport rxn in the universal
pf_model.repair()

# Rename these - keep reaction and associated info but have a new name
# the BIGG reaction with the same name is in a difference compartment or with slightly different cofactors
dict_for_new_comp = {'SBTR':'SBTRph',
'G3PDm':'G3PDm_pf',
'SELCYSS':'SELCYSS_pf',
'ACPS1':'ACPS1ap',
'DSBGGT':'DSBGGTc',
'DSBCGT':'DSBCGTc',
'DSBAO1':'DSBAO1e',
'CLPNH181pp':'CLPNH181',
'CLPNH180pp':'CLPNH180',
'CLPNH161pp':'CLPNH161',
'CLPNH160pp':'CLPNH160',
'CLPNH141pp':'CLPNH141',
'CLPNH140pp':'CLPNH140',
'CLPNH120pp':'CLPNH120',
'PLIPA2G181pp':'PLIPA2G181',
'PLIPA2G180pp':'PLIPA2G180',
'PLIPA2G161pp':'PLIPA2G161',
'PLIPA2G160pp':'PLIPA2G160',
'PLIPA2G141pp':"PLIPA2G141",
'PLIPA2G140pp':'PLIPA2G140',
'PLIPA2G120pp':'PLIPA2G120',
'PLIPA2E181pp':'PLIPA2E181',
'PLIPA2E141pp':'PLIPA2E141',
'PLIPA2E120pp':'PLIPA2E120',
'PLIPA2A181pp':'PLIPA2A181',
'PLIPA2A180pp':'PLIPA2A180',
'PLIPA2A161pp':'PLIPA2A161',
'PLIPA2A160pp':'PLIPA2A160',
'PLIPA2A141pp':'PLIPA2A141',
'PLIPA2A140pp':'PLIPA2A140',
'PLIPA2A120pp':'PLIPA2A120',
'LPLIPAL1G181pp':'LPLIPAL1G181',
'LPLIPAL1G180pp':'LPLIPAL1G180',
'LPLIPAL1G161pp':'LPLIPAL1G161',
'LPLIPAL1G160pp':'LPLIPAL1G160',
'LPLIPAL1G141pp':'LPLIPAL1G141',
'LPLIPAL1G140pp':'LPLIPAL1G140',
'LPLIPAL1G120pp':'LPLIPAL1G120',
'LPLIPAL1E141pp':'LPLIPAL1E141',
'LPLIPAL1E120pp':'LPLIPAL1E120',
'LPLIPAL1A181pp':'LPLIPAL1A181',
'LPLIPAL1A180pp':'LPLIPAL1A180',
'LPLIPAL1A161pp':'LPLIPAL1A161',
'LPLIPAL1A160pp':'LPLIPAL1A160',
'LPLIPAL1A141pp':'LPLIPAL1A141',
'LPLIPAL1A140pp':'LPLIPAL1A140',
'LPLIPAL1A120pp':'LPLIPAL1A120',
'PUNP8' : 'PUNP_pf',
'PI4P5K': 'PI4P5K_pf',
'GPIMTer_L': 'GPIMT',
'CBPSam':'CBPSac',
'DOLASNT':'DOLASNT_pf',
'PPPGOm':'PPPGOm_pf',
'LALDO2x':'LALDO2x_pf',
'GLUDy':'GLUDy_pf',
'PGK':'PGK_pf',
'DOLDPP':'DOLDPP_pf',
'RPE':'RPE_pf'}

for key, value in dict_for_new_comp.items():
    if key in [r.id for r in pf_model.reactions]:
        pf_model.reactions.get_by_id(key).id = value

for rxn in pf_model.reactions:
    if 'homo sapiens' in rxn.name:
        pf_model.reactions.get_by_id(rxn.id).name = rxn.name.replace("homo sapiens", "")
pf_model.repair()

# fix charges and formulas
pf_model = hf.fix_charge_or_formula(pf_model)

# fix food vacucole compartement
for met in pf_model.metabolites:
    if met.compartment == 'food vacuole':
        met.compartment = 'food_vacuole'

# TO DO: ask Memote to recognized EuPathDB gene IDs
for gene in pf_model.genes:
    gene.annotation['EuPathDB.genes'] = [gene.id]

# Add SBO terms
pf_model = hf.add_sbo_terms(pf_model)

pf_model.repair()

# fix a few SBO terms
for r_id in ["CHSTEROLt","CGMPt","SM_host","CAMPt"]:
    pf_model.reactions.get_by_id(r_id).annotation['sbo'] = 'SBO:0000185'
for r_id in ["perox_heme","HMGLB","CYOOm","CYOR_q8_m"]:
    pf_model.reactions.get_by_id(r_id).annotation['sbo'] = 'SBO:0000176'
pf_model.reactions.SK_fldox_ap.annotation['sbo'] = 'SBO:0000183'
pf_model.reactions.SK_trdrd_c.annotation['sbo'] = 'SBO:0000183'

# rename trna expression reactions as pseudo reactions and add expression to reaction name
for rxn in [r for r in pf_model.reactions if r.id.startswith('trna_')]:
    pf_model.reactions.get_by_id(rxn.id).name = 'trna expression ({})'.format(rxn.id)
    pf_model.reactions.get_by_id(rxn.id).annotation['sbo'] = 'SBO:0000183'
    pf_model.reactions.get_by_id(rxn.id).id = 'SK_{}'.format(rxn.reaction.split(' <=> ')[0])

for rxn in [r for r in pf_model.reactions if '_prod' in r.id]:
    pf_model.reactions.get_by_id(rxn.id).annotation['sbo'] = 'SBO:0000631' #pseudoreaction
    logger.info('added pseudoreaction SBO')
pf_model.repair()

# fix another issues identified with memote (inchi strings from bigg being listed as NaN)
for met in pf_model.metabolites:
    if 'inchi' in met.annotation.keys():
        if met.annotation['inchi'] == 'nan': del met.annotation['inchi']
    if 'reactome' in met.annotation.keys(): # no reactome ID in BiGG is compliant
        del met.annotation['reactome']
    pf_model.metabolites.get_by_id(met.id).annotation = met.annotation
for rxn in pf_model.reactions:
    if rxn.id in [r.id for r in universal_model.reactions] and 'bigg.reaction' not in rxn.annotation.keys():
        pf_model.reactions.get_by_id(rxn.id).annotation['bigg.reaction'] = rxn.id

pf_model.id = 'iPfal19_v1'
pf_model.name = 'iPfal19'
pf_model.compartments = {'cytosol': 'c', 'extracellular': 'e', 'apicoplast': 'ap', 'mitochondria': 'm', 'food_vacuole': 'fv'}
pf_model.notes = 'This model is the third iteration of the asexual blood-stage Plasmodium falciparum 3D7 genome-scale metabolic network reconstruction. The original reconstruction was generated using a custom pipeline by Plata et al (DOI: 10.1038/msb.2010.60) from P. falciparum Dd2 genome and curated to P. falciparum 3D7 and Dd2 function. Multiple rounds of curation were conducted (DOI: 10.1186/s12864-017-3905-1,10.1186/s12859-019-2756-y, and unpublished by Maureen Carey). Gene IDs can be mapped to sequences on https://plasmodb.org/ and reaction and metabolite nomenclature maps to data on http://bigg.ucsd.edu/.'
pf_model.repair()

os.chdir(model_path)
cobra.io.write_sbml_model(pf_model, "iPfal19_without_annotation.xml")

pf_model.annotation["taxonomy"] = "36329"
pf_model.annotation["genome"] = "https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/fasta/data/PlasmoDB-44_Pfalciparum3D7_Genome.fasta"
pf_model.annotation["DOI"] = "pending"
pf_model.annotation["authors"] = 'Maureen Carey, mac9jc@virginia.edu'
# cobrapy cannot currently handle the following:
#pf_model.annotation["authors"] =  [
#{"familyName": "Carey","givenName": "Maureen","organisation": "University of Virginia","email": "mac9jc@virginia.edu"},
#{"familyName": "Untaroiu","givenName": "Ana","organisation": "University of Virginia","email": "amu4pv@virginia.edu"}, 
#{"familyName": "Plata","givenName": "German","organisation": "Columbia University","email": "gap2118@columbia.edu"}]
pf_model.annotation["species"] = "Plasmodium falciparum"
pf_model.annotation["strain"] = "3D7"
pf_model.annotation["tissue"] = "parasite in the asexual blood-stage"
pf_model.annotation["terms_of_distribution"] = "CC-BY" 
pf_model.annotation["created"] = day
pf_model.annotation["sbo"] = "SBO:0000624"
pf_model.annotation["curation"] = ['DOI: 10.1038/msb.2010.60', 'DOI: 10.1186/s12864-017-3905-1', 'DOI: 10.1186/s12859-019-2756-y', 'unpublished by Maureen Carey']
pf_model.annotation["genedb"] = "Pfalciparum"

logger.info({'can the model grow?':pf_model.slim_optimize()})

cobra.io.save_json_model(pf_model, "iPfal19.json")
cobra.io.write_sbml_model(pf_model, "iPfal19.xml")
