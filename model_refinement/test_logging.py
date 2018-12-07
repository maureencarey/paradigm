import cobra
import os
from cobra import Model, Reaction, Metabolite
import pandas as pd
import requests
import logging
import argparse
from datetime import datetime


parser = argparse.ArgumentParser(description='Read in the species model')
parser.add_argument('model_file')
args = parser.parse_args()
model_fname = vars(args)['model_file']

# parse arguments for global variables
SPECIES_ID = model_fname

day = datetime.now().strftime('%d_%m_%Y')

logging.basicConfig(filename='test_logging_{}_{}.log'.format(SPECIES_ID,day), level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)

og_path = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm"
data_path = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data"
model_path = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/models"

os.chdir(model_path)
pf_model = cobra.io.read_sbml_model("iPfal17.xml")
logging.info('finished loading model')

model1 = pf_model.copy()
model1.remove_reactions([model1.reactions.get_by_id(r.id) for r in pf_model.reactions[1:99]])
pf_model.remove_reactions(pf_model.reactions[1:100])

from cobra.flux_analysis.parsimonious import add_pfba
def pfba_gapfill_implementation(input_model, universal_model_ex, objective_reaction_id):
    # objective_reaction is a reaction id

    universal_model_pfba = universal_model_ex.copy()
    add_pfba(universal_model_pfba)

    # penalize adding demand reactions
    coef = universal_model_pfba.objective.get_linear_coefficients(universal_model_pfba.variables)
    for key,value in coef.items():
        if key.name.startswith('DM_') or key.name.startswith('SK_'):
            coef[key] = 1000.
        elif key.name.startswith('EX_'):
            coef[key] = 50.
    for x in universal_model_pfba.variables:
        if x.name.startswith('DM_') or x.name.startswith('SK_') or x.name.startswith('EX_'):
            universal_model_pfba.objective.set_linear_coefficients(coef)

    rxns_to_remove = [rxn for rxn in universal_model_pfba.reactions if rxn.id \
                      in [rxn.id for rxn in input_model.reactions]]
    universal_model_pfba.remove_reactions(rxns_to_remove)
    universal_model_pfba.add_reactions([rxn for rxn in input_model.reactions])
    universal_model_pfba.reactions.get_by_id(objective_reaction_id).lower_bound = 0.05

    solution = universal_model_pfba.optimize()
    if solution.status == 'infeasible':
        logging.info('pFBA gapfilling for {} is infeasible!'.format(objective_reaction_id))
    else:
        logging.info('pFBA gapfilling - feasible')
    get_fluxes = set([r.id for r in universal_model_pfba.reactions]) - set([rxn.id for rxn in input_model.reactions]) ### ERROR IN ORIGINAL CODE
    add_reactions_to_model = [rxn for rxn in get_fluxes if abs(solution.x_dict[rxn]) > 1E-8]

    # double check
    logging.info('double checking pFBA solution')
    test_model = input_model.copy()
    add_reactions_list = [universal_model_pfba.reactions.get_by_id(r).copy() for r in add_reactions_to_model]
    test_model.add_reactions(add_reactions_list)
    sol = test_model.optimize()
    if sol.status == 'infeasible':
        logging.info('pFBA solution is infeasible!')
    return(add_reactions_to_model)

pfba_gapfill_implementation(pf_model, model1, 'biomass')


