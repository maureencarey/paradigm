import pandas as pd
from cobra.flux_analysis.parsimonious import add_pfba

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

# pFBA based gapfilling, implementation from Greg Medlock
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
    get_fluxes = set([r.id for r in universal_model_pfba.reactions]) - set([rxn.id for rxn in input_model.reactions]) ### ERROR IN ORIGINAL CODE HAD MODEL instead of INPUT_MODEL
    print(solution)
    print(dir(solution))
    add_reactions_to_model = [rxn for rxn in get_fluxes if abs(solution.fluxes[rxn]) > 1E-8]

    # double check
    logging.info('double checking pFBA solution')
        test_model = input_model.copy()
        add_reactions_list = [universal_model_pfba.reactions.get_by_id(r).copy() for r in add_reactions_to_model]
        test_model.add_reactions(add_reactions_list)
        sol = test_model.optimize()
        if sol.status == 'infeasible':
            logging.info('pFBA solution is infeasible!')
    return(add_reactions_to_model)

def intersection(rxn_list_compartments, acceptable_compartments):
    temp = set(acceptable_compartments)
    unacceptable_comp = [value for value in rxn_list_compartments if value not in temp]
    return(unacceptable_comp)
