import cobra
from cobra import Model, Reaction, Metabolite
import os

model_path = "/home/mac9jc/paradigm/models"
os.chdir(model_path)
model = cobra.io.load_json_model("iPfal19_updated.json")
print(model.slim_optimize())

print('made up met')
heme2_degraded_fv = Metabolite('heme2_degraded_fv',formula='',
    name='degraded heme',compartment='fv')
model.add_metabolites([heme2_degraded_fv])
model.add_boundary(model.metabolites.get_by_id('heme2_degraded_fv'), type = "demand")
model.objective = 'DM_heme2_degraded_fv'
model.reactions.get_by_id('DM_heme2_degraded_fv').lower_bound = 0.001
sol = model.slim_optimize()
print('sol < 0.01')
print(sol < 0.01)
print(sol)

print('3ocpalm9eACP_ap')
model.add_boundary(model.metabolites.get_by_id('3ocpalm9eACP_ap'), type = "demand")
model.objective = 'DM_3ocpalm9eACP_ap'
sol = model.slim_optimize()
print('sol < 0.01')
print(sol < 0.01)
print(sol)
