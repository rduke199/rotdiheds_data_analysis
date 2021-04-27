from pymongo import MongoClient
import pandas as pd

mc = MongoClient("mongodb://wf:wf200@128.163.204.13:27017/ocelot_wf")
db = mc["ocelot_wf"]
# Get smiles data from database (run only once...then use generated csv file)
molsconf_collection  = db["DSmolconf"]
print("Starting to collect results")
molsconf_results = molsconf_collection.find({"geometries.chrombackbone.smiles": {"$exists":"true"}})
print("Done getting results")

molsconf_data = pd.DataFrame.from_records(molsconf_results)
print("Done making data frame")
geometry_data = pd.json_normalize(molsconf_data.geometries, max_level=0)
chrombackbone_data = pd.json_normalize(geometry_data.chrombackbone, max_level=0)
chrombackbone_data.set_index(molsconf_data._id, inplace=True)
smiles_data = chrombackbone_data.smiles
smiles_data.to_csv('smiles_data.csv')