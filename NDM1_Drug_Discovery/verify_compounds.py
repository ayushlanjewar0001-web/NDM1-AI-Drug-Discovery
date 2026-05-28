import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

def verify_compounds(csv_path):
    print(f"Reading compounds from {csv_path}...")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"{csv_path} not found.")
        
    df = pd.read_csv(csv_path)
    results = []
    
    for idx, row in df.iterrows():
        comp_id = row['Compound_ID']
        smiles = row['SMILES']
        description = row['Description']
        
        mol = Chem.MolFromSmiles(smiles)
        is_valid = mol is not None
        
        comp_data = {
            "Compound_ID": comp_id,
            "SMILES": smiles,
            "Description": description,
            "RDKit_Valid": is_valid,
            "MW": None,
            "LogP": None,
            "HBD": None,
            "HBA": None,
            "TPSA": None,
            "Formula": None
        }
        
        if is_valid:
            comp_data["MW"] = round(Descriptors.MolWt(mol), 2)
            comp_data["LogP"] = round(Descriptors.MolLogP(mol), 2)
            comp_data["HBD"] = Descriptors.NumHDonors(mol)
            comp_data["HBA"] = Descriptors.NumHAcceptors(mol)
            comp_data["TPSA"] = round(Descriptors.TPSA(mol), 2)
            comp_data["Formula"] = Chem.rdMolDescriptors.CalcMolFormula(mol)
            print(f"  {comp_id}: Valid. Formula={comp_data['Formula']}, MW={comp_data['MW']}")
        else:
            print(f"  {comp_id}: INVALID SMILES!")
            
        results.append(comp_data)
        
    res_df = pd.DataFrame(results)
    res_csv = "data/verified_compounds.csv"
    res_df.to_csv(res_csv, index=False)
    print(f"Saved verified compounds to {res_csv}")
    
    # Save a JSON copy for ease of parsing by orchestrator/subagents
    res_json = "data/verified_compounds.json"
    res_df.to_json(res_json, orient="records", indent=4)
    print(f"Saved verified compounds JSON to {res_json}")
    
    return results

if __name__ == "__main__":
    csv_path = "data/input_compounds.csv"
    verify_compounds(csv_path)
