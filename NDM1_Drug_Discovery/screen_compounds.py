import os
import pandas as pd
import json
from rdkit import Chem
from rdkit.Chem import Descriptors

def screen_compounds(csv_path, output_json_path):
    print(f"Loading compounds from {csv_path}...")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Input file {csv_path} does not exist.")
        
    df = pd.read_csv(csv_path)
    screened_leads = []
    rejected_compounds = []
    
    for idx, row in df.iterrows():
        comp_id = row['Compound_ID']
        smiles = row['SMILES']
        description = row['Description']
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  {comp_id}: Invalid SMILES structure. Skipping.")
            continue
            
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Descriptors.MolLogP(mol), 2)
        
        # Lipinski criteria
        mw_pass = mw <= 500
        logp_pass = logp <= 5.0
        
        comp_data = {
            "Compound_ID": comp_id,
            "SMILES": smiles,
            "Description": description,
            "MW": mw,
            "LogP": logp,
            "MW_Pass": mw_pass,
            "LogP_Pass": logp_pass,
            "Passed_Screening": mw_pass and logp_pass
        }
        
        if mw_pass and logp_pass:
            screened_leads.append(comp_data)
            print(f"  {comp_id}: PASSED (MW={mw}, LogP={logp})")
        else:
            rejected_compounds.append(comp_data)
            print(f"  {comp_id}: FAILED (MW={mw}, LogP={logp})")
            
    # Save the screened leads to JSON
    os.makedirs(os.path.dirname(output_json_path), exist_ok=True)
    with open(output_json_path, 'w') as jf:
        json.dump(screened_leads, jf, indent=4)
        
    print(f"Saved {len(screened_leads)} screened leads to {output_json_path}")
    
    # Save a detailed log of all compounds for report generation
    with open("data/screening_log.json", 'w') as jf:
        json.dump({
            "screened_leads": screened_leads,
            "rejected_compounds": rejected_compounds
        }, jf, indent=4)
        
    return screened_leads, rejected_compounds

if __name__ == "__main__":
    csv_path = "data/input_compounds.csv"
    output_json_path = "data/screened_leads.json"
    screen_compounds(csv_path, output_json_path)
