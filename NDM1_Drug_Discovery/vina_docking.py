import os
import json
import subprocess

def run_command(cmd):
    print(f"Running command: {cmd}")
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Error executing command: {cmd}")
        print(f"Stdout:\n{result.stdout}")
        print(f"Stderr:\n{result.stderr}")
        return False, result.stdout, result.stderr
    return True, result.stdout, result.stderr

def parse_vina_score(pdbqt_path):
    if not os.path.exists(pdbqt_path):
        return None
    with open(pdbqt_path, 'r') as f:
        for line in f:
            if "REMARK VINA RESULT:" in line:
                # Example: REMARK VINA RESULT:      -5.8      0.000      0.000
                parts = line.split()
                if len(parts) >= 4:
                    return float(parts[3])
    return None

def main():
    screened_json = "data/screened_leads.json"
    receptor_pdb = "data/3RGG.pdb"
    receptor_pdbqt = "data/3RGG.pdbqt"
    final_output = "data/final_docking_results.json"
    
    print("Initializing Layer 2 Molecular Docking using AutoDock Vina...")
    
    # 1. Read screened leads
    if not os.path.exists(screened_json):
        raise FileNotFoundError(f"{screened_json} not found. Please run Layer 1 screening first.")
    with open(screened_json, 'r') as f:
        leads = json.load(f)
        
    # 2. Check and prepare receptor
    if not os.path.exists(receptor_pdbqt):
        print(f"Receptor PDBQT file {receptor_pdbqt} not found. Preparing receptor using Open Babel...")
        success, stdout, stderr = run_command(f"obabel {receptor_pdb} -O {receptor_pdbqt} -xr")
        if not success:
            raise RuntimeError("Failed to prepare receptor PDBQT.")
    else:
        print(f"Using existing receptor PDBQT: {receptor_pdbqt}")
        
    # 3. Active site coordinates (centroid of AIR ligand in 3RGG.pdb)
    # Centroid: (7.314, -38.596, 14.329)
    center_x = 7.314
    center_y = -38.596
    center_z = 14.329
    size_x = 20.0
    size_y = 20.0
    size_z = 20.0
    
    docking_results = []
    
    # 4. Dock each ligand in sequence
    for idx, lead in enumerate(leads):
        comp_id = lead["Compound_ID"]
        smiles = lead["SMILES"]
        desc = lead["Description"]
        
        print(f"\nProcessing {comp_id} ({desc})...")
        
        # Define file paths
        ligand_pdbqt = f"data/ligand_{comp_id}.pdbqt"
        docked_pdbqt = f"data/docked_{comp_id}.pdbqt"
        log_file = f"data/docking_{comp_id}.log"
        
        # 4a. Convert SMILES to 3D PDBQT using Open Babel
        print(f"  Generating 3D conformation and PDBQT for {comp_id}...")
        ob_cmd = f"obabel -:\"{smiles}\" --gen3d -O {ligand_pdbqt}"
        success, stdout, stderr = run_command(ob_cmd)
        if not success or not os.path.exists(ligand_pdbqt):
            print(f"  [ERROR] Open Babel conversion failed for {comp_id}. Skipping.")
            continue
            
        # 4b. Execute AutoDock Vina with CPU limit = 2
        print(f"  Running AutoDock Vina docking for {comp_id}...")
        vina_cmd = (
            f"vina --receptor {receptor_pdbqt} --ligand {ligand_pdbqt} "
            f"--center_x {center_x} --center_y {center_y} --center_z {center_z} "
            f"--size_x {size_x} --size_y {size_y} --size_z {size_z} "
            f"--cpu 2 --out {docked_pdbqt} --exhaustiveness 8 > {log_file} 2>&1"
        )
        success, stdout, stderr = run_command(vina_cmd)
        
        # 4c. Parse scores
        score = parse_vina_score(docked_pdbqt)
        if score is not None:
            print(f"  [SUCCESS] {comp_id} docked successfully. Affinity Score: {score} kcal/mol")
        else:
            print(f"  [ERROR] Failed to extract docking score for {comp_id}.")
            score = 999.0  # Placeholder for failed docking
            
        docking_results.append({
            "Compound_ID": comp_id,
            "SMILES": smiles,
            "Description": desc,
            "MW": lead["MW"],
            "LogP": lead["LogP"],
            "Binding_Affinity_kcal_mol": score,
            "Docked_Pose_Path": docked_pdbqt,
            "Vina_Log_Path": log_file
        })
        
        # Clean up temporary ligand preparation file
        if os.path.exists(ligand_pdbqt):
            os.remove(ligand_pdbqt)
            
    # Sort results by binding affinity (lowest energy is best)
    docking_results.sort(key=lambda x: x["Binding_Affinity_kcal_mol"])
    
    # 5. Save master results file
    with open(final_output, 'w') as f:
        json.dump(docking_results, f, indent=4)
        
    print(f"\nCompleted docking for all leads. Master results saved to {final_output}")

if __name__ == "__main__":
    main()
