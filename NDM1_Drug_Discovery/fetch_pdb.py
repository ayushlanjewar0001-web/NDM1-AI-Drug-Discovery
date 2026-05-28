import os
import urllib.request

def download_pdb(pdb_id, output_path):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"Downloading {pdb_id} from {url}...")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    urllib.request.urlretrieve(url, output_path)
    print(f"Saved to {output_path}")

def parse_active_site(pdb_path):
    print(f"Parsing {pdb_path}...")
    zinc_coords = []
    coordinating_residues = []
    resolution = "Unknown"
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("REMARK   2 RESOLUTION."):
                parts = line.split()
                if len(parts) >= 4:
                    resolution = parts[3]
            elif line.startswith("HETATM") and "ZN" in line:
                # Extract ZN coordinate details
                # PDB HETATM format:
                # 31-38: X coordinate
                # 39-46: Y coordinate
                # 47-54: Z coordinate
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                res_seq = int(line[22:26].strip())
                zinc_coords.append({
                    "atom": atom_name,
                    "res_name": res_name,
                    "chain": chain_id,
                    "res_seq": res_seq,
                    "coords": (x, y, z)
                })
            # Let's also look for coordinating residues if we can find them
            # or print out some nearby residues.
            
    # Compute centroid of Zn atoms
    if zinc_coords:
        xs = [z["coords"][0] for z in zinc_coords]
        ys = [z["coords"][1] for z in zinc_coords]
        zs = [z["coords"][2] for z in zinc_coords]
        centroid = (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))
    else:
        centroid = None

    print(f"Resolution: {resolution} A")
    print(f"Found {len(zinc_coords)} Zinc atoms:")
    for z in zinc_coords:
        print(f"  {z['atom']} ({z['res_name']} {z['res_seq']} Chain {z['chain']}): {z['coords']}")
    if centroid:
        print(f"Centroid of Zinc atoms: {centroid}")
        
    return resolution, zinc_coords, centroid

if __name__ == "__main__":
    pdb_id = "3RGG"
    output_path = "data/3RGG.pdb"
    if not os.path.exists(output_path):
        download_pdb(pdb_id, output_path)
    resolution, zincs, centroid = parse_active_site(output_path)
    
    # Save the coordinates details to a file for subsequent steps
    import json
    with open("data/active_site_coords.json", "w") as jf:
        json.dump({
            "pdb_id": pdb_id,
            "resolution": resolution,
            "zinc_atoms": zincs,
            "active_site_centroid": centroid
        }, jf, indent=4)
    print("Saved active site coordinate data to data/active_site_coords.json")
