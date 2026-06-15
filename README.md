# NDM-1 Superbug AI Drug Discovery Pipeline

An industrial-grade computational drug discovery pipeline built to screen inhibitors against **New Delhi metallo-beta-lactamase 1 (NDM-1)** using automated multi-agent architecture and python runtimes.

# 🛠️ Architecture & Tech Stack
- **Framework:** Multi-Agent Automation Framework (Python Sandbox)
- **Cheminformatics:** RDKit (Structure validation & Lipinski filtering)
- **Physics Engine:** AutoDock Vina (High-throughput molecular docking)
- **Target Protein:** NDM-1 enzyme (PDB ID: 3RGG)

# 🔄 Project Workflow
1. **Target Identification:** Automated extraction of binding pocket coordinates for 3RGG via RCSB PDB API.
2. **Layer 1 Filter:** Data validation of input SMILES strings and structural filtration violating Lipinski's Rule of 5.
3. **Layer 2 Docking:** Execution of local high-throughput docking with resource throttling (`--cpu 2`) to optimize host machine compute.
4. **Optimization:** Generation of structured Lead Optimization logs tracking precise thermodynamic binding affinity scores ($kcal/mol$).
