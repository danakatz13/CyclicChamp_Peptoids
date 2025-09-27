import os
import glob
from openbabel import openbabel
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit import RDLogger

input_folder = "./incomplete"
output_folder = "."
nonmatches_file = "incorrect.txt"
smarts_pattern = Chem.MolFromSmarts(
    '[C;H3;X4][C](=[O])[N:1]([C])[C:2][C:3](=[O])[N:9]([C:10])[C:11]')

os.makedirs(output_folder, exist_ok=True)
RDLogger.DisableLog('rdApp.*')  # suppress RDKit warnings

conv = openbabel.OBConversion()
conv.SetInAndOutFormats("gamout", "pdb")
nonmatching_files = []

for file_path in glob.glob(os.path.join(input_folder, "*.out")):
    mol_ob = openbabel.OBMol()
    if not conv.ReadFile(mol_ob, file_path):
        print(f"Failed to read {file_path}")
        continue

    pdb_filename = os.path.splitext(os.path.basename(file_path))[0] + ".pdb"
    pdb_path = os.path.join(output_folder, pdb_filename)
    conv.WriteFile(mol_ob, pdb_path)

    mol = Chem.MolFromPDBFile(pdb_path, sanitize=False, removeHs=False)
    if mol is None:
        print(f"Failed to load {pdb_filename} in RDKit")
        continue

    matches = mol.GetSubstructMatches(smarts_pattern)
    if matches:
        continue
    else:
        print(f"{pdb_filename} does NOT match the SMARTS pattern")
        nonmatching_files.append(pdb_filename)
        '''
        img = Draw.MolToImage(mol)
        img.show()
        '''
if nonmatching_files:
    with open(nonmatches_file, "w") as f:
        for filename in nonmatching_files:
            f.write(filename + "\n")

print(f"\n {len(nonmatching_files)} files did NOT match. Saved to {nonmatches_file}.")

