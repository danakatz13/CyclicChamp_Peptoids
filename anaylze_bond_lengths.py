import sys
import os
import csv
from io import BytesIO
from openbabel import openbabel
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem.rdMolTransforms as MolTransforms
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image


already_shown = False

def draw_and_show_molecule(mol, filename):
    global already_shown
    if already_shown:
        return  # skip for all other PDBs
    already_shown = True

    print(f"  -> Showing molecule for {filename}...")
    pil_image = None
    try:
        drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
        drawer.drawOptions().addAtomIndices = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        img_bytes = drawer.GetDrawingText()
        pil_image = Image.open(BytesIO(img_bytes))
    except Exception as e:
        print(f"     Failed to draw molecule: {e}")
    
    if pil_image:
        pil_image.show()

def get_bond_lengths_from_pdb(mol, conf):
    bond_lengths = {}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        length = MolTransforms.GetBondLength(conf, a1, a2)
        bond_lengths[f"{a1}-{a2}"] = length
    return bond_lengths


def get_bond_angles_from_pdb(mol, conf):
    angles = {}
    for atom in mol.GetAtoms():
        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                a1, a2, a3 = neighbors[i], atom.GetIdx(), neighbors[j]
                angle = MolTransforms.GetAngleDeg(conf, a1, a2, a3)
                angles[f"{a1}-{a2}-{a3}"] = angle
    return angles


def get_dihedrals_from_pdb(mol, conf):
    torsions = {}
    for bond in mol.GetBonds():
        a2 = bond.GetBeginAtom()
        a3 = bond.GetEndAtom()
        for nbr1 in a2.GetNeighbors():
            if nbr1.GetIdx() == a3.GetIdx():
                continue
            for nbr2 in a3.GetNeighbors():
                if nbr2.GetIdx() == a2.GetIdx():
                    continue
                a1, a4 = nbr1.GetIdx(), nbr2.GetIdx()
                dihedral = MolTransforms.GetDihedralDeg(conf, a1, a2.GetIdx(), a3.GetIdx(), a4)
                torsions[f"{a1}-{a2.GetIdx()}-{a3.GetIdx()}-{a4}"] = dihedral
    return torsions


def write_bond_and_angle_table(output_csv, pdb_file_list, all_bond_data, all_angle_data):
    """
    Writes both bond lengths and bond angles into one pivoted table.
    Rows = bond or angle identifiers
    Columns = pdb structures
    """
    # Collect all keys
    all_bonds = sorted({k for d in all_bond_data.values() for k in d})
    all_angles = sorted({k for d in all_angle_data.values() for k in d})

    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)

        # Header row: Type + ID + PDB filenames
        header = ["Type", "ID"] + [os.path.basename(p) for p in pdb_file_list]
        writer.writerow(header)

        # Write bond lengths
        for bond in all_bonds:
            row = ["Bond", bond]
            for pdb_path in pdb_file_list:
                value = all_bond_data.get(pdb_path, {}).get(bond, "")
                if value != "":
                    value = f"{value:.3f}"
                row.append(value)
            writer.writerow(row)

        # Write bond angles
        for angle in all_angles:
            row = ["Angle", angle]
            for pdb_path in pdb_file_list:
                value = all_angle_data.get(pdb_path, {}).get(angle, "")
                if value != "":
                    value = f"{value:.3f}"
                row.append(value)
            writer.writerow(row)

    print(f"\n Bond-lengths and angles table written to {output_csv}")


import os
import sys
import glob
import csv
from rdkit import Chem, RDLogger


# ensure output folder exists
output_folder = "./output"
os.makedirs(output_folder, exist_ok=True)


# suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')


# set up OpenBabel conversion: GAMESS output (.out) → PDB
conv = openbabel.OBConversion()
conv.SetInAndOutFormats("gamout", "pdb")


nonmatching_files = []


def convert_out_to_pdb(input_folder, output_folder):
    pdb_files = []
    for file_path in glob.glob(os.path.join(input_folder, "*.out")):
        mol_ob = obabel.OBMol()
        if not conv.ReadFile(mol_ob, file_path):
            print(f"Failed to read {file_path}")
            continue


        pdb_filename = os.path.splitext(os.path.basename(file_path))[0] + ".pdb"
        pdb_path = os.path.join(output_folder, pdb_filename)
        if not conv.WriteFile(mol_ob, pdb_path):
            print(f"Failed to write {pdb_filename}")
            continue


        pdb_files.append(pdb_path)
        return pdb_files
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
            "Usage: python analyze_bond_lengths.py <output_csv_file> <pdb_file_1> <pdb_file_2> ...",
            file=sys.stderr,
        )
        sys.exit(1)

    output_csv = sys.argv[1]
    pdb_file_list = sys.argv[2:]

    all_bond_data = {}
    all_angle_data = {}
    all_dihedral_data = {}

    for pdb_path in pdb_file_list:
        print(f"\n--- Processing: {os.path.basename(pdb_path)} ---")
        if not os.path.exists(pdb_path):
            print(f"Error: PDB file '{pdb_path}' not found.", file=sys.stderr)
            continue
        mol = Chem.MolFromPDBFile(pdb_path, sanitize=True, removeHs=False)
        if mol is None:
            print(f"Error: Could not read molecule from {pdb_path}.", file=sys.stderr)
            continue

        mol = Chem.AddHs(mol)

        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, randomSeed=0xf00d)

        draw_and_show_molecule(mol, os.path.basename(pdb_path))

        if mol.GetNumConformers() == 0:
            print(
                f"Warning: No 3D conformation in {pdb_path}. Skipping measurements.",
                file=sys.stderr,
            )
            all_bond_data[pdb_path] = {}
            all_angle_data[pdb_path] = {}
            all_dihedral_data[pdb_path] = {}
            continue

        conf = mol.GetConformer()
        all_bond_data[pdb_path] = get_bond_lengths_from_pdb(mol, conf)
        all_angle_data[pdb_path] = get_bond_angles_from_pdb(mol, conf)
        all_dihedral_data[pdb_path] = get_dihedrals_from_pdb(mol, conf)

    # ✅ This check happens AFTER the loop
    if not any([all_bond_data, all_angle_data]):
        print("\nNo data was extracted. CSV file will not be created.", file=sys.stderr)
        sys.exit(1)

    # Write combined output
    write_bond_and_angle_table(output_csv, pdb_file_list, all_bond_data, all_angle_data)

