#!/usr/bin/env python3
import os
import sys
import glob
from openbabel import openbabel
from rdkit import RDLogger

def convert_out_to_pdb(input_folder, output_folder, nonmatches_file="incorrect.txt"):
    os.makedirs(output_folder, exist_ok=True)
    RDLogger.DisableLog('rdApp.*')  # suppress RDKit warnings

    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats("gamout", "pdb")

    nonmatching_files = []

    for file_path in glob.glob(os.path.join(input_folder, "*.out")):
        mol_ob = openbabel.OBMol()
        if not conv.ReadFile(mol_ob, file_path):
            print(f"Failed to read {file_path}")
            nonmatching_files.append(file_path)
            continue

        pdb_filename = os.path.splitext(os.path.basename(file_path))[0] + ".pdb"
        pdb_path = os.path.join(output_folder, pdb_filename)

        if not conv.WriteFile(mol_ob, pdb_path):
            print(f"Failed to write {pdb_filename}")
            nonmatching_files.append(file_path)
            continue

        #print(f"Converted {file_path} → {pdb_path}")

    if nonmatching_files:
        with open(nonmatches_file, "w") as f:
            for item in nonmatching_files:
                f.write(f"{item}\n")
        print(f"\nSome files could not be converted. See {nonmatches_file}")
    else:
        print("\n All files converted successfully!")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: convert_out_to_pdb.py <input_folder> <output_folder>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    convert_out_to_pdb(input_folder, output_folder)

