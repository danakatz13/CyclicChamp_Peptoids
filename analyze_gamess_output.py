import os
import sys
import re
import csv
import math

def normalize_angle(angle):
    """
    Normalizes an angle to be within the range (-180, 180].
    """
    normalized = math.fmod(angle, 360.0)
    if normalized > 180.0:
        normalized -= 360.0
    elif normalized <= -180.0:
        normalized += 360.0
    return normalized

def calculate_periodic_deviation(angle1, angle2):
    """
    Calculates the minimum deviation between two angles, considering periodic boundaries.
    Angles are expected to be in degrees.
    """
    norm_angle1 = normalize_angle(angle1)
    norm_angle2 = normalize_angle(angle2)

    diff = abs(norm_angle1 - norm_angle2)
    return min(diff, 360.0 - diff)


def parse_initial_dihedrals_from_filename(filename):
    """
    Parses initial PHI, PSI, OMEGA values from a filename like
    'phip040_psim160_omegap170.out'.
    Assumes 'p' for positive, 'm' for negative.
    Returns a dictionary with 'phi', 'psi', 'omega' as floats, or None if not found.
    """
    initial_dihedrals = {}

    phi_match = re.search(r'phi([pm])(\d+)', filename, re.IGNORECASE)
    psi_match = re.search(r'psi([pm])(\d+)', filename, re.IGNORECASE)
    omega_match = re.search(r'omega([pm])(\d+)', filename, re.IGNORECASE)

    if phi_match:
        sign = 1 if phi_match.group(1).lower() == 'p' else -1
        initial_dihedrals['phi'] = float(phi_match.group(2)) * sign
    else:
        initial_dihedrals['phi'] = None

    if psi_match:
        sign = 1 if psi_match.group(1).lower() == 'p' else -1
        initial_dihedrals['psi'] = float(psi_match.group(2)) * sign
    else:
        initial_dihedrals['psi'] = None

    if omega_match:
        sign = 1 if omega_match.group(1).lower() == 'p' else -1
        initial_dihedrals['omega'] = float(omega_match.group(2)) * sign
    else:
        initial_dihedrals['omega'] = None
    
    return initial_dihedrals


def extract_and_evaluate_file(filepath, deviation_threshold):
    """
    Processes a single .out file, extracts final dihedrals and energy,
    parses initial dihedrals from filename, compares them, and determines if flagged.
    
    Returns a tuple: (extracted_data_dict, flagged_data_dict_or_None)
    - extracted_data_dict: Contains filename, final PHI, PSI, OMEGA, ENERGY.
    - flagged_data_dict_or_None: Contains detailed deviation info if flagged, else None.
    """
    # Initialize all values to None
    phi_val = None
    psi_val = None
    omega_val = None
    final_energy = None

    file_basename = os.path.basename(filepath) 

    try:
        with open(filepath, 'r') as f:
            in_restraint_block = False
            for line in f:
                stripped_line = line.strip()

                # Search for the TOTAL ENERGY IN SOLVENT line anywhere in the file.
                if "TOTAL ENERGY IN SOLVENT" in line:
                    energy_match = re.search(r"=\s*(-?\d+\.\d+)", line)
                    if energy_match:
                        # --- MODIFICATION: Convert Hartrees to kcal/mol ---
                        energy_in_hartrees = float(energy_match.group(1))
                        # Conversion factor: 1 Hartree = 627.5095 kcal/mol
                        final_energy = energy_in_hartrees * 627.5095

                # Original logic to find dihedral values
                if "RESTRAINTS STATUS" in line:
                    in_restraint_block = True
                    continue 
                elif in_restraint_block:
                    if stripped_line.startswith("DIH"):
                        parts = stripped_line.split()
                        if len(parts) >= 5:
                            idx1, idx2, idx3, idx4 = parts[1:5]
                            try:
                                value = float(parts[5]) 
                            except ValueError:
                                continue 

                            if idx1 == '2' and idx2 == '4' and idx3 == '6' and idx4 == '7':
                                phi_val = value
                            elif idx1 == '4' and idx2 == '6' and idx3 == '7' and idx4 == '9':
                                psi_val = value
                            elif idx1 == '1' and idx2 == '2' and idx3 == '4' and idx4 == '6':
                                omega_val = value
                    # The block to stop reading restraints is still useful
                    elif "MAXIMUM GRADIENT" in line or "RMS GRADIENT" in line:
                        in_restraint_block = False
                        
    except Exception as e:
        print(f"Error parsing file {filepath}: {e}", file=sys.stderr)
        return None, None # Return None for both if parsing fails

    # Prepare extracted data for the first output file
    extracted_data = {
        'filename': file_basename,
        'phi': phi_val,
        'psi': psi_val,
        'omega': omega_val,
        'energy': final_energy
    }

    # If GAMESS-calculated values for dihedrals are not all found, we cannot perform comparison
    # but we still return the extracted_data (which will have None for missing values)
    if phi_val is None or psi_val is None or omega_val is None:
        print(f"Warning: Missing PHI/PSI/OMEGA values in {filepath}. Skipping deviation check.", file=sys.stderr)
        return extracted_data, None

    # Now, perform the comparison for flagging
    initial_dihedrals = parse_initial_dihedrals_from_filename(file_basename)
    initial_phi = initial_dihedrals.get('phi')
    initial_psi = initial_dihedrals.get('psi')
    initial_omega = initial_dihedrals.get('omega')

    deviations = {}
    is_flagged = False

    # Compare PHI with periodic boundary conditions
    if initial_phi is not None:
        dev_phi = calculate_periodic_deviation(phi_val, initial_phi)
        deviations['phi'] = dev_phi
        if dev_phi > deviation_threshold:
            is_flagged = True
    else:
        deviations['phi'] = "N/A (initial PHI not found)"

    # Compare PSI with periodic boundary conditions
    if initial_psi is not None:
        dev_psi = calculate_periodic_deviation(psi_val, initial_psi)
        deviations['psi'] = dev_psi
        if dev_psi > deviation_threshold:
            is_flagged = True
    else:
        deviations['psi'] = "N/A (initial PSI not found)"

    # Compare OMEGA with periodic boundary conditions
    if initial_omega is not None:
        dev_omega = calculate_periodic_deviation(omega_val, initial_omega)
        deviations['omega'] = dev_omega
        if dev_omega > deviation_threshold:
            is_flagged = True
    else:
        deviations['omega'] = "N/A (initial OMEGA not found)"
    
    flagged_data = None
    if is_flagged:
        flagged_data = {
            'filename': file_basename,
            'initial_phi': initial_phi,
            'initial_psi': initial_psi,
            'initial_omega': initial_omega,
            'final_phi': phi_val,
            'final_psi': psi_val,
            'final_omega': omega_val,
            'final_energy': final_energy, # This will be the converted energy
            'deviation_phi': f"{deviations['phi']:.2f}" if isinstance(deviations['phi'], float) else deviations['phi'],
            'deviation_psi': f"{deviations['psi']:.2f}" if isinstance(deviations['psi'], float) else deviations['psi'],
            'deviation_omega': f"{deviations['omega']:.2f}" if isinstance(deviations['omega'], float) else deviations['omega']
        }
    
    return extracted_data, flagged_data


def main(target_directory=".", 
         all_data_output_filename="extracted_restraint_values.csv",
         flagged_output_filename="flagged_deviations.csv",
         deviation_threshold=10.0):
    """
    Main function to find .out files, extract all data, compare dihedrals,
    and write all extracted data (sorted) and flagged files to separate CSVs.
    """
    if not os.path.isdir(target_directory):
        print(f"Error: Target directory '{target_directory}' not found.", file=sys.stderr)
        sys.exit(1)

    print(f"Starting analysis of .out files in {target_directory}")
    print(f"Non-flagged data will be written to: {all_data_output_filename}")
    print(f"Flagged files (deviation > {deviation_threshold} deg) will be written to: {flagged_output_filename}")
    print("----------------------------------------------------")

    all_out_files = []
    for root, _, files in os.walk(target_directory):
        for file in files:
            if file.endswith(".out"):
                all_out_files.append(os.path.join(root, file))

    if not all_out_files:
        print(f"No .out files found in '{target_directory}'.")
        return

    non_flagged_results = []
    flagged_files_data = []

    for filepath in all_out_files:
        extracted_entry, flagged_entry = extract_and_evaluate_file(filepath, deviation_threshold)
        
        # If an entry is flagged, it only goes to the flagged list.
        # Otherwise, it goes to the main results list.
        if flagged_entry:
            flagged_files_data.append(flagged_entry)
            print(f"Flagged: {flagged_entry['filename']} (Deviations: PHI={flagged_entry['deviation_phi']}, PSI={flagged_entry['deviation_psi']}, OMEGA={flagged_entry['deviation_omega']})")
        elif extracted_entry:
            # Only add to the main list if it was not flagged
            non_flagged_results.append(extracted_entry)

    # --- Write the non-flagged results (sorted by OMEGA) ---
    if non_flagged_results:
        # Sort the data by OMEGA value (increasing order)
        try:
            # Filter out entries where 'omega' might be None before sorting
            sortable_data = [d for d in non_flagged_results if d.get('omega') is not None]
            non_sortable_data = [d for d in non_flagged_results if d.get('omega') is None]

            sortable_data.sort(key=lambda x: x['omega'])
            all_data_sorted = sortable_data + non_sortable_data # Append non-sortable at the end

            print(f"\nWriting {len(all_data_sorted)} non-flagged data entries to {all_data_output_filename}")
            with open(all_data_output_filename, 'w', newline='') as outfile:
                fieldnames = ['filename', 'omega', 'phi', 'psi', 'energy']
                writer = csv.DictWriter(outfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(all_data_sorted)
            print(f"Successfully wrote non-flagged data to {all_data_output_filename}")
        except Exception as e:
            print(f"Error writing non-flagged data to '{all_data_output_filename}': {e}", file=sys.stderr)
    else:
        print("\nNo valid, non-flagged data was found to write.")

    # --- Write flagged files data ---
    if flagged_files_data:
        print(f"\nWriting {len(flagged_files_data)} flagged files to {flagged_output_filename}")
        try:
            with open(flagged_output_filename, mode='w', newline='') as outfile:
                fieldnames = ['filename', 'initial_phi', 'initial_psi', 'initial_omega',
                              'final_phi', 'final_psi', 'final_omega', 'final_energy',
                              'deviation_phi', 'deviation_psi', 'deviation_omega']
                writer = csv.DictWriter(outfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(flagged_files_data)
        except Exception as e:
            print(f"Error writing flagged files to '{flagged_output_filename}': {e}", file=sys.stderr)
    else:
        print("\nNo files exceeded the deviation threshold.")

    print("----------------------------------------------------")
    print("Analysis complete.")

if __name__ == "__main__":
    # Command-line argument usage:
    # python combined_dihedral_analysis.py [target_directory] [all_data_output_file] [flagged_output_file] [deviation_threshold]

    target_dir_arg = "." # Default to current directory
    all_data_output_arg = "extracted_restraint_values.csv" # Default output for all data
    flagged_output_arg = "flagged_deviations.csv" # Default output for flagged data
    deviation_limit = 10.0 # Default deviation threshold

    if len(sys.argv) > 1:
        target_dir_arg = sys.argv[1]
    if len(sys.argv) > 2:
        all_data_output_arg = sys.argv[2]
    if len(sys.argv) > 3:
        flagged_output_arg = sys.argv[3]
    if len(sys.argv) > 4:
        try:
            deviation_limit = float(sys.argv[4])
        except ValueError:
            print(f"Warning: Invalid deviation threshold '{sys.argv[4]}'. Using default 10.0.", file=sys.stderr)

    main(target_dir_arg, all_data_output_arg, flagged_output_arg, deviation_limit)
