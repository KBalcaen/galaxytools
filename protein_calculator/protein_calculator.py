################################################################################################
# This tool was developed by Kevin Balcaen for usage by the VIB Protein Core staff.            #
# Grateful acknowledgment to Robbe Fonteyn for the error checking on the code and proofreading.#
# The tool was updated to a Galaxy tool wrapper for usegalaxy.be by Boris Depoortere           #
################################################################################################

import argparse
import pandas as pd
from json import dumps
from plotly import utils, express as px
from Calculate_protein_properties import (get_isoelectric_point, calculate_dn_dc, calculate_total_masses,
                                          calculate_extinction_coefficient, get_net_charge)
from Sequence_functions import normalize_sequence, check_protein_sequence, format_sequence, letter_count
from references import amino_acid_data
from jinja2 import Environment, FileSystemLoader
import base64
import os

VERSION = '1.0.3'

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process a protein sequence or FASTA file.')

    # Definitie van de argumenten
    parser.add_argument('--sequence', type=str, required=True,
                        help='Protein sequence or path to a FASTA file')

    parser.add_argument('--name', type=str, required=False, help='Name of the protein (optional if FASTA has headers)')

    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {VERSION}',
                        help='Show the version of the tool and exit')

    return parser.parse_args()

def read_sequences_from_fasta(fasta_path):
    """Leest meerdere sequenties uit een FASTA-bestand en retourneert een lijst van (naam, sequentie)-tuples."""
    sequences = []
    with open(fasta_path, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header and seq_lines:
                    sequences.append((header, ''.join(seq_lines).upper()))
                header = line[1:].strip().split()[0]  # Neem eerste woord als naam
                seq_lines = []
            else:
                seq_lines.append(line)
        if header and seq_lines:
            sequences.append((header, ''.join(seq_lines).upper()))
    return sequences

def resolve_sequences(sequence_input):
    """Geeft een lijst van (name, sequence)-tuples terug, ongeacht of het een string of FASTA-bestand is."""
    if os.path.isfile(sequence_input):
        return read_sequences_from_fasta(sequence_input)
    else:
        return [("input_sequence", sequence_input.upper())]

def main():
    args = parse_arguments()
    sequences = resolve_sequences(args.sequence)

    for i, (fasta_name, sequence) in enumerate(sequences):
        if args.name and len(sequences) == 1:
            name = args.name
        else:
            name = fasta_name

        print(f"Processing: {name}")

        sequence = normalize_sequence(sequence)
        error_message = check_protein_sequence(sequence)
        if error_message:
            print(f"Error in {name}: {error_message}")
            continue

        sequence_formatted = format_sequence(sequence, show_residue_number=True, line_length=55)

        frequency = letter_count(sequence)
        amino_acid_composition = []
        total_monoisotopic_mass = 18.0152
        total_average_mass = 18.0152
        total_extinction_coefficient_cystines = 0
        total_extinction_coefficient_reduced = 0
        total_count = sum(frequency.values())

        for amino_acid in sorted(frequency.keys()):
            count = frequency[amino_acid]
            monoisotopic_mass, average_mass = calculate_total_masses(amino_acid, count)
            total_monoisotopic_mass += monoisotopic_mass
            total_average_mass += average_mass
            total_extinction_coefficient_cystines += calculate_extinction_coefficient(amino_acid, count, True)
            total_extinction_coefficient_reduced += calculate_extinction_coefficient(amino_acid, count, False)
            data = amino_acid_data.get(amino_acid)
            percentage = (count / total_count) * 100
            amino_acid_composition.append({
                'amino_acid': amino_acid,
                'long_name': data.get('Long', ''),
                'mono_weight': round(data.get('Monoisotopic Weight (Da)', ''), 2),
                'avg_weight': round(data.get('Average Weight (Da)', ''), 2),
                'count': count,
                'percentage': f"{percentage:.2f}%"
            })

        number_info = {'number': total_count}

        if total_monoisotopic_mass > 9999:
            temp_str = f"{total_monoisotopic_mass:,.2f}"
            rounded_total_monoisotopic_mass = temp_str.replace(",", " ")
        else:
            rounded_total_monoisotopic_mass = f"{total_monoisotopic_mass:.2f}"

        if total_average_mass > 9999:
            temp_str = f"{total_average_mass:,.2f}"
            rounded_total_average_mass = temp_str.replace(",", " ")
        else:
            rounded_total_average_mass = f"{total_average_mass:.2f}"

        molecular_weight_info = {
            'mono_weight': rounded_total_monoisotopic_mass,
            'avg_weight': rounded_total_average_mass
        }

        molar_absorbance_info = {
            'extinction_coefficient_cystines': total_extinction_coefficient_cystines,
            'extinction_coefficient_reduced': total_extinction_coefficient_reduced,
            'absorbance_mono_cystines': round(total_extinction_coefficient_cystines / total_monoisotopic_mass, 4) if total_monoisotopic_mass != 0 else 0,
            'absorbance_avg_cystines': round(total_extinction_coefficient_cystines / total_average_mass, 4) if total_average_mass != 0 else 0,
            'absorbance_mono_reduced': round(total_extinction_coefficient_reduced / total_monoisotopic_mass, 4) if total_monoisotopic_mass != 0 else 0,
            'absorbance_avg_reduced': round(total_extinction_coefficient_reduced / total_average_mass, 4) if total_average_mass != 0 else 0
        }

        pI = round(get_isoelectric_point(sequence), 2)
        net_charge_at_different_pH = []
        pH_values = []
        net_charges = []
        charges = get_net_charge(sequence)

        for pH, net_charge in charges.items():
            net_charge_at_different_pH.append({
                'pH': pH,
                'net_charge': round(net_charge, 2)
            })
            pH_values.append(pH)
            net_charges.append(net_charge)

        df = pd.DataFrame({"pH": pH_values, "Net Charges": net_charges})
        colors = ['#3ABBBA', '#5A2A82', '#FF681E']
        titration_curve = px.line(df, x="pH", y="Net Charges", color_discrete_sequence=colors)
        titration_curve.update_layout(title_text="Titration Curve", title_x=0.5)

        file_prefix = name.replace(" ", "_")
        html_file = f"{file_prefix}_report.html"
        png_file = f"{file_prefix}_plot.png"
        json_file = f"{file_prefix}_titration.json"

        titration_curve.write_html(html_file)
        titration_curve.write_image(png_file, format="png")

        with open(png_file, "rb") as png_file_obj:
            titration_png_base64 = base64.b64encode(png_file_obj.read()).decode("utf-8")

        titration_json = dumps(titration_curve, cls=utils.PlotlyJSONEncoder)
        dn_dc_value = round(calculate_dn_dc(sequence, amino_acid_data), 6)

        env = Environment(loader=FileSystemLoader('templates'))
        template = env.get_template('results.html')
        data = {
            "name": name,
            "sequence": sequence_formatted,
            "sequence2": sequence,
            "amino_acid_composition": amino_acid_composition,
            "number_info": number_info,
            "molecular_weight_info": molecular_weight_info,
            "molar_absorbance_info": molar_absorbance_info,
            "pI": pI,
            "net_charge_at_different_pH": net_charge_at_different_pH,
            "titration_json": titration_json,
            "titration_png_base64": titration_png_base64,
            "titration_curve": titration_curve,
            "dn_dc_value": dn_dc_value
        }

        output = template.render(data)
        with open(html_file, 'w', encoding="utf-8") as f:
            f.write(output)

        # Print Results
        print(f'Name: {name}')
        print(f'Sequence: {sequence_formatted}')
        print(f'Amino Acid Composition:')
        print(
            f'{"Amino Acid (Short)":<20} {"Amino Acid (Long)":<20} {"Monoisotopic Weight (Da)":<25} {"Average Weight (Da)":<25} {"# Counts":<10} {"% of Total":<10}')
        for aa in amino_acid_composition:
            print(
                f'{aa["amino_acid"]:<20} {aa["long_name"]:<20} {aa["mono_weight"]:<25} {aa["avg_weight"]:<25} {aa["count"]:<10} {aa["percentage"]:<10}')
        print(f'Molecular Weight Info: {molecular_weight_info}')
        print(f'Molar Absorbance Info: {molar_absorbance_info}')
        print(f'pI: {pI}')
        print(f'Net Charge at Different pH: {net_charge_at_different_pH}')
        print(f'Titration Curve JSON: {titration_json}')
        print(f'dn/dc Value: {dn_dc_value}')

        print(f"HTML rendered and saved to {html_file}")

if __name__ == '__main__':
    main()
