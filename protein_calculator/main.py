################################################################################################
# This tool was developed by Kevin Balcaen for usage by the VIB Protein Core staff.            #
# Grateful acknowledgment to Robbe Fonteyn for the error checking on the code and proofreading.#
################################################################################################

from json import dumps

import pandas as pd
from flask import Flask, render_template, request, session, redirect, url_for
from plotly import utils, express as px

from Calculate_protein_properties import (get_isoelectric_point, calculate_dn_dc, calculate_total_masses,
                                          calculate_extinction_coefficient, get_net_charge)
from Sequence_functions import normalize_sequence, check_protein_sequence, format_sequence, letter_count
from references import pK_values, amino_acid_data

# Flask-application initialisation
app = Flask(__name__)


# Route for homepage with form
@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        name = request.form['name']
        sequence = request.form['sequence'].upper()

        # Normalise sequence to remove spaces and enters and check validity of sequence
        sequence = normalize_sequence(sequence)

        # Store necessary information in session
        session['name'] = name
        session['sequence'] = sequence

        # Redirect to the results page
        return redirect(url_for('results'))

    return render_template('index.html')


# Results route
@app.route('/results', methods=["GET"])
def results():
    name = session.get('name')
    sequence = session.get('sequence')

    error_message = check_protein_sequence(sequence)

    # Error Checking
    if error_message:
        session["error_message"] = error_message
        return redirect(url_for('calculator_error'))

    sequence_formatted = format_sequence(sequence, show_residue_number=True, line_length=55)

    # Perform calculations
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

    number_info = {
        'number': total_count
    }

    # Clean up of total_monoisotopic_mass and total_average_mass
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
        'absorbance_mono_cystines': round(total_extinction_coefficient_cystines / total_monoisotopic_mass,
                                          4) if total_monoisotopic_mass != 0 else 0,
        'absorbance_avg_cystines': round(total_extinction_coefficient_cystines / total_average_mass,
                                         4) if total_average_mass != 0 else 0,
        'absorbance_mono_reduced': round(total_extinction_coefficient_reduced / total_monoisotopic_mass,
                                         4) if total_monoisotopic_mass != 0 else 0,
        'absorbance_avg_reduced': round(total_extinction_coefficient_reduced / total_average_mass,
                                        4) if total_average_mass != 0 else 0
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

    # Plotly Titration Curve
    df = pd.DataFrame({"pH": pH_values, "Net Charges": net_charges})
    colors = ['#3ABBBA', '#5A2A82', '#FF681E']
    titration_curve = px.line(df, x="pH", y="Net Charges", color_discrete_sequence=colors)
    titration_curve.update_layout(title_text="Titration Curve", title_x=0.5)

    titration_json = dumps(titration_curve, cls=utils.PlotlyJSONEncoder)

    dn_dc_value = round(calculate_dn_dc(sequence, amino_acid_data), 6)

    # Put right variables in session
    session["name"] = name
    session['sequence'] = sequence
    session['molecular_weight_info'] = molecular_weight_info
    session['molar_absorbance_info'] = molar_absorbance_info
    session['dn_dc_value'] = dn_dc_value
    session['pI'] = pI
    session['net_charge_at_different_pH'] = net_charge_at_different_pH

    return render_template('results.html',
                           name=name,
                           sequence=sequence_formatted,
                           sequence2=sequence,
                           amino_acid_composition=amino_acid_composition,
                           number_info=number_info,
                           molecular_weight_info=molecular_weight_info,
                           molar_absorbance_info=molar_absorbance_info,
                           pI=pI,
                           net_charge_at_different_pH=net_charge_at_different_pH,
                           titration_json=titration_json,
                           dn_dc_value=dn_dc_value)


# Export view page
@app.route('/export_view', methods=["GET"])
def export_data():
    # Get session variables
    name = session.get('name')
    sequence = session.get('sequence')
    molecular_weight_info = session.get('molecular_weight_info')
    molar_absorbance_info = session.get('molar_absorbance_info')
    dn_dc_value = session.get('dn_dc_value')
    pI = session.get('pI')
    net_charge_at_different_pH = session.get('net_charge_at_different_pH')

    return render_template("export_view.html",
                           name=name,
                           sequence=sequence,
                           molecular_weight_info=molecular_weight_info,
                           molar_absorbance_info=molar_absorbance_info,
                           dn_dc_value=dn_dc_value,
                           pI=pI,
                           net_charge_at_different_pH=net_charge_at_different_pH)


@app.route('/calculator_error', methods=["GET"])
def calculator_error():
    error_message = session.get("error_message")
    return render_template('protein_calc_error.html', error_message=error_message)


# Start the Flask-application
if __name__ == '__main__':
    app.secret_key = "protein core"
    app.config["SESSION_TYPE"] = "filesystem"

    app.run(debug=False, port=8000, host="localhost")
