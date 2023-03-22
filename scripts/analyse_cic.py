import argparse
import os
import subprocess
import sys
import yaml
import papermill as pm

sys.path.append('.')
sys.path.append('..')

from pyMultiOmics.common import create_if_not_exist


def load_config_file(config_file):
    with open(config_file) as f:
        config = yaml.safe_load(f)
    return config


def main(args):
    input_excel = os.path.abspath(args.input_excel)
    output_dir = os.path.abspath(args.output_dir)

    print('Input Excel:', input_excel)
    print('Output Directory:', output_dir)

    # Check if template file exists in parent directory
    template_parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    template_file = os.path.abspath(
        os.path.join(template_parent_dir, 'templates', 'cic_template.ipynb'))

    if os.path.isfile(template_file):
        print('Template Notebook:', template_file)
    else:
        # Check if template file exists in current directory
        template_file = os.path.abspath(os.path.join('.', 'templates', 'cic_template.ipynb'))

        if os.path.isfile(template_file):
            print('Template Notebook:', template_file)
        else:
            print('Error: template notebook not found')
            return

    output_notebook = os.path.abspath(os.path.join(output_dir, 'cic_analysis.ipynb'))
    output_pdf = os.path.abspath(os.path.join(output_dir, 'cic_analysis.pdf'))

    print('Output Notebook:', output_notebook)
    print('Output PDF:', output_pdf)

    config_file = os.path.abspath(args.config_file)
    processing_params = load_config_file(config_file)

    print('Pre-processing parameters:')
    print(yaml.dump(processing_params, indent=4))

    params = {}
    params['file_name'] = input_excel
    params.update(processing_params)

    create_if_not_exist(output_dir)
    _ = pm.execute_notebook(template_file, output_notebook, report_mode=True, parameters=params)

    command = ['jupyter', 'nbconvert', '--TemplateExporter.exclude_input=True', '--output',
               output_pdf,
               '--to', 'pdf', output_notebook]
    result = subprocess.run(command, capture_output=True, text=True)

    if result.returncode == 0:
        print(f'\nPDF output saved to {output_pdf}\n')
    else:
        print(f'\nFailed to create PDF: {result.stderr}\n')


if __name__ == '__main__':
    print('***********************************')
    print('* CIC Analysis                    *')
    print('***********************************')
    print()

    parser = argparse.ArgumentParser(
        description='Runs a notebook-based analysis on input data and saves the results as a PDF.')
    parser.add_argument('--input_excel', required=True, help='Path to the input Excel file.')
    parser.add_argument('--output_dir', required=True, help='Path to the output directory.')
    parser.add_argument('--config_file', required=True, help='Path to the YAML configuration file.')
    args = parser.parse_args()
    main(args)
