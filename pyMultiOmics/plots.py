import json
import shutil
from urllib.parse import quote

import pandas as pd
import requests
from loguru import logger


def download_reactome_diagram(stId, out_file, species=None, expression_dict=None):
    """
    Download a reaction or pathway diagram from Reactome.
    If provided, the downloaded diagram can be overlaid with the expression data.

    Args:
        stId: the Reactome's stable ID of the reaction or pathway to download, e.g. 'R-HSA-70921'
        out_file: the name of output file containing the image
        species: the name of the species. Useful constants are defined in pyMultiOmics.constants.
        expression_dict: a dictionary of expression values to colour in the diagram, e.g.
                        {
                            'C00135': 1.0,           # kegg
                            'C01342': 2.0,           # kegg
                            '30817': 3.0,            # ChEBI
                            'HAL': 4.0,              # gene name
                            'ENSG00000172508': 5.0,  # ensembl
                            'Q96NU7': 6.0            # uniprot
                        }

    Returns: None

    """

    # if expression data is available, send it to Reactome Analysis service and get back a token.
    reactome_token = None
    if expression_dict is not None:
        assert species is not None
        expression_df = pd.DataFrame.from_dict(expression_dict.items()).set_index(0)
        expression_df = expression_df.rename(columns={1: 'value'})
        expression_df.index.name = '#id'

        # Submit expression data to Reactome Analysis service
        expression_data = expression_df.to_csv(sep='\t', header=True, index_label='#id',
                                               float_format='%.15f')
        encoded_species = quote(species)

        status_code, json_response = send_reactome_expression_data(expression_data,
                                                                   encoded_species)
        if status_code == 200:
            pathways_df, reactome_url, reactome_token = parse_reactome_json(json_response)

    # Generate diagram by calling Reactome Exporter service.
    # The API is documented here:
    # https://reactome.org/ContentService/#/exporter/diagramImageUsingGET_1
    image_url = 'https://reactome.org/ContentService/exporter/diagram/%s.png?quality=8' \
                '&diagramProfile=standard&analysisProfile=strosobar' % stId

    # If we have an analysis token, then include it for diagram export too.
    if reactome_token is not None:
        image_url += '&token=%s&resource=TOTAL&expColumn=0' % reactome_token

    # Save image_url to file
    r = requests.get(image_url, stream=True)
    if r.status_code == 200:
        with open(out_file, 'wb') as f:
            r.raw.decode_content = True
            shutil.copyfileobj(r.raw, f)
            logger.debug('Image saved to %s' % out_file)


def send_reactome_expression_data(data, encoded_species):
    """
    Send expression data to Reactome Analysis service.
    API doc: https://reactome.org/AnalysisService/#/identifiers/getPostTextUsingPOST
    The data format is described here: https://reactome.org/dev/analysis

    Args:
        data: the data to send in tsv format
        encoded_species: the name of the species

    Returns: HTTP status code and a JSON response that can be parsed by parse_reactome_json()

    """

    url = 'https://reactome.org/AnalysisService/identifiers/?interactors=false&species=' + \
          encoded_species + '&sortBy=ENTITIES_PVALUE&order=ASC&resource=TOTAL&pValue=1&includeDisease=true'
    logger.debug('POSTing expression data to Reactome Analysis Service: ' + url)

    response = requests.post(url, headers={'Content-Type': 'text/plain'},
                             data=data.encode('utf-8'))
    logger.debug('Received HTTP status code: %d' % response.status_code)

    status_code = response.status_code
    if status_code == 200:
        json_response = json.loads(response.text)
    else:
        json_response = None
    return status_code, json_response


def parse_reactome_json(json_response):
    """
    Parse the JSON returned by Reactome Analysis service
    Args:
        json_response: the JSON response

    Returns: a pathway dataframe containing p-values of pathways from over-representation analysis,
             the reactome URL used to submit the data to, and an analysis token

    """
    # see https://reactome.org/userguide/analysis for results explanation
    token = json_response['summary']['token']
    pathways = json_response['pathways']

    reactome_url = 'https://reactome.org/PathwayBrowser/#DTAB=AN&ANALYSIS=' + token
    logger.debug('Received expression analysis token: ' + token)

    # https://stackoverflow.com/questions/6027558/flatten-nested-dictionaries-compressing-keys
    pathways_df = pd.json_normalize(pathways, sep='_')
    return pathways_df, reactome_url, token
