from collections import defaultdict

import numpy as np
import pandas as pd
import requests
import shutil

def download_reactome_plot(stId, out_file, species=None, expression_dict=None):

    # Calls https://reactome.org/ContentService/#/exporter/diagramImageUsingGET_1
    image_url = 'https://reactome.org/ContentService/exporter/diagram/%s.png?quality=8' \
                '&diagramProfile=standard&analysisProfile=strosobar' % stId

    # if expression data is available
    expression_df = pd.DataFrame.from_dict(expression_dict.items()).set_index(0)
    expression_df = expression_df.rename(columns={1: 'log_FC'})
    expression_df.index.name = '#id'

    expression_data = expression_df.to_csv(sep='\t', header=True, index_label='#id', float_format='%.15f')
    encoded_species = quote(species)

    r = requests.get(image_url, stream=True)
    if r.status_code == 200:
        with open(out_file, 'wb') as f:
            r.raw.decode_content = True
            shutil.copyfileobj(r.raw, f)

            # if token is not None:
    #     image_url += '&token=%s&resource=TOTAL&expColumn=0' % token


# def send_expression_data(ds, case, control, species):
#     # send expression data to reactome for diagram exporter
#     int_df = ds.change_zero_peak_ints()
#     annot_df = ds.get_annotations()
#     design = ds.get_experimental_design()
#
#     case_cols = design['groups'][case]
#     control_cols = design['groups'][control]
#
#     df = annot_df.join(int_df).set_index('entity_id').sort_index()
#     case_df = np.log2(df[case_cols])
#     control_df = np.log2(df[control_cols])
#     lfcs = np.mean(case_df, axis=1) - np.mean(control_df, axis=1)
#
#     # for each compound, compute the average fold changes if there are multiple values
#     # TODO: we need a better way to do this
#     temp = defaultdict(list)
#     for idx, lfc in lfcs.iteritems():
#         temp[idx].append(lfc)
#
#     fold_changes = {}
#     for idx in temp:
#         mean = np.mean(temp[idx])
#         # logger.debug(idx, mean)
#         fold_changes[idx] = mean
#
#     # create expression dataframe to send to reactome
#     expression_df = pd.DataFrame.from_dict(fold_changes.items()).set_index(0)
#     expression_df = expression_df.rename(columns={1: 'log_FC'})
#     expression_df.index.name = '#id'
#
#     expression_data = expression_df.to_csv(sep='\t', header=True, index_label='#id', float_format='%.15f')
#     encoded_species = quote(species)
#
#     status_code, json_response = send_reactome_expression_data(expression_data, encoded_species)
#     if status_code == 200:
#         pathways_df, reactome_url, reactome_token = parse_reactome_json(json_response)
#         return reactome_token
#     else:
#         st.warning('Failed to submit expression data to Reactome.org (status_code=%d)' % status_code)
#         return None
#
# def show_reactome_diagram(df, json_response, pw_name, stId, token):
#     # st.subheader(pw)
#     label = '%s: %s' % (stId, pw_name)
#     info_url = 'https://reactome.org/content/detail/%s' % stId
#     header_markdown = '### %s [[info]](%s)' % (label, info_url)
#     if token is not None:
#         viewer_url = 'https://reactome.org/PathwayBrowser/#/%s&DTAB=AN&ANALYSIS=%s' % (stId, token)
#         header_markdown += ' [[viewer]](%s)' % viewer_url
#     st.write(header_markdown)
#
#     row = df.loc[stId]
#     st.write(row)
#     for summation in json_response['summation']:
#         summary = summation['text']
#         st.write(summary)
#
#     image_url = 'https://reactome.org/ContentService/exporter/diagram/%s.png?quality=8' \
#                 '&diagramProfile=standard&analysisProfile=strosobar' % stId
#     if token is not None:
#         image_url += '&token=%s&resource=TOTAL&expColumn=0' % token
#     logger.debug('image_url = %s' % image_url)
#     st.image(image_url, use_column_width=True)