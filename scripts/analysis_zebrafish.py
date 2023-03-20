import os, sys

import pandas as pd
from loguru import logger

from pyMultiOmics.base import SingleOmicsData, MultiOmicsData

sys.path.append('..')

from pyMultiOmics.constants import *
from pyMultiOmics.mapping import Mapper
from pyMultiOmics.common import set_log_level_info
from pyMultiOmics.analysis import *
from pyMultiOmics.query import *
from pyMultiOmics.pipelines import *


DATA_FOLDER = os.path.abspath(os.path.join('..', 'notebooks', 'zebrafish_data'))
gene_data = pd.read_csv(os.path.join(DATA_FOLDER, 'gene_data_combined.csv'), index_col='Identifier')
gene_design = pd.read_csv(os.path.join(DATA_FOLDER, 'gene_design.csv'), index_col='sample')
protein_data = pd.read_csv(os.path.join(DATA_FOLDER, 'protein_data.csv'), index_col='Uniprot')
protein_design = pd.read_csv(os.path.join(DATA_FOLDER, 'protein_design.csv'), index_col='sample')
compound_data = pd.read_csv(os.path.join(DATA_FOLDER, 'compound_data_kegg.csv'), index_col='Identifier')
compound_design = pd.read_csv(os.path.join(DATA_FOLDER, 'compound_design.csv'), index_col='sample')

set_log_level_info()

gene_so = SingleOmicsData(GENES, gene_data, gene_design)
prot_so = SingleOmicsData(GENES, protein_data, protein_design)
comp_so = SingleOmicsData(GENES, compound_data, compound_design)

mo = MultiOmicsData()
mo.add_data([gene_so, prot_so, comp_so])

m = Mapper(mo, DANIO_RERIO, metabolic_pathway_only=True) \
        .build()

ap = AnalysisPipeline(mo, m)

# method = INFERENCE_DESEQ
method = INFERENCE_T_TEST
ap.run_de(method, GENES, 'Distal', 'Proximal')
ap.run_de(method, GENES, 'Distal', 'Middle')
ap.run_de(method, GENES, 'Proximal', 'Middle')

method = INFERENCE_T_TEST
ap.run_de(method, PROTEINS, 'Distal', 'Proximal')
ap.run_de(method, PROTEINS, 'Distal', 'Middle')
ap.run_de(method, PROTEINS, 'Proximal', 'Middle')
ap.run_de(method, COMPOUNDS, 'Distal', 'Proximal')
ap.run_de(method, COMPOUNDS, 'Distal', 'Middle')
ap.run_de(method, COMPOUNDS, 'Proximal', 'Middle')

# Retrieve a single node
node_id = '15366'
res = QueryBuilder(ap) \
        .add(Entity(node_id)) \
        .run()
logger.info(res)

# Retrieve multiple nodes
node_id = ['15366', 'ENSDARG00000037781', 'F1QAA7']
res = QueryBuilder(ap) \
        .add(Entity(node_id)) \
        .run()
logger.info(res)

# Retrieve nodes connected to a single node
query_id = 'F1QAA7'
res = QueryBuilder(ap) \
        .add(Entity(query_id)) \
        .add(Connected()) \
        .run()
logger.info(res)

# Retrieve top-10 significantly changing genes
case = 'Distal'
control = 'Proximal'
pval = 0.05
fc_lte = -2
fc_gte = 2
N = 20
res = QueryBuilder(ap) \
        .add(Select(GENES)) \
        .add(SignificantDE(case, control, pval, fc_lte=fc_lte, fc_gte=fc_gte, N=N)) \
        .run()
logger.info(res)

# Find the compounds that are connected to the DE genes above

res = QueryBuilder(ap) \
        .add(Select(GENES)) \
        .add(SignificantDE(case, control, pval, fc_lte=fc_lte, fc_gte=fc_gte, N=N)) \
        .add(Connected(data_type=COMPOUNDS)) \
        .run()
logger.info(res)

# Plot some heatmap using Plotly

res = QueryBuilder(ap) \
        .add(Select(GENES)) \
        .add(SignificantDE(case, control, pval, fc_lte=fc_lte, fc_gte=fc_gte, N=N)) \
        .run()
logger.info(res)

data_type = GENES
analysis = ap.get_de_analysis(data_type, case, control)
wi = analysis.wi
data_df, design_df = wi.data_df, wi.design_df

case_group = design_df[design_df['group'] == case].index.values.tolist()
control_group = design_df[design_df['group'] == control].index.values.tolist()
selection = case_group + control_group

heatmap_df = wi.data_df.loc[res.index.values]
heatmap_df = heatmap_df[selection]
heatmap_df

from plotly import express as px
px.imshow(wi._normalise_df(heatmap_df, log=True))