import os, sys

import pandas as pd
from loguru import logger

sys.path.append('..')

from pyMultiOmics.constants import *
from pyMultiOmics.mapping import Mapper
from pyMultiOmics.common import set_log_level_info
from pyMultiOmics.analysis import *
from pyMultiOmics.query import *
from pyMultiOmics.pipelines import *


DATA_FOLDER = os.path.abspath(os.path.join('..', 'notebooks', 'test_data', 'fly_data'))
gene_data = pd.read_csv(os.path.join(DATA_FOLDER, 'flyatlas_data.csv'), index_col='Identifier')
gene_design = pd.read_csv(os.path.join(DATA_FOLDER, 'flyatlas_design.csv'), index_col='sample')

compound_data = pd.read_csv(os.path.join(DATA_FOLDER, 'fly_metabolomics_no_dupes.csv'), index_col='Identifier')
compound_design = pd.read_csv(os.path.join(DATA_FOLDER, 'fly_df_design.csv'), index_col='sample')

set_log_level_info()

m = Mapper(DROSOPHILA_MELANOGASTER, metabolic_pathway_only=True, include_related_chebi=True) \
        .set_gene(gene_data, gene_design) \
        .set_compound(compound_data, compound_design) \
        .build()

ap = AnalysisPipeline(m)