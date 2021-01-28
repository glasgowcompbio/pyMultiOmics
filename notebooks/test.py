import os
import sys

import pandas as pd

sys.path.append('..')

from pyMultiOmics.constants import DANIO_RERIO
from pyMultiOmics.mapping import Mapper
from pyMultiOmics.common import set_log_level_info

DATA_FOLDER = os.path.abspath(os.path.join('test_data', 'zebrafish_data'))

gene_data = pd.read_csv(os.path.join(DATA_FOLDER, 'gene_data_combined.csv'), index_col='Identifier')
gene_design = pd.read_csv(os.path.join(DATA_FOLDER, 'gene_design.csv'), index_col='sample')

protein_data = pd.read_csv(os.path.join(DATA_FOLDER, 'protein_data.csv'), index_col='Uniprot')
protein_design = pd.read_csv(os.path.join(DATA_FOLDER, 'protein_design.csv'), index_col='sample')

compound_data = pd.read_csv(os.path.join(DATA_FOLDER, 'compound_data_kegg.csv'), index_col='Identifier')
compound_design = pd.read_csv(os.path.join(DATA_FOLDER, 'compound_design.csv'), index_col='sample')

set_log_level_info()

m = Mapper(DANIO_RERIO, metabolic_pathway_only=True) \
    .set_gene(gene_data, gene_design) \
    .set_protein(protein_data, protein_design) \
    .set_compound(compound_data, compound_design) \
    .build()
