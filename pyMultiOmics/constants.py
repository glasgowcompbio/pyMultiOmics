import os


def get_data_path(data_dir, file_name):
    return os.path.abspath(os.path.join(data_dir, '%s' % file_name))


BASE_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
EXTERNAL_COMPOUND_NAMES = get_data_path(DATA_DIR, 'compound_names.p')
EXTERNAL_KEGG_TO_CHEBI = get_data_path(DATA_DIR, 'kegg_to_chebi.p')
EXTERNAL_GENE_NAMES = get_data_path(DATA_DIR, 'gene_names.p')
EXTERNAL_GO_DATA = get_data_path(DATA_DIR, 'go_data.p')

CHEBI_BFS_RELATION_DICT = get_data_path(DATA_DIR, 'chebi_bfs_relation_dict.pkl')
CHEBI_RELATION_TSV = get_data_path(DATA_DIR, 'relation.tsv')

ARABIDOPSIS_THALIANA = 'Arabidopsis thaliana'
BOS_TAURUS = 'Bos taurus'
CAENORHABDITIS_ELEGANS = 'Caenorhabditis elegans'
CANIS_LUPUS_FAMILIARIS = 'Canis lupus familiaris'
DANIO_RERIO = 'Danio rerio'
DICTYOSTELIUM_DISCOIDEUM = 'Dictyostelium discoideum'
DROSOPHILA_MELANOGASTER = 'Drosophila melanogaster'
GALLUS_GALLUS = 'Gallus gallus'
HOMO_SAPIENS = 'Homo sapiens'
MUS_MUSCULUS = 'Mus musculus'
ORYZA_SATIVA = 'Oryza sativa'
RATTUS_NORVEGICUS = 'Rattus norvegicus'
SACCHAROMYCES_CEREVISIAE = 'Saccharomyces cerevisiae'
SUS_SCROFA = 'Sus scrofa'

NA = '-'
ALL = 'ALL'

GENE_PK = 'gene_pk'
PROTEIN_PK = 'protein_pk'
COMPOUND_PK = 'compound_pk'
REACTION_PK = 'reaction_pk'
PATHWAY_PK = 'pathway_pk'

GENE_ID = 'gene_id'
PROTEIN_ID = 'protein_id'
COMPOUND_ID = 'compound_id'
REACTION_ID = 'reaction_id'
PATHWAY_ID = 'pathway_id'

GENOMICS = 0
TRANSCRIPTOMICS = 1
PROTEOMICS = 2
METABOLOMICS = 3
REACTIONS = 4
PATHWAYS = 5
GENES_TO_PROTEINS = 6
PROTEINS_TO_REACTIONS = 7
COMPOUNDS_TO_REACTIONS = 8
REACTIONS_TO_PATHWAYS = 9
MULTI_OMICS = 10

# aliases
GENES = GENOMICS
TRANSCRIPTS = TRANSCRIPTOMICS
PROTEINS = PROTEOMICS
COMPOUNDS = METABOLOMICS

DataType = (
    (GENOMICS, 'Genomics'),
    (TRANSCRIPTOMICS, 'Transcriptomics'),
    (PROTEOMICS, 'Proteomics'),
    (METABOLOMICS, 'Metabolomics'),
    (REACTIONS, 'Reactions'),
    (PATHWAYS, 'Pathways')
)
DataTypeDict = dict(DataType)

MAPPING = {
    GENOMICS: 'genes',
    PROTEOMICS: 'proteins',
    METABOLOMICS: 'compounds',
    REACTIONS: 'reactions',
    PATHWAYS: 'pathways',
    GENES_TO_PROTEINS: 'gene_proteins',
    PROTEINS_TO_REACTIONS: 'protein_reactions',
    COMPOUNDS_TO_REACTIONS: 'compound_reactions',
    REACTIONS_TO_PATHWAYS: 'reaction_pathways'
}

PKS = {
    GENOMICS: GENE_PK,
    PROTEOMICS: PROTEIN_PK,
    METABOLOMICS: COMPOUND_PK,
    REACTIONS: REACTION_PK,
    PATHWAYS: PATHWAY_PK
}

IDS = {
    GENOMICS: GENE_ID,
    PROTEOMICS: PROTEIN_ID,
    METABOLOMICS: COMPOUND_ID,
    REACTIONS: REACTION_ID,
    PATHWAYS: PATHWAY_ID
}

COMPOUND_DATABASE_KEGG = 'KEGG'
COMPOUND_DATABASE_CHEBI = 'ChEBI'

IDENTIFIER_COL = 'Identifier'
PADJ_COL_PREFIX = 'padj_'
FC_COL_PREFIX = 'FC_'
SAMPLE_COL = 'sample'
GROUP_COL = 'group'
FACTOR_COL = 'factor'
DEFAULT_GROUP_NAME = 'default'

PIMP_PEAK_ID_COL = 'Peak id'
SMALL = 1E-6

# Constants used in the Inference page

INFERENCE_LOADED = 0
INFERENCE_T_TEST = 1
INFERENCE_CORRELATION = 2
INFERENCE_PCA = 3
INFERENCE_PALS = 4
INFERENCE_ORA = 5
INFERENCE_GSEA = 6
INFERENCE_DESEQ = 7
INFERENCE_LIMMA = 8
INFERENCE_REACTOME = 9

T_TEST_THRESHOLD = 0.05

QUERY_ENTITY_ID = 'entity_id'
QUERY_OBSERVED = 'observed'
QUERY_DATA_TYPE = 'data_type'
QUERY_NODE_ID = 'node_id'
QUERY_DISPLAY_NAME = 'display_name'
QUERY_SOURCE_ID = 'source_id'

TRUNCATE_LIMIT = 10000

MEASUREMENT_DF_LABEL = 'measurement dataframe'
DESIGN_DF_LABEL = 'sample metadata dataframe'
