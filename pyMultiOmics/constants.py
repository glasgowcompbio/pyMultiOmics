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

GENES = 'genes'
PROTEINS = 'proteins'
COMPOUNDS = 'compounds'
REACTIONS = 'reactions'
PATHWAYS = 'pathways'
GENES_TO_PROTEINS = 'gene_proteins'
PROTEINS_TO_REACTIONS = 'protein_reactions'
COMPOUNDS_TO_REACTIONS = 'compound_reactions'
REACTIONS_TO_PATHWAYS = 'reaction_pathways'
MULTI_OMICS = 'multi_omics'

PKS = {
    GENES: GENE_PK,
    PROTEINS: PROTEIN_PK,
    COMPOUNDS: COMPOUND_PK,
    REACTIONS: REACTION_PK,
    PATHWAYS: PATHWAY_PK
}

IDS = {
    GENES: GENE_ID,
    PROTEINS: PROTEIN_ID,
    COMPOUNDS: COMPOUND_ID,
    REACTIONS: REACTION_ID,
    PATHWAYS: PATHWAY_ID
}

COMPOUND_DATABASE_KEGG = 'KEGG'
COMPOUND_DATABASE_CHEBI = 'ChEBI'

IDENTIFIER_COL = 'Identifier'
PVALUE_COL_PREFIX = 'pvalue_'
PADJ_COL_PREFIX = 'padj_' # adjusted pvalue column, unused?
FC_COL_PREFIX = 'FC_'
SAMPLE_COL = 'sample'
GROUP_COL = 'group'
FACTOR_COL = 'factor'
DEFAULT_GROUP_NAME = 'default'

PIMP_PEAK_ID_COL = 'Peak id'
SMALL = 1E-6

# Constants used in the Inference page

INFERENCE_LOADED = 'loaded'
INFERENCE_T_TEST = 't_test'
INFERENCE_CORRELATION = 'correlation'
INFERENCE_PCA = 'pca'
INFERENCE_PALS = 'pals'
INFERENCE_ORA = 'ora'
INFERENCE_GSEA = 'gsea'
INFERENCE_DESEQ = 'deseq'
INFERENCE_LIMMA = 'limma'
INFERENCE_REACTOME = 'reactome'

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
