import os

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
EXTERNAL_COMPOUND_NAMES = os.path.join(BASE_DIR, 'data', 'compound_names.p')
EXTERNAL_KEGG_TO_CHEBI = os.path.join(BASE_DIR, 'data', 'kegg_to_chebi.p')
EXTERNAL_GENE_NAMES = os.path.join(BASE_DIR, 'data', 'gene_names.p')
EXTERNAL_GO_DATA = os.path.join(BASE_DIR, 'data', 'go_data.p')

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

GENOMICS, TRANSCRIPTOMICS, PROTEOMICS, METABOLOMICS, REACTIONS, PATHWAYS, MULTI_OMICS = range(0, 7)
DataType = (
    (GENOMICS, 'Genomics'),
    (TRANSCRIPTOMICS, 'Transcriptomics'),
    (PROTEOMICS, 'Proteomics'),
    (METABOLOMICS, 'Metabolomics'),
    (REACTIONS, 'Reactions'),
    (PATHWAYS, 'Pathways')
)
DataTypeDict = dict(DataType)