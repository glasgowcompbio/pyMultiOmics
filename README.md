# pyMultiOmics

pyMultiOmics is a Python package for multi-omics data integration and analysis. 
It uses the Reactome database to map entities (genes, transcripts, proteins, compounds) to
their reactions and pathways. The results is then shown as a network graph. Various analyses such as
differential analysis, pathway activity analysis can be performed on this network graph, with the results 
overlaid on the graph.

### Installation

Simply run:
```
pip install pyMultiOmics
```

### Usage

Example basic usage is shown below:
```
from pyWebOmics.mapping import Mapper

species_name = DANIO_RERIO
m = Mapper(species_name) \
        .set_gene(gene_data, gene_design) \
        .set_protein(protein_data, protein_design) \
        .set_compound(compound_data, compound_design) \
        .build()
```

`m` contains a mapper object, which can be interogated to obtain the data integration results. 
Please refer to this notebook for a demo.