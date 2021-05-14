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
Please refer to [this notebook](https://github.com/glasgowcompbio/pyMultiOmics/blob/main/notebooks/analysis_zebrafish_chebi.ipynb) for a demo.

Once mapping is completed, further analysis can be done.
For example, finding DE entities:
```
from pyMultiOmics.analysis import AnalysisPipeline

ap = AnalysisPipeline(m)

method = INFERENCE_DESEQ
ap.run_de(method, GENES, 'Distal', 'Proximal')
ap.run_de(method, GENES, 'Distal', 'Middle')
ap.run_de(method, GENES, 'Proximal', 'Middle')

method = INFERENCE_LIMMA
ap.run_de(method, PROTEINS, 'Distal', 'Proximal')
ap.run_de(method, PROTEINS, 'Distal', 'Middle')
ap.run_de(method, PROTEINS, 'Proximal', 'Middle')

method = INFERENCE_T_TEST
ap.run_de(method, COMPOUNDS, 'Distal', 'Proximal')
ap.run_de(method, COMPOUNDS, 'Distal', 'Middle')
ap.run_de(method, COMPOUNDS, 'Proximal', 'Middle')
```

Various queries can now be performed on the pipeline:

1. Retrieve a single node
```
from pyMultiOmics.query import QueryBuilder

node_id = '15366'
res = QueryBuilder(ap) \
        .add(Entity(node_id)) \
        .run()
res
```

2. Retrieve multiple nodes
```
node_id = ['15366', 'ENSDARG00000037781', 'F1QAA7']
res = QueryBuilder(ap) \
        .add(Entity(node_id)) \
        .run()
res
```

3. Retrieve nodes connected to a single node
```
query_id = 'F1QAA7'
res = QueryBuilder(ap) \
        .add(Entity(query_id)) \
        .add(Connected()) \
        .run()
res
```

4. Retrieve top-10 significantly changing genes
```
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
res
```

5. Find the compounds that are connected to the DE genes above

```
res = QueryBuilder(ap) \
        .add(Select(GENES)) \
        .add(SignificantDE(case, control, pval, fc_lte=fc_lte, fc_gte=fc_gte, N=N)) \
        .add(Connected(data_type=COMPOUNDS)) \
        .run()
res
```

6. Retrieve entity info
```
res = QueryBuilder(ap) \
        .add(Select(GENES)) \
        .add(SignificantDE(case, control, pval, fc_lte=fc_lte, fc_gte=fc_gte, N=N)) \
        .add(Connected()) \
        .add(Info()) \
        .run()
res
```