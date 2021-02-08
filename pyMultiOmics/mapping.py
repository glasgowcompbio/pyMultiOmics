import json

import networkx as nx
import pandas as pd
from loguru import logger

from .common import as_list
from .constants import REACTIONS, PROTEOMICS, METABOLOMICS, GENOMICS, TRANSCRIPTOMICS, PATHWAYS, DataTypeDict, \
    COMPOUND_DATABASE_CHEBI, MAPPING, PKS, IDS, NA, GENES_TO_PROTEINS, PROTEINS_TO_REACTIONS, COMPOUNDS_TO_REACTIONS, \
    REACTIONS_TO_PATHWAYS, GENES, PROTEINS, COMPOUNDS, TRANSCRIPTS
from .functions import reactome_mapping


class Mapper():
    def __init__(self, species, metabolic_pathway_only=True, compound_database_str=COMPOUND_DATABASE_CHEBI):
        self.species_list = [species]
        self.metabolic_pathway_only = metabolic_pathway_only
        self.compound_database_str = compound_database_str

        self.G = None

        self.gene_df = None
        self.gene_design = None

        self.protein_df = None
        self.protein_design = None

        self.compound_df = None
        self.compound_design = None

    def set_gene(self, gene_df, gene_design):
        self.gene_df = gene_df
        self.gene_design = gene_design
        return self

    def set_protein(self, protein_df, protein_design):
        self.protein_df = protein_df
        self.protein_design = protein_design
        return self

    def set_compound(self, compound_df, compound_design):
        self.compound_df = compound_df
        self.compound_design = compound_design
        return self

    def get_dfs(self, data_type):
        data_df = None
        design_df = None
        if data_type == GENOMICS:
            data_df = self.gene_df
            design_df = self.gene_design
        elif data_type == PROTEOMICS:
            data_df = self.protein_df
            design_df = self.protein_design
        elif data_type == METABOLOMICS:
            data_df = self.compound_df
            design_df = self.compound_design
        return data_df, design_df

    def build(self):

        # map different omics entities to reactome
        results = reactome_mapping(self.gene_df, self.protein_df, self.compound_df, self.compound_database_str,
                                   self.species_list, self.metabolic_pathway_only)
        designs = {
            GENOMICS: self.gene_design,
            PROTEOMICS: self.protein_design,
            METABOLOMICS: self.compound_design
        }

        # create network graphs, can be retrieved using self._get_graph()
        self.G = nx.Graph()

        # add nodes and edges from mapping results to the graph
        self._add_nodes(results, designs)
        self._add_edges(results)
        self._cleanup()

        logger.info('Created a multi-omics network with %d nodes and %d edges' %
                    (len(self.G.nodes), len(self.G.edges)))
        logger.info('node_counts = %s' % self.num_nodes())
        return self

    def get_node(self, node_id):
        graph = self._get_graph()
        return graph.nodes[node_id]

    def get_node_type(self, node_id):
        node = self.get_node(node_id)
        return node['data_type']

    def get_display_name(self, node_id):
        node = self.get_node(node_id)
        return node['display_name']

    def is_observed(self, node_id):
        node = self.get_node(node_id)
        return node['observed']

    def get_nodes(self, types=None):
        graph = self._get_graph()
        nodes = list(graph.nodes(data=True))
        if types is None:
            return nodes
        types = as_list(types)

        # need to filter
        # below, node[0] is the id, while node[1] is the data dict
        results = filter(lambda node: node[1]['data_type'] in types, nodes)
        return list(results)

    def num_nodes(self, types=None):
        if types is None:
            types = [GENOMICS, TRANSCRIPTOMICS, PROTEOMICS, METABOLOMICS, REACTIONS, PATHWAYS]
        else:
            types = as_list(types)
        results = {DataTypeDict[t]: len(self.get_nodes(types=t)) for t in types}
        return results

    def get_neighbours(self, query_id, include_types=None, exclude_types=None):
        graph = self._get_graph()
        results = list(graph.neighbors(query_id))

        # filter for types to be included
        if include_types is not None:
            include_types = as_list(include_types)
            results = filter(lambda node: self.get_node_type(node) in include_types, results)

        # filter for types to be excluded
        if exclude_types is not None:
            exclude_types = as_list(exclude_types)
            results = filter(lambda node: self.get_node_type(node) not in exclude_types, results)

        return list(results)

    def get_connected(self, query_id, dest_type=None, observed=None):
        dest_types = as_list(dest_type) if dest_type is not None else [GENES, TRANSCRIPTS,
                                                                       PROTEINS, COMPOUNDS,
                                                                       REACTIONS, PATHWAYS]
        possible_obs = as_list(observed) if observed is not None else [True, False]

        data = []
        for data_type in dest_types:
            connected = self._search_graph(query_id, data_type)
            for node_id in connected:
                if node_id == query_id:
                    continue
                node_info = self.get_node(node_id)
                display_name = node_info['display_name']
                obs = self.is_observed(node_id)
                row = [node_id, display_name, MAPPING[data_type], obs]

                if obs is None:  # reactions, pathways don't have any observed value
                    data.append(row)
                elif obs in possible_obs:  # everything else
                    data.append(row)

        df = pd.DataFrame(data, columns=['entity_id', 'display_name', 'data_type', 'observed'])
        df = df.set_index('entity_id')
        return df

    def _search_graph(self, query_id, dest_type):
        visited = []  # list of visited nodes
        queue = []  # nodes to process next
        seen_types = set()  # omics type seen so far

        visited.append(query_id)
        queue.append(query_id)

        results = []
        while queue:
            s = queue.pop(0)
            node_type = self.get_node_type(s)
            seen_types.add(node_type)

            if node_type == dest_type:
                # terminate search for s
                results.append(s)
                continue
            else:
                # keep searching in the neighbours of s
                # but exclude nodes of type that we've seen before
                neighbours = self.get_neighbours(s, exclude_types=seen_types)
                for node in neighbours:
                    if node not in visited:
                        visited.append(node)
                        queue.append(node)
        return results

    def _get_graph(self):
        return self.G

    def _json_str_to_df(self, json_str):
        res = json.loads(json_str)
        df = pd.DataFrame(res)
        return df

    def _filter_series(self, s, includes):
        if len(includes) == 0:
            return None

        new_dict = {}
        for index, value in s.items():
            if index in includes:
                new_dict[index] = value
        new_s = pd.Series(new_dict)
        return new_s

    def _insert_edge_if_node_exists(self, source_id, dest_id):
        graph = self._get_graph()
        if source_id in graph.nodes and dest_id in graph.nodes:
            graph.add_edge(source_id, dest_id)

    def _add_nodes(self, results, designs):
        graph = self._get_graph()
        data_types = [GENOMICS, PROTEOMICS, METABOLOMICS, REACTIONS, PATHWAYS]
        for data_type in data_types:

            display_data_type = MAPPING[data_type]
            logger.info('Processing nodes: %s' % display_data_type)

            # get identifier and display name columns
            id_col = PKS[data_type]
            name_col = IDS[data_type]

            # convert the dataframes that we need
            data_df = self._json_str_to_df(results[data_type])
            try:
                design_df = designs[data_type]
                sample_names = design_df.index.values.tolist()
            except KeyError:
                sample_names = []
            except AttributeError:
                sample_names = []

            # insert each row in dataframe into a network graph G
            for idx, row in data_df.iterrows():

                # get node identifier
                identifier = row[id_col]
                if identifier == NA:
                    continue

                # select only sample names columns to get the measurements
                observed = row['obs']
                measurements = None
                if observed:
                    measurements = self._filter_series(row, sample_names)

                # create node info dictionary
                node_info = {
                    'observed': observed,
                    'display_name': row[name_col],
                    'data_type': data_type,
                    'display_data_type': display_data_type,
                    'measurements': measurements,
                }

                # insert node into two separate graphs: G and G_complete
                # G_complete contains both actual and inferred nodes
                # G contains only the actual nodes
                graph.add_node(identifier, **node_info)

    def _add_edges(self, results):
        data_types = [GENES_TO_PROTEINS, PROTEINS_TO_REACTIONS, COMPOUNDS_TO_REACTIONS, REACTIONS_TO_PATHWAYS]
        for data_type in data_types:
            display_data_type = MAPPING[data_type]
            logger.info('Processing edges: %s' % display_data_type)

            df = self._json_str_to_df(results[data_type])
            for idx, row in df.iterrows():
                source_id, dest_id = row.values.tolist()

                # only add edges between nodes that already exist (have been added before)
                self._insert_edge_if_node_exists(source_id, dest_id)

    def _cleanup(self):

        # get all reactions nodes from the graph
        reactions = self.get_nodes(types=REACTIONS)
        to_delete = []
        for reaction_id, reaction_data in reactions:
            # FIXME: we shouldn't need to do this?
            # delete a reaction where the display_name is the same as the id
            # this would be the case when this reaction isn't found in the query to map reaction -> pathways
            # i.e. when only metabolic pathways are selected
            # in this case, this reaction actually shouldn't be included in the JSON results, but they are.
            # as a workaround we'd delete them here
            reaction_name = reaction_data['display_name']
            if reaction_id == reaction_name:
                to_delete.append(reaction_id)

        logger.debug('Deleting %d reactions that should not be there' % len(to_delete))
        logger.debug(to_delete[0:100])
        graph = self._get_graph()
        graph.remove_nodes_from(to_delete)
