import json

import networkx as nx
import pandas as pd
from loguru import logger

from .common import as_list
from .constants import REACTIONS, PROTEOMICS, METABOLOMICS, GENOMICS, TRANSCRIPTOMICS, PATHWAYS, DataTypeDict, \
    COMPOUND_DATABASE_CHEBI, MAPPING, PKS, IDS, NA, GENES_TO_PROTEINS, PROTEINS_TO_REACTIONS, COMPOUNDS_TO_REACTIONS, \
    REACTIONS_TO_PATHWAYS
from .functions import reactome_mapping


class Mapper():
    def __init__(self, species, metabolic_pathway_only=True, compound_database_str=COMPOUND_DATABASE_CHEBI):
        self.species_list = [species]
        self.metabolic_pathway_only = metabolic_pathway_only
        self.compound_database_str = compound_database_str

        self.G = None
        self.G_complete = None

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

    def build(self):

        # map different omics entities to reactome
        results = reactome_mapping(self.gene_df, self.protein_df, self.compound_df, self.compound_database_str,
                                   self.species_list, self.metabolic_pathway_only)
        designs = {
            GENOMICS: self.gene_design,
            PROTEOMICS: self.protein_design,
            METABOLOMICS: self.compound_design
        }

        # create network graphs
        # G_complete contains both types of nodes:
        # - actual: only the genes, proteins and compounds actually observed in the data
        # - inferred: other entities (genes, proteins and compounds) not observed in the data, but linked through
        #   their reactions to the actual nodes
        self.G = nx.Graph()
        self.G_complete = nx.Graph()

        # add nodes and edges from mapping results to the graph
        self._add_nodes(self.G, self.G_complete, results, designs)
        self._add_edges(self.G, self.G_complete, results)
        self._cleanup(self.G, self.G_complete)

        logger.info('Created a multi-omics network with %d nodes and %d edges' %
                    (len(self.G.nodes), len(self.G.edges)))
        logger.info('node_counts = %s' % self.num_nodes())
        return self

    def get_node(self, node_id, complete=False):
        graph = self._get_graph(complete)
        return graph.nodes[node_id]

    def get_node_type(self, node_id, complete=False):
        node = self.get_node(node_id, complete=complete)
        return node['data_type']

    def get_nodes(self, types=None, complete=False):
        graph = self._get_graph(complete)
        nodes = list(graph.nodes(data=True))
        if types is None:
            return nodes
        types = as_list(types)

        # need to filter
        # below, node[0] is the id, while node[1] is the data dict
        results = filter(lambda node: node[1]['data_type'] in types, nodes)
        return list(results)

    def num_nodes(self, types=None, complete=False):
        if types is None:
            types = [GENOMICS, TRANSCRIPTOMICS, PROTEOMICS, METABOLOMICS, REACTIONS, PATHWAYS]
        else:
            types = as_list(types)
        results = {DataTypeDict[t]: len(self.get_nodes(types=t, complete=complete)) for t in types}
        return results

    def get_neighbours(self, query_id, include_types=None, exclude_types=None, complete=False):
        graph = self._get_graph(complete)
        results = list(graph.neighbors(query_id))

        # filter for types to be included
        if include_types is not None:
            include_types = as_list(include_types)
            results = filter(lambda node: self.get_node_type(node, complete=complete) in include_types, results)

        # filter for types to be excluded
        if exclude_types is not None:
            exclude_types = as_list(exclude_types)
            results = filter(lambda node: self.get_node_type(node, complete=complete) not in exclude_types, results)

        return list(results)

    def get_connected(self, query_id, dest_type, complete=False):
        visited = []  # list of visited nodes
        queue = []  # nodes to process next
        seen_types = set()  # omics type seen so far

        visited.append(query_id)
        queue.append(query_id)

        results = []
        while queue:
            s = queue.pop(0)
            node_type = self.get_node_type(s, complete=complete)
            seen_types.add(node_type)

            if node_type == dest_type:
                # terminate search for s
                results.append(s)
                continue
            else:
                # keep searching in the neighbours of s
                # but exclude nodes of type that we've seen before
                neighbours = self.get_neighbours(s, exclude_types=seen_types, complete=complete)
                for node in neighbours:
                    if node not in visited:
                        visited.append(node)
                        queue.append(node)
        return results

    def _get_graph(self, complete):
        graph = self.G_complete if complete else self.G
        return graph

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

    def _insert_edge_if_node_exists(self, graph, source_id, dest_id):
        if source_id in graph.nodes and dest_id in graph.nodes:
            graph.add_edge(source_id, dest_id)

    def _add_nodes(self, G, G_complete, results, designs):
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
                measurements = self._filter_series(row, sample_names)

                # create node info dictionary
                observed = row['obs']
                node_info = {
                    'obs': observed,
                    'display_name': row[name_col],
                    'data_type': data_type,
                    'display_data_type': display_data_type,
                    'measurements': measurements,
                }

                # insert node into two separate graphs: G and G_complete
                # G_complete contains both actual and inferred nodes
                # G contains only the actual nodes
                G_complete.add_node(identifier, **node_info)
                if data_type in [REACTIONS, PATHWAYS]:
                    # for reactions, pathways, observed is always None
                    # so we just insert them anyway
                    G.add_node(identifier, **node_info)
                else:
                    # for other data type, e.g. genes, proteins, compounds
                    # make sure it's observed before adding them to G
                    if observed:
                        G.add_node(identifier, **node_info)

    def _add_edges(self, G, G_complete, results):
        data_types = [GENES_TO_PROTEINS, PROTEINS_TO_REACTIONS, COMPOUNDS_TO_REACTIONS, REACTIONS_TO_PATHWAYS]
        for data_type in data_types:
            display_data_type = MAPPING[data_type]
            logger.info('Processing edges: %s' % display_data_type)

            df = self._json_str_to_df(results[data_type])
            for idx, row in df.iterrows():
                source_id, dest_id = row.values.tolist()

                # only add edges between nodes that already exist (have been added before)
                self._insert_edge_if_node_exists(G, source_id, dest_id)
                self._insert_edge_if_node_exists(G_complete, source_id, dest_id)

    def _cleanup(self, G, G_complete):

        # get all reactions nodes from the incomplete graph
        complete = False
        reactions = self.get_nodes(types=REACTIONS, complete=complete)

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
        graph = self._get_graph(True)
        graph.remove_nodes_from(to_delete)
        graph = self._get_graph(False)
        graph.remove_nodes_from(to_delete)

        self._delete_disconnected_nodes(True)
        self._delete_disconnected_nodes(False)

    def _delete_disconnected_nodes(self, complete):
        logger.debug('complete = %s' % complete)
        for type in [PATHWAYS, REACTIONS, PROTEOMICS, METABOLOMICS, TRANSCRIPTOMICS, GENOMICS]:
            nodes = self.get_nodes(types=type, complete=complete)
            to_delete = []
            for node_id, node_data in nodes:
                if len(self.get_neighbours(node_id, complete=complete)) == 0:
                    to_delete.append(node_id)

            logger.debug('Deleted %d disconnected nodes of type %d' % (len(to_delete), type))
            graph = self._get_graph(complete)
            graph.remove_nodes_from(to_delete)
