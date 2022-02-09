import json

import networkx as nx
import pandas as pd
from loguru import logger

from .analysis import AnalysisPipeline
from .common import as_list
from .constants import REACTIONS, PROTEINS, COMPOUNDS, GENES, PATHWAYS, COMPOUND_DATABASE_CHEBI, PKS, IDS, NA, \
    GENES_TO_PROTEINS, PROTEINS_TO_REACTIONS, COMPOUNDS_TO_REACTIONS, REACTIONS_TO_PATHWAYS
from .functions import reactome_mapping
from .query import QueryBuilder, Connected, Entity


class Mapper():
    def __init__(self, multi_omics_data, species,
                 metabolic_pathway_only=True,
                 compound_database_str=COMPOUND_DATABASE_CHEBI,
                 include_related_chebi=False):
        self.multi_omics_data = multi_omics_data
        self.gene_df, self.gene_design = self.multi_omics_data.get_dfs(GENES)
        self.protein_df, self.protein_design = self.multi_omics_data.get_dfs(PROTEINS)
        self.compound_df, self.compound_design = self.multi_omics_data.get_dfs(COMPOUNDS)

        self.species_list = [species]
        self.metabolic_pathway_only = metabolic_pathway_only
        self.compound_database_str = compound_database_str
        self.include_related_chebi = include_related_chebi

        self.G = None

    def build(self):

        # map different omics entities to reactome
        results = reactome_mapping(self.gene_df, self.protein_df, self.compound_df, self.compound_database_str,
                                   self.species_list, self.metabolic_pathway_only, self.include_related_chebi)
        designs = {
            GENES: self.gene_design,
            PROTEINS: self.protein_design,
            COMPOUNDS: self.compound_design
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
            types = [GENES, PROTEINS, COMPOUNDS, REACTIONS, PATHWAYS]
        else:
            types = as_list(types)
        results = {t: len(self.get_nodes(types=t)) for t in types}
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
        """
        Retrieve mapped entities that are connected to query_id
        :param query_id: the node ID to query
        :param dest_type: connected entity types, e.g. GENES, PROTEINS, COMPOUNDS, REACTIONS, PATHWAYS
        :param observed: True to get entities observed in the data, False to retrieve everything mapped by Reactome
        :return: a dataframe of connected entities
        """
        ap = AnalysisPipeline(self.multi_omics_data, self)
        res = QueryBuilder(ap) \
            .add(Entity(query_id)) \
            .add(Connected(data_type=dest_type, observed=observed)) \
            .run()
        return res

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
        data_types = [GENES, PROTEINS, COMPOUNDS, REACTIONS, PATHWAYS]
        for data_type in data_types:

            display_data_type = data_type
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
            display_data_type = data_type
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
