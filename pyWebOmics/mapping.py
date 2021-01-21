import itertools

import networkx as nx
from loguru import logger

from .common import as_list
from .constants import REACTIONS, PROTEOMICS, METABOLOMICS, GENOMICS, TRANSCRIPTOMICS, PATHWAYS, DataTypeDict
from .reactome import uniprot_to_reaction, compound_to_reaction


class Mapper():
    def __init__(self, species):
        self.species_list = [species]
        self.G = None

        self.has_transcript = False
        self.transcript_df = None
        self.transcript_design = None

        self.has_protein = False
        self.protein_df = None
        self.protein_design = None

        self.has_compound = False
        self.compound_df = None
        self.compound_design = None

    def set_transcript(self, transcript_df, transcript_design):
        self.transcript_df = transcript_df
        self.transcript_design = transcript_design
        self.has_transcript = True
        return self

    def set_protein(self, protein_df, protein_design):
        self.protein_df = protein_df
        self.protein_design = protein_design
        self.has_protein = True
        return self

    def set_compound(self, compound_df, compound_design):
        self.compound_df = compound_df
        self.compound_design = compound_design
        self.has_compound = True
        return self

    def build(self):
        self.G = nx.Graph()
        all_reactions = []

        protein_results = None
        if self.has_protein:
            protein_ids = list(self.protein_df.index.values)
            protein_results, _ = uniprot_to_reaction(protein_ids, self.species_list)
            all_reactions.extend(list(protein_results.values()))

        compound_results = None
        if self.has_compound:
            compound_ids = list(self.compound_df.index.values)
            compound_results, _ = compound_to_reaction(compound_ids, self.species_list)
            all_reactions.extend(list(compound_results.values()))

        # https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
        all_reactions = list(itertools.chain.from_iterable(all_reactions))

        # add nodes for reactions into G
        for reaction in all_reactions:
            reaction_id = reaction['reaction_id']
            reaction_name = reaction['reaction_name']
            self.G.add_node(reaction_id, name=reaction_name, type=REACTIONS)

        # add nodes and edges for proteins
        self._add_reaction_edges(protein_results, PROTEOMICS)

        # add nodes and edges for compounds
        self._add_reaction_edges(compound_results, METABOLOMICS)

        logger.info('Created a multi-omics network with %d nodes and %d edges' %
                    (len(self.G.nodes), len(self.G.edges)))
        logger.info('node_counts = %s' % self.num_nodes())
        return self

    def get_nodes(self, types=None):
        nodes = list(self.G.nodes(data=True))
        if types is None:
            return nodes
        types = as_list(types)

        # need to filter
        # below, node[0] is the id, while node[1] is the data dict
        results = filter(lambda node: node[1]['type'] in types, nodes)
        return list(results)

    def num_nodes(self, types=None):
        if types is None:
            types = [GENOMICS, TRANSCRIPTOMICS, PROTEOMICS, METABOLOMICS, REACTIONS, PATHWAYS]
        results = {DataTypeDict[t]: len(self.get_nodes(types=t)) for t in types}
        return results

    def get_neighbours(self, query_id, types=None):
        nbs = list(self.G.neighbors(query_id))
        if types is None:
            return nbs
        types = as_list(types)

        # need to filter by filter_types
        results = filter(lambda nb: self.G.nodes[nb]['type'] in types, nbs)
        return list(results)

    def get_connected(self, query_id, source_type, target_type):
        assert source_type in [PROTEOMICS, METABOLOMICS, REACTIONS]
        assert target_type in [PROTEOMICS, METABOLOMICS, REACTIONS]
        results = []

        # get anything that is directly connected to reactions (single-hop in the graph)
        if source_type == REACTIONS or target_type == REACTIONS:
            results = self.get_neighbours(query_id, types=target_type)

        # get connections between entites that need to go through reactions
        elif source_type in [PROTEOMICS, METABOLOMICS] and target_type in [PROTEOMICS, METABOLOMICS]:

            # get neighbouring reactions of the query_id node
            nb_reactions = self.get_neighbours(query_id, types=REACTIONS)

            # get neighbouring entities of all the reactions that are of target_type
            for reaction in nb_reactions:
                nb_entities = self.get_neighbours(reaction, types=target_type)
                results.extend(nb_entities)

        # return unique results
        return list(set(results))

    def _add_reaction_edges(self, query_results, node_type):
        for entity_id in query_results:
            self.G.add_node(entity_id, type=node_type)

            entity_reactions = query_results[entity_id]
            for reaction in entity_reactions:
                reaction_id = reaction['reaction_id']
                self.G.add_edge(entity_id, reaction_id)
