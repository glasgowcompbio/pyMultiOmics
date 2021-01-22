import itertools

import networkx as nx
from loguru import logger

from .common import as_list
from .constants import REACTIONS, PROTEOMICS, METABOLOMICS, GENOMICS, TRANSCRIPTOMICS, PATHWAYS, DataTypeDict
from .reactome import uniprot_to_reaction, compound_to_reaction, ensembl_to_uniprot


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

        transcript_results = None
        if self.has_transcript:
            transcript_ids = list(self.transcript_df.index.values)
            self._add_nodes(transcript_ids, TRANSCRIPTOMICS, observed=True)
            transcript_results, _ = ensembl_to_uniprot(transcript_ids, self.species_list)

        protein_results = None
        if self.has_protein:
            protein_ids = list(self.protein_df.index.values)
            self._add_nodes(protein_ids, PROTEOMICS, observed=True)
            protein_results, _ = uniprot_to_reaction(protein_ids, self.species_list)
            all_reactions.extend(list(protein_results.values()))

        compound_results = None
        if self.has_compound:
            compound_ids = list(self.compound_df.index.values)
            self._add_nodes(compound_ids, METABOLOMICS, observed=True)
            compound_results, _ = compound_to_reaction(compound_ids, self.species_list)
            all_reactions.extend(list(compound_results.values()))

        # https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
        all_reactions = list(itertools.chain.from_iterable(all_reactions))

        # add nodes for reactions into G
        for reaction in all_reactions:
            reaction_id = reaction['reaction_id']
            reaction_name = reaction['reaction_name']
            self.G.add_node(reaction_id, name=reaction_name, type=REACTIONS)

        # add edges based on mapping results
        self._add_edges(transcript_results)
        self._add_edges(protein_results, dest_key='reaction_id')
        self._add_edges(compound_results, dest_key='reaction_id')

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
        for node in nodes:
            if 'type' not in node[1]:
                print(node)
        results = filter(lambda node: node[1]['type'] in types, nodes)
        return list(results)

    def num_nodes(self, types=None):
        if types is None:
            types = [GENOMICS, TRANSCRIPTOMICS, PROTEOMICS, METABOLOMICS, REACTIONS, PATHWAYS]
        results = {DataTypeDict[t]: len(self.get_nodes(types=t)) for t in types}
        return results

    def get_neighbours(self, query_id, types_include=None, types_exclude=None):
        results = list(self.G.neighbors(query_id))

        # filter for types to be included
        if types_include is not None:
            types_include = as_list(types_include)
            results = filter(lambda node: self.get_node_type(node) in types_include, results)

        # filter for types to be excluded
        if types_exclude is not None:
            types_exclude = as_list(types_exclude)
            results = filter(lambda node: self.get_node_type(node) not in types_exclude, results)

        return list(results)

    def get_node_type(self, node):
        return self.G.nodes[node]['type']

    def get_connected(self, query_id, dest_type):
        visited = []  # list of visited nodes
        queue = []  # nodes to process next
        seen_types = set()  # omics type seen so far

        visited.append(query_id)
        queue.append(query_id)

        results = []
        path = []
        while queue:
            s = queue.pop(0)
            node_type = self.get_node_type(s)
            seen_types.add(node_type)
            path.append((s, node_type))

            if node_type == dest_type:
                # terminate search for s
                results.append(s)
                continue
            else:
                # keep searching in the neighbours of s
                # but exclude nodes of type that we've seen before
                neighbours = self.get_neighbours(s, types_exclude=seen_types)
                for node in neighbours:
                    if node not in visited:
                        visited.append(node)
                        queue.append(node)

        return {
            'results': results,
            'path': path
        }

    def _add_nodes(self, entity_ids, node_type, observed=True):
        for entity_id in entity_ids:
            self.G.add_node(entity_id, type=node_type, obs=observed)

    def _add_edges(self, query_results, dest_key=None):
        for source_id in query_results:
            destinations = query_results[source_id]
            for dest in destinations:
                dest_id = dest[dest_key] if dest_key is not None else dest

                # only add edges between nodes that already exist (have been added before)
                if source_id in self.G.nodes and dest_id in self.G.nodes:
                    self.G.add_edge(source_id, dest_id)