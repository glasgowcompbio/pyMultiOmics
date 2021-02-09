import pandas as pd

from .common import as_list
from .constants import REACTIONS, PATHWAYS, MAPPING, GENES, PROTEINS, COMPOUNDS, TRANSCRIPTS


class QueryBuilder():
    def __init__(self, mapping):
        self.mapping = mapping
        self.tasks = []

    def add(self, task):
        self.tasks.append(task)
        return self

    def run(self):
        res = None
        for task in self.tasks:
            res = task.run(self.mapping, res)
        return res


class Query():
    def run(self, mapping, previous_res):
        raise NotImplementedError


class SingleEntity(Query):
    def __init__(self, node_id):
        self.node_id = node_id

    def run(self, mapping, previous_res):
        node_id = self.node_id
        node_info = mapping.get_node(node_id)
        display_name = node_info['display_name']
        obs = mapping.is_observed(node_id)
        data_type = mapping.get_node_type(node_id)
        row = [node_id, display_name, MAPPING[data_type], obs]
        data = [row]
        df = pd.DataFrame(data, columns=['entity_id', 'display_name', 'data_type', 'observed'])
        df = df.set_index('entity_id')
        return df


class Connected(Query):
    def __init__(self, dest_type=None, observed=None):
        self.dest_type = dest_type
        self.observed = observed

    def run(self, mapping, previous_res):
        node_ids = previous_res.index.values
        dest_types = as_list(self.dest_type) if self.dest_type is not None else [GENES, TRANSCRIPTS,
                                                                                 PROTEINS, COMPOUNDS,
                                                                                 REACTIONS, PATHWAYS]
        possible_obs = as_list(self.observed) if self.observed is not None else [True, False]

        data = []
        for query_id in node_ids:
            for data_type in dest_types:
                connected = self._search_graph(mapping, query_id, data_type)
                for node_id in connected:
                    if node_id == query_id:
                        continue
                    node_info = mapping.get_node(node_id)
                    display_name = node_info['display_name']
                    obs = mapping.is_observed(node_id)
                    row = [node_id, display_name, MAPPING[data_type], obs]

                    if obs is None:  # reactions, pathways don't have any observed value
                        data.append(row)
                    elif obs in possible_obs:  # everything else
                        data.append(row)

        df = pd.DataFrame(data, columns=['entity_id', 'display_name', 'data_type', 'observed'])
        df = df.set_index('entity_id')
        return df

    def _search_graph(self, mapping, node_id, dest_type):
        visited = []  # list of visited nodes
        queue = []  # nodes to process next
        seen_types = set()  # omics type seen so far

        visited.append(node_id)
        queue.append(node_id)

        results = []
        while queue:
            s = queue.pop(0)
            node_type = mapping.get_node_type(s)
            seen_types.add(node_type)

            if node_type == dest_type:
                # terminate search for s
                results.append(s)
                continue
            else:
                # keep searching in the neighbours of s
                # but exclude nodes of type that we've seen before
                neighbours = mapping.get_neighbours(s, exclude_types=seen_types)
                for node in neighbours:
                    if node not in visited:
                        visited.append(node)
                        queue.append(node)
        return results
