import numpy as np
import pandas as pd
from loguru import logger

from .common import as_list
from .constants import REACTIONS, PATHWAYS, GENES, PROTEINS, COMPOUNDS, QUERY_DISPLAY_NAME, QUERY_NODE_ID, \
    QUERY_DATA_TYPE, QUERY_OBSERVED, QUERY_ENTITY_ID, QUERY_SOURCE_ID
from .info import get_info


class QueryBuilder():
    def __init__(self, pipeline):
        self.pipeline = pipeline
        self.tasks = []
        self.results = []

    def add(self, task):
        self.tasks.append(task)
        return self

    def run(self):
        self.results = []
        prev_task = None
        for task in self.tasks:
            # run task and get its result
            task.run(self.pipeline, prev_task)
            res = task.get_result()

            # track all results generated so far
            self.results.append(res)
            prev_task = task

        # return the last result, if any
        if len(self.results) > 0:
            return self.results[-1]
        else:
            return None


class Query():
    def __init__(self):
        self.result = None

    def get_result(self):
        return self.result

    def run(self, pipeline, previous_query):
        raise NotImplementedError

    def _get_node_data(self, mapping, node_id):
        node_info = mapping.get_node(node_id)
        display_name = node_info[QUERY_DISPLAY_NAME]
        obs = mapping.is_observed(node_id)
        data_type = mapping.get_node_type(node_id)
        return {
            QUERY_NODE_ID: node_id,
            QUERY_DISPLAY_NAME: display_name,
            QUERY_DATA_TYPE: data_type,
            QUERY_OBSERVED: obs
        }

    def _to_row(self, node_data):
        node_id = node_data[QUERY_NODE_ID]
        display_name = node_data[QUERY_DISPLAY_NAME]
        data_type = node_data[QUERY_DATA_TYPE]
        obs = node_data[QUERY_OBSERVED]
        row = [node_id, display_name, data_type, obs]
        return row

    def _to_df(self, data):
        columns = [QUERY_ENTITY_ID, QUERY_DISPLAY_NAME, QUERY_DATA_TYPE, QUERY_OBSERVED]
        df = pd.DataFrame(data, columns=columns)
        df = df.set_index(QUERY_ENTITY_ID)
        return df


class Entity(Query):
    def __init__(self, node_ids):
        super().__init__()
        self.node_ids = as_list(node_ids)

    def run(self, pipeline, _):
        mapping = pipeline.mapping
        data = []
        for node_id in self.node_ids:
            try:
                node_data = self._get_node_data(mapping, node_id)
                row = self._to_row(node_data)
                data.append(row)
            except KeyError:
                logger.warning(" %s wasn't mapped, it may not be present in Reactome" % node_id)
        df = self._to_df(data)
        self.result = df


class Select(Query):
    def __init__(self, data_types):
        super().__init__()
        self.data_types = as_list(data_types)

    def run(self, pipeline, previous_query):
        mapping = pipeline.mapping
        # get all nodes of data_type from the mapping object
        res = mapping.get_nodes(types=self.data_types)
        node_ids = [node_id for node_id, node_data in res]

        # select entities having node_ids, and save the result
        e = Entity(node_ids)
        e.run(pipeline, previous_query)
        self.result = e.get_result()


class Connected(Query):
    def __init__(self, data_type=None, observed=None):
        super().__init__()
        self.data_types = as_list(data_type) if data_type is not None else [GENES, PROTEINS, COMPOUNDS,
                                                                            REACTIONS, PATHWAYS]
        self.observed = observed

    def run(self, pipeline, previous_query):
        mapping = pipeline.mapping
        previous_res = previous_query.get_result()
        node_ids = previous_res.index.values
        possible_obs = as_list(self.observed) if self.observed is not None else [True, False]

        data = []
        for node_id in node_ids:
            for data_type in self.data_types:
                connected = self._search_graph(mapping, node_id, data_type)
                for connected_id in connected:
                    if connected_id == node_id:
                        continue

                    node_data = self._get_node_data(mapping, connected_id)
                    row = self._to_row(node_data) + [node_id]

                    obs = node_data['observed']
                    if obs is None:  # reactions, pathways don't have any observed value
                        data.append(row)
                    elif obs in possible_obs:  # everything else
                        data.append(row)

        df = self._to_df(data).drop_duplicates()
        self.result = df

    def _to_df(self, data):
        columns = [QUERY_ENTITY_ID, QUERY_DISPLAY_NAME, QUERY_DATA_TYPE, QUERY_OBSERVED, QUERY_SOURCE_ID]
        df = pd.DataFrame(data, columns=columns)
        df = df.set_index(QUERY_ENTITY_ID)
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


class SignificantDE(Query):
    def __init__(self, case, control, pval, fc_lte=None, fc_gte=None, N=None, ascending=None):
        self.case = case
        self.control = control
        self.pval = pval
        self.fc_lte = fc_lte
        self.fc_gte = fc_gte
        self.N = N
        self.ascending = ascending

    def run(self, pipeline, previous_query):
        # get data type from the previous query, ensure it's only for one type
        data_types = previous_query.data_types
        if len(data_types) > 1:
            raise ValueError('SignificantDE can only be used with one data type')
        data_type = data_types[0]

        # get the DE result for the case vs control comparison
        de_res = pipeline.get_de_results(data_type, self.case, self.control)

        # filter results
        pval_col = 'padj_%s_vs_%s' % (self.case, self.control)
        fc_col = 'FC_%s_vs_%s' % (self.case, self.control)
        if self.fc_lte is None and self.fc_gte is None:
            # filter by p-value only
            filtered = de_res[de_res[pval_col] <= self.pval]
        else:
            # filter by p-value and fold-change
            fc_lte = -np.inf if self.fc_lte is None else self.fc_lte
            fc_gte = np.inf if self.fc_gte is None else self.fc_gte
            filtered = de_res[
                (de_res[pval_col] <= self.pval) & ((de_res[fc_col] <= fc_lte) | (de_res[fc_col] >= fc_gte))]

        # merge with previous selection result
        previous_res = previous_query.get_result()
        merged = pd.merge(previous_res, filtered, how='inner', left_index=True, right_index=True)

        # if the sort order is not specified, take N/2 from the top and bottom halves of FC result, and combine them
        if self.ascending is None:
            N = int(self.N / 2)
            top = merged.sort_values(fc_col, ascending=False).head(N)
            bottom = merged.sort_values(fc_col, ascending=True).head(N)
            combined = pd.concat([top, bottom], axis=0)
            final_res = combined.sort_values(fc_col, ascending=False)
        else:
            # else take the top N result according to the sort order
            final_res = merged.sort_values(fc_col, ascending=self.ascending).head(self.N)
        self.result = final_res


class Info(Query):
    def run(self, pipeline, previous_query):
        previous_res = previous_query.get_result()
        copy_df = previous_res.copy()
        node_info = []
        node_images = []
        node_links = []
        for node_id, node_data in copy_df.iterrows():
            data_type = node_data[QUERY_DATA_TYPE]
            res = get_info(node_id, data_type)
            node_info.append(res['infos'])
            node_images.append(res['images'])
            node_links.append(res['links'])

        copy_df['infos'] = node_info
        copy_df['images'] = node_images
        copy_df['links'] = node_links
        self.result = copy_df
