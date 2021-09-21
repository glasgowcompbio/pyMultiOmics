from collections import OrderedDict

import pandas as pd
from loguru import logger

from .constants import MEASUREMENT_DF_LABEL, DESIGN_DF_LABEL, PADJ_COL_PREFIX, PVALUE_COL_PREFIX, FC_COL_PREFIX
from .info import get_info


class SingleOmicsData():
    def __init__(self, data_type, measurement_df, design_df, feature_annot_df=None, significant_df=None):
        self.data_type = data_type

        # extract columns containing p-values and fold changes into its own dataframe
        if significant_df is not None:
            self.significant_df = significant_df
        else:  # try to extract from the data
            keep_df, drop_df = self._get_significant_df(measurement_df)
            if not drop_df.empty:
                self.significant_df = drop_df
                measurement_df = keep_df

        # clean data
        cleaned_measurement_df, cleaned_design_df = self._clean_data(measurement_df, design_df)
        msg = 'cleaned_measurement_df = %s, cleaned_design_df = %s' % (
            cleaned_measurement_df.shape,
            cleaned_design_df.shape
        )
        assert cleaned_measurement_df.shape[1] == cleaned_design_df.shape[0], msg

        self.initial_measurement_df = cleaned_measurement_df
        self.design_df = cleaned_design_df
        self.feature_annot_df = feature_annot_df

        # An ordered dict of dataframes as it's processed through the pipeline.
        # The ordering is important since data is processed sequentially
        # The last entry is the current result to use
        self.processed_dfs = OrderedDict()

    @property
    def data_df(self):
        od = self.processed_dfs
        if len(od) > 0:
            return od[next(reversed(od))]
        else:
            return self.initial_measurement_df

    def get_initial_measurement_df(self):
        return self.initial_measurement_df.copy()

    def _get_significant_df(self, df):
        # old data in GraphOmics was using PADJ_COL_PREFIX rather than PVALUE_COL_PREFIX
        df = df.copy()
        df.columns = df.columns.str.replace(PADJ_COL_PREFIX, PVALUE_COL_PREFIX)

        # find significant information in columns starting with FC_COL_PREFIX, PVALUE_COL_PREFIX
        drop_cols = [FC_COL_PREFIX, PVALUE_COL_PREFIX]
        drop_cols = tuple([x.lower() for x in drop_cols])
        to_drop = list(filter(lambda x: x.lower().startswith(drop_cols), df.columns))
        to_keep = [col for col in df.columns if col not in to_drop]

        keep_df = df[to_keep]
        drop_df = df[to_drop]
        return keep_df, drop_df

    def _clean_data(self, measurement_df, design_df):

        # drop duplicate rows and columns by values
        # measurement_df = self._drop_dupes_by_values(measurement_df, MEASUREMENT_DF_LABEL)
        # design_df = self._drop_dupes_by_values(design_df, DESIGN_DF_LABEL) # don't do this!

        # drop duplicate rows and columns by sample names
        measurement_df = self._drop_dupes_by_colnames(measurement_df, MEASUREMENT_DF_LABEL)
        design_df = self._drop_dupes_by_colnames(design_df.transpose(), DESIGN_DF_LABEL).transpose()

        # keep common samples having both measurements and metadata
        measurement_df, design_df = self._keep_common_samples(measurement_df, design_df)
        return measurement_df, design_df

    def _drop_dupes_by_values(self, df, label):
        # drop duplicate rows, keep no duplicates
        no_dupe_rows = df.drop_duplicates(keep='first')

        # drop duplicate columns, keep no duplicates
        no_dupe = no_dupe_rows.transpose().drop_duplicates(keep='first').transpose()

        # print message is something has been dropped
        if df.shape != no_dupe.shape:
            logger.warning('Dropped duplicate from %s by values: %d rows and %d cols' % (
                label,
                df.shape[0] - no_dupe.shape[0],
                df.shape[1] - no_dupe.shape[1]
            ))
        return no_dupe

    def _drop_dupes_by_colnames(self, df, label):
        # find columns that have the same name
        # https://stackoverflow.com/questions/14984119/python-pandas-remove-duplicate-columns
        cleaned_df = df.loc[:, ~df.columns.duplicated()]

        n_cols_initial = df.shape[1]
        n_cols_cleaned = cleaned_df.shape[1]
        diff = n_cols_initial - n_cols_cleaned
        if diff > 0:
            logger.warning('Dropped %d duplicate sample names from %s' % (diff, label))
        return cleaned_df

    def _keep_common_samples(self, measurement_df, design_df):
        # find common sample names (rows and columns) in measurement and design dfs
        cols = measurement_df.columns.values
        rows = design_df.index.values
        common = set(cols).intersection(set(rows))

        # select the common row and col names
        selected_cols = [col for col in cols if col in common]
        selected_rows = [row for row in rows if row in common]
        cleaned_measurement_df = measurement_df[selected_cols]
        cleaned_design_df = design_df.loc[selected_rows]

        diff = measurement_df.shape[1] - cleaned_measurement_df.shape[1]
        if diff > 0:
            logger.warning('Dropped %d columns from measurement dataframe due to missing metadata' % diff)

        diff = design_df.shape[0] - cleaned_design_df.shape[0]
        if diff > 0:
            logger.warning('Dropped %d columns from sample metadata due to missing measurements' % diff)

        return cleaned_measurement_df, cleaned_design_df

    def __repr__(self):
        dtype_str = self.data_type
        shape = self.data_df.shape
        return '%s data with (%d, %d) measurements' % (dtype_str, shape[0], shape[1])


class MultiOmicsData():
    def __init__(self, publication=None, url=None):
        self.views = {}
        self.publication = publication
        self.url = url

    def add_data(self, omics_data):
        # check if list
        if not isinstance(omics_data, list):
            # if not, put it in a list
            items = [omics_data]
        else:
            items = omics_data

        # push list items to dictionary
        for item in items:
            view_name = item.data_type
            self.views[view_name] = item

    def has_data(self, data_type):
        return data_type in self.views

    def get_data(self, data_type):
        try:
            return self.views[data_type]
        except KeyError:
            return None

    def get_dfs(self, data_type):
        if self.has_data(data_type):
            data = self.get_data(data_type)
            return data.data_df, data.design_df
        else:
            return None, None

    def get_info(self, entity_id, data_type):
        return get_info(entity_id, data_type)

    def to_mofa(self):
        res = pd.DataFrame()
        for v in self.views:
            data = self.views[v]
            df = data.data_df
            df['feature'] = df.index
            df = df.melt(id_vars='feature',var_name='sample')
            df['view'] = data.data_type
            df = df.join(data.design_df, on = 'sample')
            res = res.append(df)
        return res

    def __repr__(self):
        msg = 'Multi-omics data container'
        if self.publication is not None:
            msg += '\n- publication: %s' % self.publication
        if self.url is not None:
            msg += '\n- URL: %s' % self.url
        if len(self.views) > 0:
            msg += '\n- Views: %d modalities' % len(self.views)
            for v in self.views:
                msg += '\n\t - %s' % self.views[v]
        return msg
