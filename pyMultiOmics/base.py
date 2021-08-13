from collections import OrderedDict

import pandas as pd
from loguru import logger

from .constants import DataTypeDict, MEASUREMENT_DF_LABEL, DESIGN_DF_LABEL


class SingleOmicsData():
    def __init__(self, data_type, measurement_df, design_df, feature_annot_df=None):
        self.data_type = data_type
        cleaned_measurement_df, cleaned_design_df = self._clean_data(measurement_df, design_df)
        msg = 'cleaned_measurement_df = %s, cleaned_design_df = %s' % (
            cleaned_measurement_df.shape,
            cleaned_design_df.shape
        )
        assert cleaned_measurement_df.shape[1] == cleaned_design_df.shape[0], msg

        self.initial_measurement_df = cleaned_measurement_df
        self.design_df = cleaned_design_df
        self.feature_annot_df = feature_annot_df
        self.processed_dfs = OrderedDict()

    @property
    def data_df(self):
        od = self.processed_dfs
        if len(od) > 0:
            return od[next(reversed(od))]
        else:
            return self.initial_measurement_df

    @property
    def data_label(self):
        try:
            dtype_str = DataTypeDict[self.data_type]
        except KeyError:
            dtype_str = self.data_type
        return dtype_str

    def get_initial_measurement_df(self):
        return self.initial_measurement_df.copy()

    def _clean_data(self, measurement_df, design_df):

        # drop duplicate rows and columns by values
        measurement_df = self._drop_dupes_by_values(measurement_df, MEASUREMENT_DF_LABEL)
        design_df = self._drop_dupes_by_values(design_df, DESIGN_DF_LABEL)

        # drop duplicate rows and columns by sample names
        measurement_df = self._drop_dupes_by_colnames(measurement_df, MEASUREMENT_DF_LABEL)
        design_df = self._drop_dupes_by_colnames(design_df.transpose(), DESIGN_DF_LABEL).transpose()

        # keep common samples having both measurements and metadata
        measurement_df, design_df = self._keep_common_samples(measurement_df, design_df)
        return measurement_df, design_df

    def _drop_dupes_by_values(self, df, label):
        no_dupe_rows = df.drop_duplicates(keep=False)
        no_dupe = no_dupe_rows.transpose().drop_duplicates(keep=False).transpose()
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
        dtype_str = self.data_label
        shape = self.data_df.shape
        return '%s data with (%d, %d) measurements' % (dtype_str, shape[0], shape[1])


class MultiOmicsData():
    def __init__(self, publication=None, url=None):
        self.views = {}
        self.publication = publication
        self.url = url

    def add_data(self, omics_data):
        if not isinstance(omics_data, list):
            items = [omics_data]
        else:
            items = omics_data

        for item in items:
            view_name = item.data_type
            self.views[view_name] = item

    def to_mofa(self):
        res = pd.DataFrame()
        for v in self.views:
            data = self.views[v]
            df = data.data_df
            df = df.melt(var_name='feature')
            df['view'] = data.data_label
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
