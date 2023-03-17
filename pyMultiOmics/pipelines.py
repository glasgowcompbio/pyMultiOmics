import numpy as np
import pandas as pd
import pylab as plt
import seaborn as sns
from loguru import logger
from scipy import stats
from sklearn import preprocessing
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
from statsmodels.sandbox.stats.multicomp import multipletests

from pyMultiOmics.constants import GROUP_COL, IDS, PKS, SMALL


class Inference(object):
    def __init__(self, data_df, design_df, data_type, remove_cols=None, min_value=SMALL,
                 replace_mean=True):

        data_df = data_df.copy()

        # remove all the default columns from dataframe if nothing provided
        try:
            if remove_cols is None:
                remove_cols = ['padj_', 'FC_', 'significant_', 'obs', PKS[data_type],
                               IDS[data_type]]
            remove_cols = tuple([x.lower() for x in remove_cols])
            to_drop = list(filter(lambda x: x.lower().startswith(remove_cols), data_df.columns))
            df = data_df.drop(to_drop, axis=1)
        except KeyError:
            df = data_df

        # remove rows that are all NAs and all 0s
        df = df.dropna(how='all')
        df = df.loc[~(df == 0).all(axis=1)]

        self.data_df = df
        self.design_df = design_df
        self.data_type = data_type

        # data imputation:
        # - if all zero in group then replace with min_value
        # - replace any other zeros with mean of group
        self._impute_data(min_value, replace_mean=replace_mean)

    def _impute_data(self, min_value, replace_mean=True):
        if self.design_df is not None:
            grouping = self.design_df.groupby('group')
            for group, samples in grouping.groups.items():

                # If all zero in group then replace with minimum
                temp = self.data_df.loc[:, samples]
                temp = (temp == 0).all(axis=1)
                self.data_df.loc[temp, samples] = min_value

                if replace_mean:
                    # replace any other zeros with mean of group
                    subset_df = self.data_df.loc[:, samples]
                    self.data_df.loc[:, samples] = subset_df.mask(subset_df == 0,
                                                                  subset_df.mean(axis=1), axis=0)

        # replace all 0s with min_value
        self.data_df = self.data_df.replace(0, min_value)

    def normalise(self, kind, log, std_method):
        # df is features X samples, and it's normalised along the feature axis
        # if doing PCA on the samples, transpose the df so it goes the right way
        df = self._normalise_df(self.data_df, log=log, method=std_method)
        df = df.transpose() if kind == 'samples' else df
        return df

    def _normalise_df(self, data_df, log=False, method='standard'):
        if data_df.empty or method is None:
            return data_df

        data_arr = data_df.values
        if log:
            data_arr = np.log(data_arr)

        assert method in ['standard', 'minmax']
        if method == 'standard':
            # center data to have 0 mean and unit variance for heatmap and pca
            # data_arr is features X samples, so axis=1 is to normalise along the features
            scaled_data = preprocessing.scale(data_arr, axis=1)

        elif method == 'minmax':
            # normalise features to be within 0 to 1
            scaler = MinMaxScaler()
            scaled_data = scaler.fit_transform(data_arr.transpose())
            scaled_data = scaled_data.transpose()

        # set the values back to the dataframe
        sample_names = data_df.columns
        data_df[sample_names] = scaled_data
        return data_df

    def run_deseq_features(self, keep_threshold, case, control):
        logger.info('DeSEQ2 support is not available, please install rpy2 extra package')

    def run_limma_features(self, case, control):
        logger.info('limma support is not available, please install rpy2 extra package')

    def run_ttest_features(self, case, control, log=False):
        logger.info('t-test case is %s, control is %s' % (case, control))
        count_data = self.data_df
        col_data = self.design_df
        sample_group = col_data[col_data[GROUP_COL] == case]
        case_data = count_data[sample_group.index]
        sample_group = col_data[col_data[GROUP_COL] == control]
        control_data = count_data[sample_group.index]

        nrow, _ = count_data.shape
        pvalues = []
        lfcs = []
        indices = []

        # T-test for the means of two independent samples
        for i in range(nrow):

            case = case_data.iloc[i, :].values
            control = control_data.iloc[i, :].values
            idx = count_data.index[i]

            # remove 0 values, which were originally NA when exported from PiMP
            case = case[case != 0]
            control = control[control != 0]

            # log the data if it isn't already logged
            if log:
                case_log = np.log2(case)
                control_log = np.log2(control)
            else:
                case_log = case
                control_log = control

            statistics, pvalue = stats.ttest_ind(case_log, control_log)
            if not np.isnan(pvalue):
                lfc = np.mean(case_log) - np.mean(control_log)
                pvalues.append(pvalue)
                lfcs.append(lfc)
                indices.append(idx)

        # correct p-values
        reject, pvals_corrected, _, _ = multipletests(pvalues, method='fdr_bh')
        result_df = pd.DataFrame({
            'padj': pvals_corrected,
            'log2FoldChange': lfcs
        }, index=indices)
        return result_df

    def PCA(self, normalise=None, log=False, n_components=10, hue=None, style=None,
            palette='bright', return_fig=False, kind='samples'):
        assert kind in ['samples', 'features']

        df = self.normalise(kind, log, normalise)

        pca = PCA(n_components=n_components)
        pcs = pca.fit_transform(df)
        pc1_values = pcs[:, 0]
        pc2_values = pcs[:, 1]

        sns.set_context('poster')
        fig = plt.figure(figsize=(10, 5))

        palette = None if hue is None else palette
        g = sns.scatterplot(x=pc1_values, y=pc2_values, hue=hue, style=style, palette=palette)
        if hue is not None:
            g.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fontsize=12)

        print('PCA explained variance', pca.explained_variance_ratio_.cumsum())

        if return_fig:
            return pc1_values, pc2_values, fig
        else:
            return pc1_values, pc2_values

    def heatmap(self, N=None, normalise=None, log=False, kind='samples'):
        assert kind in ['samples', 'features']

        df = self.normalise(kind, log, normalise)

        # select the first N rows for heatmap
        if N is not None:
            df = df[0:N]

        sns.set_context('poster')
        plt.figure(figsize=(10, 10))
        sns.heatmap(df)

    def cluster(self, n_clusters, normalise=None, log=False, kind='samples'):
        assert kind in ['samples', 'features']

        df = self.normalise(kind, log, normalise)

        # initialize the KMeans model with the number of clusters
        kmeans = KMeans(n_clusters=n_clusters)

        # fit the data to the model
        kmeans.fit(df)
        labels = kmeans.labels_
        centroids = kmeans.cluster_centers_

        return labels, centroids
