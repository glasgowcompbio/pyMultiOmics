import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pylab as plt
import seaborn as sns
from adjustText import adjust_text
from kneed import KneeLocator
from loguru import logger
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_samples, silhouette_score

from .constants import PROTEINS, COMPOUNDS, REACTIONS, PATHWAYS, SMALL, GENES, \
    INFERENCE_T_TEST, INFERENCE_DESEQ, INFERENCE_LIMMA, CIC_COMPOUNDS
from .pipelines import Inference


class AnalysisPipeline(object):
    def __init__(self, multi_omics_data, mapping):
        self.multi_omics_data = multi_omics_data
        self.mapping = mapping

        # stores the results from DE analyses
        self.de_analyses = {
            GENES: [],
            PROTEINS: [],
            COMPOUNDS: [],
            REACTIONS: [],
            PATHWAYS: [],
            CIC_COMPOUNDS: []
        }

    def get_mapping(self):
        return self.mapping

    def run_de(self, selected_method, data_type, case, control):
        method_to_run = {
            INFERENCE_T_TEST: TTestAnalysis,
            INFERENCE_DESEQ: DESeq2Analysis,
            INFERENCE_LIMMA: LimmaAnalysis
        }

        # choose one analysis, depending on the selected method
        chosen_analysis = method_to_run[selected_method]

        # check if already exists
        found = False
        for analysis in self.de_analyses[data_type]:
            if analysis.case == case and analysis.control == control:
                found = True

        if not found:
            data_df, design_df = self.multi_omics_data.get_dfs(data_type)
            if data_df is not None and design_df is not None:
                # actually instantiate the DE Analysis here
                analysis = chosen_analysis(data_df, design_df, data_type, case, control)
                analysis.run()
                self.de_analyses[data_type].append(analysis)
        else:
            logger.warning('DE analysis for data_type=%s case=%s control=%s already exists' % (
                data_type, case, control
            ))

    def get_de_analysis(self, data_type, case, control, method=None):
        if method is not None:
            assert method in [INFERENCE_T_TEST, INFERENCE_DESEQ, INFERENCE_LIMMA]

        found = None
        for analysis in self.de_analyses[data_type]:
            # select only analysis for a particular case and control, if specified
            if analysis.case != case or analysis.control != control:
                continue

            # select method, if specified
            if method is not None and analysis.get_analysis_type() != method:
                continue

            # return the first result we found
            found = analysis
        return found

    def get_de_results(self, data_type, case, control, method=None):
        # get results dataframe for the DE analysis
        analysis = self.get_de_analysis(data_type, case, control, method=method)
        res = analysis.get_results()

        # rename the column 'log2FoldChange' to 'FC'
        res = res.rename(columns={
            'log2FoldChange': 'FC'
        })

        # add suffix to 'padj' and 'FC' columns
        suffix = '_%s_vs_%s' % (analysis.case, analysis.control)
        res = res.add_suffix(suffix)
        return res.sort_index()

    def de_sort_and_filter(self, df, p_value_colname, p_value_thresh, fc_colname,
                           fc_sort_order='asc', top_n=None, fc_iqr_thresh=1.5):

        # sort the DataFrame by "FC" column in ascending or descending order
        assert fc_sort_order in ['asc', 'desc']
        ascending = True if fc_sort_order == 'asc' else False
        df_sorted = df.sort_values(fc_colname, ascending=ascending)

        # filter the DataFrame to only include rows with a p-value less than the given threshold
        # and with log-fold change outliers removed
        sig_mask, fc_mask, combined_mask = self.filter_de_df(
            df_sorted, fc_colname, p_value_colname, p_value_thresh,
            fc_iqr_thresh=fc_iqr_thresh)
        df_filtered = df_sorted[combined_mask]

        # display only the top-N rows, if requested
        if top_n is not None:
            df_filtered = df_filtered.head(top_n)

        return df_filtered

    def PCA(self, data_type, normalise=None, log=False, n_components=10, hue=None, style=None,
            palette='bright', return_fig=False, kind='samples'):
        assert kind in ['samples', 'features']

        data_df, design_df = self.multi_omics_data.get_dfs(data_type)
        min_replace = SMALL
        wi = Inference(data_df, design_df, data_type, min_value=min_replace)
        df = wi.normalise(kind, log, normalise)

        pca = PCA(n_components=n_components)
        pcs = pca.fit_transform(df)
        pc1_values = pcs[:, 0]
        pc2_values = pcs[:, 1]
        # print('PCA explained variance', pca.explained_variance_ratio_.cumsum())

        if return_fig:
            sns.set_context('poster')
            fig = plt.figure(figsize=(10, 5))
            palette = None if hue is None else palette
            g = sns.scatterplot(x=pc1_values, y=pc2_values, hue=hue, style=style, palette=palette)
            if hue is not None:
                g.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fontsize=12)

        return pcs

    def scatter_PCA(self, data_type, normalise=None, log=False, n_components=10, hue=None,
                    style=None,
                    palette='bright', return_fig=False, kind='samples'):
        assert kind in ['samples', 'features']

        data_df, design_df = self.multi_omics_data.get_dfs(data_type)
        min_replace = SMALL
        wi = Inference(data_df, design_df, data_type, min_value=min_replace)
        df = wi.normalise(kind, log, normalise)

        pca = PCA(n_components=n_components)
        pcs = pca.fit_transform(df)
        print('PCA explained variance', pca.explained_variance_ratio_.cumsum())

        if return_fig:
            sns.set_context('poster')
            palette = None if hue is None else palette

            # Set the number of rows and columns for the subplots
            n_rows = n_components - 1
            n_cols = n_components - 1
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))

            for i in range(n_rows):
                for j in range(i + 1, n_cols + 1):
                    ax = axes[i, j - 1]
                    g = sns.scatterplot(x=pcs[:, i], y=pcs[:, j], hue=hue, style=style,
                                        palette=palette, ax=ax, s=20)
                    ax.set_xlabel(f'PC{i + 1}')
                    ax.set_ylabel(f'PC{j + 1}')

                    # Move the legend outside the plot for the last subplot in the row
                    if j == n_cols:
                        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fontsize=12)
                    else:
                        ax.legend().set_visible(False)

            plt.tight_layout()
            plt.show()

        return pcs

    def heatmap(self, data_type, N=None, normalise=None, log=False, return_fig=False,
                kind='samples', show_ticks=True, selected_cluster=None, cluster_labels=None,
                cluster_order=None):
        assert kind in ['samples', 'features']

        data_df, design_df = self.multi_omics_data.get_dfs(data_type)
        min_replace = SMALL
        wi = Inference(data_df, design_df, data_type, min_value=min_replace)
        df = wi.normalise(kind, log, normalise)

        # select the first N rows for heatmap
        if N is not None:
            df = df[0:N]

        # reorder rows according to the provided cluster ordering
        if cluster_order is not None:
            df = df.iloc[cluster_order, :]

        # select a cluster to display in the heatmap
        if selected_cluster is not None and cluster_labels is not None:
            # select indices of cluster_labels that are the same as selected_cluster
            selected_indices = np.where(cluster_labels == selected_cluster)[0]

            # then select rows in df that correspond to the indices above
            df = df.iloc[selected_indices]

        if return_fig:
            sns.set_context('poster')
            fig = plt.figure(figsize=(10, 10))
            sns.heatmap(df, cmap='seismic')
            if kind == 'samples':
                plt.xlabel('Analytes')
                plt.ylabel('Samples')
            else:
                plt.xlabel('Samples')
                plt.ylabel('Analytes')
            if not show_ticks:
                plt.xticks([])
                plt.yticks([])

        return df

    def cluster(self, data_type, normalise=None, log=False, kind='samples', return_fig=False,
                k=None, pcs=None):
        assert kind in ['samples', 'features']

        if pcs is None:
            data_df, design_df = self.multi_omics_data.get_dfs(data_type)
            min_replace = SMALL
            wi = Inference(data_df, design_df, data_type, min_value=min_replace)
            df = wi.normalise(kind, log, normalise)
            data = df.values
        else:
            data = pcs

        # use our own k
        silhouette_scores = []
        if k is not None:
            n_clusters = k
        else:  # use elbow method to choose n_clusters
            n_clusters, silhouette_scores = self.choose_n_cluster(data, pcs, return_fig)

        # final clustering using best_n_clusters
        kmeans = KMeans(n_clusters=n_clusters, random_state=0)
        kmeans.fit(data)
        labels = kmeans.labels_
        centroids = kmeans.cluster_centers_
        order = np.argsort(labels)
        silhouette_scores = np.array(silhouette_scores)

        return labels, order, centroids, silhouette_scores

    def choose_n_cluster(self, data, pcs, return_fig):

        possible_clusters = range(2, 11)
        silhouette_scores = self.compute_silhouette_scores(data, possible_clusters)
        n_clusters = self.detect_elbow(possible_clusters, silhouette_scores)

        # plot the scores and silhouette analysis
        if return_fig:
            # Plot silhouette scores line plot
            self.silhouette_lineplot(n_clusters, possible_clusters, silhouette_scores)

            # Plot silhouette analysis for each n_clusters_to_test
            self.silhouette_analysis(data, n_clusters, pcs, silhouette_scores)
        return n_clusters, silhouette_scores

    def compute_silhouette_scores(self, data, possible_clusters):
        silhouette_scores = []
        for n_clusters_to_test in possible_clusters:
            kmeans = KMeans(n_clusters=n_clusters_to_test, random_state=0)
            labels = kmeans.fit_predict(data)
            silhouette_avg = silhouette_score(data, labels)
            silhouette_scores.append(silhouette_avg)
        # print(silhouette_scores)
        return silhouette_scores

    def detect_elbow(self, possible_clusters, silhouette_scores):

        # Create KneeLocator instances for all combinations of curve shapes and directions
        combinations = [
            ('convex', 'decreasing'),
            ('concave', 'increasing'),
            ('convex', 'increasing'),
            ('concave', 'decreasing')
        ]
        elbow_candidates = []
        for curve, direction in combinations:
            knee_locator = KneeLocator(x=possible_clusters, y=silhouette_scores,
                                       curve=curve, direction=direction)
            elbow = knee_locator.knee
            if elbow is not None:
                elbow_candidates.append(
                    (elbow, silhouette_scores[possible_clusters.index(elbow)]))

        # Choose the elbow point with the highest silhouette score
        if elbow_candidates:
            # print(elbow_candidates)
            sorted_elbow_candidates = sorted(elbow_candidates, key=lambda x: x[1], reverse=True)
            top_two = sorted_elbow_candidates[:2]

            # Prioritize ('convex', 'decreasing') if it's in the top 2
            convex_decreasing_elbow = [elbow for elbow, _ in elbow_candidates if
                                       KneeLocator(x=possible_clusters, y=silhouette_scores,
                                                   curve='convex',
                                                   direction='decreasing').knee == elbow]
            if convex_decreasing_elbow and convex_decreasing_elbow[0] in [elbow for elbow, _ in
                                                                          top_two]:
                n_clusters = convex_decreasing_elbow[0]
            else:
                n_clusters, _ = top_two[0]
        else:
            # Fallback/default value for n_clusters in case no elbow point is detected
            n_clusters = 3
        return n_clusters

    def silhouette_lineplot(self, n_clusters, possible_clusters, silhouette_scores):
        silhouette_df = pd.DataFrame({'Score': silhouette_scores})
        fig = plt.figure(figsize=(10, 5))
        ax = sns.lineplot(data=silhouette_df)
        plt.title('Best number of clusters = %d' % n_clusters)
        ax.set_xticks(np.arange(len(possible_clusters)))
        ax.set_xticklabels(possible_clusters)
        plt.show()

    def silhouette_analysis(self, data, n_clusters, pcs, silhouette_scores):
        if n_clusters <= 2:
            to_plot = [n_clusters, n_clusters + 1, n_clusters + 2]
        else:
            to_plot = [n_clusters - 1, n_clusters, n_clusters + 1]
        for n_clusters_to_test in to_plot:
            # Create a subplot with 1 row and 2 columns
            fig, (ax1, ax2) = plt.subplots(1, 2)
            fig.set_size_inches(18, 7)

            # The 1st subplot is the silhouette plot
            ax1.set_xlim([-0.1, 1])
            ax1.set_ylim([0, len(data) + (n_clusters_to_test + 1) * 10])

            clusterer = KMeans(n_clusters=n_clusters_to_test, random_state=0)
            cluster_labels = clusterer.fit_predict(data)

            silhouette_avg = silhouette_scores[n_clusters_to_test - 2]

            # Compute the silhouette scores for each sample
            sample_silhouette_values = silhouette_samples(data, cluster_labels)

            y_lower = 10
            for i in range(n_clusters_to_test):
                ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]
                ith_cluster_silhouette_values.sort()

                size_cluster_i = ith_cluster_silhouette_values.shape[0]
                y_upper = y_lower + size_cluster_i

                color = cm.nipy_spectral(float(i) / n_clusters_to_test)
                ax1.fill_betweenx(np.arange(y_lower, y_upper), 0,
                                  ith_cluster_silhouette_values,
                                  facecolor=color, edgecolor=color, alpha=0.7)

                ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
                y_lower = y_upper + 10

            ax1.set_title("The silhouette plot for the various clusters.")
            ax1.set_xlabel("The silhouette coefficient values")
            ax1.set_ylabel("Cluster label")

            ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

            ax1.set_yticks([])  # Clear the yaxis labels / ticks
            ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

            # 2nd Plot showing the actual clusters formed
            colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters_to_test)
            first_feature = data[:, 0]
            second_feature = data[:, 1]
            ax2.scatter(first_feature, second_feature, marker=".", s=100, lw=0, alpha=0.7,
                        c=colors, edgecolor="k")

            # Labeling the clusters
            centers = clusterer.cluster_centers_
            # Draw white circles at cluster centers
            ax2.scatter(centers[:, 0], centers[:, 1], marker="o", c="white", alpha=1, s=200,
                        edgecolor="k")

            for i, c in enumerate(centers):
                ax2.scatter(c[0], c[1], marker="$%d$" % i, alpha=1, s=50, edgecolor="k")

            ax2.set_title("The visualization of the clustered data.")
            if pcs is None:
                ax2.set_xlabel("Feature space for the 1st feature")
                ax2.set_ylabel("Feature space for the 2nd feature")
            else:
                ax2.set_xlabel("PC1")
                ax2.set_ylabel("PC2")

            plt.suptitle(
                "Silhouette analysis for KMeans clustering on sample data with n_clusters = %d"
                % n_clusters_to_test,
                fontsize=14,
                fontweight="bold",
            )

            plt.show()

    def volcano(self, df, p_value_colname, p_value_thresh, fc_colname, fc_iqr_thresh=1.5,
                top_n=None):
        """
        This function generates a volcano plot of the given DataFrame using matplotlib.

        Args:
            df (pandas.DataFrame): DataFrame to plot
            p_value_thresh (float, optional): P-value threshold for identifying significant points.
                Default is 0.05.

        Returns:
            None
        """
        sig_mask, fc_mask, combined_mask = self.filter_de_df(
            df, fc_colname, p_value_colname, p_value_thresh, fc_iqr_thresh=fc_iqr_thresh)

        plt.figure(figsize=(15, 15))

        # create a scatter plot of the FC vs. -log10(p-value), removing FC outliers
        fc_mask_fc = df.loc[fc_mask, fc_colname]
        fc_mask_p_value = -1 * np.log10(df.loc[fc_mask, p_value_colname])
        plt.scatter(fc_mask_fc, fc_mask_p_value, color='grey', alpha=0.5)

        # create a scatter plot of the FC vs. -log10(p-value), removing FC outliers and
        # keeping only those with significant p-values
        combined_mask_fc = df.loc[combined_mask, fc_colname]
        combined_mask_p = -1 * np.log10(df.loc[combined_mask, p_value_colname])
        plt.scatter(combined_mask_fc, combined_mask_p, color='red', alpha=0.5)

        if top_n is not None:
            selected = top_n
            self._label_volcano(True, combined_mask_fc, combined_mask_p, df, selected)
            self._label_volcano(False, combined_mask_fc, combined_mask_p, df, selected)

        # add labels and title to the plot
        plt.xlabel('Fold Change (log2)')
        plt.ylabel('-log10(P-value)')
        plt.title('Volcano Plot')

        # add a vertical line at FC = 0
        plt.axvline(x=0, color='black', linestyle='--')

        # display the plot
        plt.show()

    def _label_volcano(self, ascending, combined_mask_fc, combined_mask_p, df, top_n):

        # Sort the combined_mask_fc and combined_mask_p series by -log10(P-value)
        sorted_indices = combined_mask_fc.sort_values(ascending=ascending).head(
            top_n).index

        # Get the positions of the top-N indices in the original DataFrame
        top_n_positions = [df.index.get_loc(idx) for idx in sorted_indices]

        # Create an empty list for storing text objects
        texts = []

        # Add labels to the top-N entries
        for position, idx in zip(top_n_positions, sorted_indices):
            label = df.index[position]
            x = combined_mask_fc.loc[idx]
            y = combined_mask_p.loc[idx]
            texts.append(plt.text(x, y, label, fontsize=18))

        # Adjust the labels to avoid overlapping
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', lw=1.0))

    def filter_de_df(self, df, fc_colname, p_value_colname, p_value_thresh,
                     fc_iqr_thresh=1.5):

        # create a boolean mask for significant p-values
        sig_mask = (df[p_value_colname] < p_value_thresh)

        # calculate the IQR and quartiles of the FC values
        q1, q3 = np.percentile(df[fc_colname], [25, 75])
        iqr = q3 - q1

        # define the lower and upper bounds for identifying outliers
        fc_lower_bound = q1 - fc_iqr_thresh * iqr
        fc_upper_bound = q3 + fc_iqr_thresh * iqr

        # filter out the outlier point(s) with low FC value
        fc_mask = (df[fc_colname] > fc_lower_bound) & (df[fc_colname] < fc_upper_bound)
        combined_mask = sig_mask & fc_mask

        return sig_mask, fc_mask, combined_mask


class CaseControlAnalysis():  # FIXME: make this an abstract class
    def __init__(self, data_df, design_df, data_type, case, control):
        self.data_df = data_df
        self.design_df = design_df
        self.data_type = data_type
        self.case = case
        self.control = control
        self.results = None
        self.wi = None

    def __repr__(self):
        return 'CaseControlAnalysis %s: %s_vs_%s' % (self.data_type, self.case, self.control)

    def get_results(self):
        return self.results

    def run(self):
        raise NotImplementedError

    def get_analysis_type(self):
        raise NotImplementedError


class TTestAnalysis(CaseControlAnalysis):
    def __init__(self, data_df, design_df, data_type, case, control):
        super().__init__(data_df, design_df, data_type, case, control)

    def run(self):
        min_replace = SMALL
        self.wi = Inference(self.data_df, self.design_df, self.data_type, min_value=min_replace)
        self.results = self.wi.run_ttest(self.case, self.control)

    def get_analysis_type(self):
        return INFERENCE_T_TEST

    def __repr__(self):
        return 't-test %s: %s_vs_%s' % (self.data_type, self.case, self.control)


class DESeq2Analysis(CaseControlAnalysis):
    def __init__(self, data_df, design_df, data_type, case, control):
        super().__init__(data_df, design_df, data_type, case, control)

    def run(self):
        assert self.data_type == GENES
        MIN_VALUE = 1
        KEEP_THRESHOLD = 10
        REPLACE_MEAN = False
        self.wi = Inference(self.data_df, self.design_df, self.data_type, min_value=MIN_VALUE,
                            replace_mean=REPLACE_MEAN)
        try:
            pd_df, rld_df, res_ordered = self.wi.run_deseq(KEEP_THRESHOLD, self.case, self.control)
            self.results = pd_df[['padj', 'log2FoldChange']]
        except Exception as e:
            logger.warning('Failed to run DESeq2: %s' % str(e))

    def get_analysis_type(self):
        return INFERENCE_DESEQ

    def __repr__(self):
        return 'DESeq2 %s: %s_vs_%s' % (self.data_type, self.case, self.control)


class LimmaAnalysis(CaseControlAnalysis):
    def __init__(self, data_df, design_df, data_type, case, control):
        super().__init__(data_df, design_df, data_type, case, control)

    def run(self):
        min_replace = SMALL
        try:
            self.wi = Inference(self.data_df, self.design_df, self.data_type,
                                min_value=min_replace)
            self.results = self.wi.run_limma(self.case, self.control)
        except Exception as e:
            logger.warning('Failed to run limma: %s' % str(e))

    def get_analysis_type(self):
        return INFERENCE_LIMMA

    def __repr__(self):
        return 'limma %s: %s_vs_%s' % (self.data_type, self.case, self.control)
