from loguru import logger

from .constants import GENES, PROTEINS, COMPOUNDS, REACTIONS, PATHWAYS, SMALL, GENES, INFERENCE_T_TEST, \
    INFERENCE_DESEQ, INFERENCE_LIMMA
from .pipelines import WebOmicsInference


class AnalysisPipeline(object):
    def __init__(self, multi_omics_data, mapping):
        self.multi_omics_data = multi_omics_data
        self.mapping = mapping

        # stores analyses and their results
        self.analyses = {
            GENES: [],
            PROTEINS: [],
            COMPOUNDS: [],
            REACTIONS: [],
            PATHWAYS: []
        }

    def run_de(self, selected_method, data_type, case, control):
        method_to_run = {
            INFERENCE_T_TEST: TTestAnalysis,
            INFERENCE_DESEQ: DESeq2Analysis,
            INFERENCE_LIMMA: LimmaAnalysis
        }
        M = method_to_run[selected_method]

        # check if already exists
        found = False
        for analysis in self.analyses[data_type]:
            if analysis.case == case and analysis.control == control:
                found = True

        if not found:
            data_df, design_df = self.multi_omics_data.get_dfs(data_type)
            if data_df is not None and design_df is not None:
                analysis = M(data_df, design_df, data_type, case, control)
                analysis.run()
                self.analyses[data_type].append(analysis)
        else:
            logger.warning('DE analysis for data_type=%s case=%s control=%s already exists' % (
                data_type, case, control
            ))

    def get_de_analysis(self, data_type, case, control, method=None):
        if method is not None:
            assert method in [INFERENCE_T_TEST, INFERENCE_DESEQ, INFERENCE_LIMMA]

        found = None
        for analysis in self.analyses[data_type]:
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


class CaseControlAnalysis():
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
        self.wi = WebOmicsInference(self.data_df, self.design_df, self.data_type, min_value=min_replace)
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
        self.wi = WebOmicsInference(self.data_df, self.design_df, self.data_type, min_value=MIN_VALUE,
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
            self.wi = WebOmicsInference(self.data_df, self.design_df, self.data_type, min_value=min_replace)
            self.results = self.wi.run_limma(self.case, self.control)
        except Exception as e:
            logger.warning('Failed to run limma: %s' % str(e))

    def get_analysis_type(self):
        return INFERENCE_LIMMA

    def __repr__(self):
        return 'limma %s: %s_vs_%s' % (self.data_type, self.case, self.control)
