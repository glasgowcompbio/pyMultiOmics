import pandas as pd
from loguru import logger

from .constants import GENES, PROTEINS, COMPOUNDS, REACTIONS, PATHWAYS, SMALL, GENOMICS, INFERENCE_T_TEST, \
    INFERENCE_DESEQ, INFERENCE_LIMMA
from .pipelines import WebOmicsInference


class AnalysisPipeline(object):
    def __init__(self, mapping):
        self.mapping = mapping
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
            data_df, design_df = self.mapping.get_dfs(data_type)
            analysis = M(data_df, design_df, data_type, case, control)
            analysis.run()
            self.analyses[data_type].append(analysis)
        else:
            logger.warning('DE analysis for data_type=%s case=%s control=%s already exists' % (
                data_type, case, control
            ))

    def de_results(self, data_type, case, control, method=None):
        if method is not None:
            assert method in [INFERENCE_T_TEST, INFERENCE_DESEQ, INFERENCE_LIMMA]

        dfs = []
        for analysis in self.analyses[data_type]:
            # select only analysis for a particular case and control, if specified
            if analysis.case != case or analysis.control != control:
                continue

            # select method, if specified
            if method is not None and analysis.get_analysis_type() != method:
                continue

            # get results dataframe for the DE analysis
            res = analysis.get_results()

            # rename the column 'log2FoldChange' to 'FC'
            res = res.rename(columns={
                'log2FoldChange': 'FC'
            })

            # add suffix to all columns: 'padj' and 'FC'
            suffix = '_%s_vs_%s' % (analysis.case, analysis.control)
            res = res.add_suffix(suffix)
            dfs.append(res)

        # combine all results
        if len(dfs) > 0:
            df = pd.concat(dfs, axis=1)
            return df
        else:
            return pd.DataFrame()

class CaseControlAnalysis():
    def __init__(self, data_df, design_df, data_type, case, control):
        self.data_df = data_df
        self.design_df = design_df
        self.data_type = data_type
        self.case = case
        self.control = control
        self.results = None

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
        wi = WebOmicsInference(self.data_df, self.design_df, self.data_type, min_value=min_replace)
        self.results = wi.run_ttest(self.case, self.control)

    def get_analysis_type(self):
        return INFERENCE_T_TEST

    def __repr__(self):
        return 't-test %s: %s_vs_%s' % (self.data_type, self.case, self.control)


class DESeq2Analysis(CaseControlAnalysis):
    def __init__(self, data_df, design_df, data_type, case, control):
        super().__init__(data_df, design_df, data_type, case, control)

    def run(self):
        assert self.data_type == GENOMICS
        wi = WebOmicsInference(self.data_df, self.design_df, self.data_type)
        min_replace = 10
        try:
            pd_df, rld_df, res_ordered = wi.run_deseq(min_replace, self.case, self.control)
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
            wi = WebOmicsInference(self.data_df, self.design_df, self.data_type, min_value=min_replace)
            self.results = wi.run_limma(self.case, self.control)
        except Exception as e:
            logger.warning('Failed to run limma: %s' % str(e))

    def get_analysis_type(self):
        return INFERENCE_LIMMA

    def __repr__(self):
        return 'limma %s: %s_vs_%s' % (self.data_type, self.case, self.control)
