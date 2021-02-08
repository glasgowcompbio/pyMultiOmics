from loguru import logger

from .constants import SMALL, GENOMICS
from .pipelines import WebOmicsInference


class CaseControlAnalysis():
    def __init__(self, data_df, design_df, data_type, case, control):
        self.data_df = data_df
        self.design_df = design_df
        self.data_type = data_type
        self.case = case
        self.control = control
        self.result_df = None

    def run(self):
        raise NotImplementedError

    def get_display_name(self):
        display_name = 'CaseControlAnalysis: %s_vs_%s' % (self.case, self.control)
        return display_name


class TTestAnalysis(CaseControlAnalysis):
    def __init__(self, data_df, design_df, data_type, case, control):
        super().__init__(data_df, design_df, data_type, case, control)

    def run(self):
        min_replace = SMALL
        wi = WebOmicsInference(self.data_df, self.design_df, self.data_type, min_value=min_replace)
        self.result_df = wi.run_ttest(self.case, self.control)
        return self.result_df

    def get_display_name(self):
        display_name = 't-test: %s_vs_%s' % (self.case, self.control)
        return display_name


class DESeq2Analysis(CaseControlAnalysis):
    def __init__(self, data_df, design_df, data_type, case, control):
        super().__init__(data_df, design_df, data_type, case, control)

    def run(self):
        assert self.data_type == GENOMICS
        wi = WebOmicsInference(self.data_df, self.design_df, self.data_type)
        min_replace = 10
        try:
            pd_df, rld_df, res_ordered = wi.run_deseq(min_replace, self.case, self.control)
            self.result_df = pd_df[['padj', 'log2FoldChange']]
            return self.result_df
        except Exception as e:
            logger.warning('Failed to run DESeq2: %s' % str(e))

    def get_display_name(self):
        display_name = 'DESeq2: %s_vs_%s' % (self.case, self.control)
        return display_name


class LimmaAnalysis(CaseControlAnalysis):
    def __init__(self, data_df, design_df, data_type, case, control):
        super().__init__(data_df, design_df, data_type, case, control)

    def run(self):
        min_replace = SMALL
        try:
            wi = WebOmicsInference(self.data_df, self.design_df, self.data_type, min_value=min_replace)
            self.result_df = wi.run_limma(self.case, self.control)
            return self.result_df
        except Exception as e:
            logger.warning('Failed to run limma: %s' % str(e))

    def get_display_name(self):
        display_name = 'limma: %s_vs_%s' % (self.case, self.control)
        return display_name
