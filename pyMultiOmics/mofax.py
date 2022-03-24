from typing import Union

import matplotlib.pyplot as plt
import mofax as mfx
from mofapy2.run.entry_point import entry_point


class MofaPipeline():
    def __init__(self, MultiOmicsData_obj: object = None, modelPath: str = None):
        self.MultiOmicsData_obj = MultiOmicsData_obj
        self.data = None
        self.model = None
        self.filepath = modelPath
        self.mofa = None

    def set_data(self):
        self.data = self.MultiOmicsData_obj.to_mofa()

    def load_mofa(self, modelPath):
        self.filepath = modelPath

    def training(self,
                 scale_groups: bool = False,
                 scale_views: bool = False,
                 nFactors: int = 10,
                 spikeslab_weights: bool = True,
                 ard_factors: bool = True,
                 ard_weights: bool = True,
                 nIter: int = 100,
                 convergence_mode: str = "fast",
                 startELBO: int = 1,
                 freqELBO: int = 1,
                 dropR2=None,
                 gpu_mode: bool = False,
                 verbose: bool = False,
                 seed: int = 42):

        ## (1) initialise the entry point ##
        ent = entry_point()

        ## (2) Set data options ##
        # - scale_groups: if groups have significantly different ranges, it is good practice to scale each group to unit variance
        # - scale_views: if views have significantly different ranges, it is good practice to scale each view to unit variance
        ent.set_data_options(
            scale_groups=scale_groups,
            scale_views=scale_views)

        self.data.drop_duplicates(subset=["group", "view", "feature", "sample"],
                                  keep='first', inplace=True)
        ent.set_data_df(self.data)

        # Advanced (using personalised values)
        ent.set_model_options(
            factors=nFactors,
            spikeslab_weights=spikeslab_weights,
            ard_factors=ard_factors,
            ard_weights=ard_weights)

        # Advanced (using personalised values)
        ent.set_train_options(
            iter=nIter,
            convergence_mode=convergence_mode,
            startELBO=startELBO,
            freqELBO=freqELBO,
            dropR2=dropR2,
            gpu_mode=gpu_mode,
            verbose=verbose,
            seed=seed)

        ####################################
        ## Build and train the MOFA model ##
        ####################################

        # Build the model 
        ent.build()

        # Run the model
        ent.run()

        self.model = ent

    def save_model(self):
        self.model.save(self.filepath, save_data=True)

    def build_mofa(self):
        self.mofa = mfx.mofa_model(self.filepath)

    def plot_top_features(self,
                          factors: int = None,
                          views: Union[str, int] = None,
                          n_features: int = None,
                          clip_threshold: float = None,
                          scale: bool = False,
                          absolute_values: bool = False,
                          only_positive: bool = False,
                          only_negative: bool = False,
                          per_view: bool = True):

        df = self.mofa.get_top_features(factors=factors,
                                        views=views,
                                        n_features=n_features,
                                        clip_threshold=clip_threshold,
                                        scale=scale,
                                        absolute_values=absolute_values,
                                        only_positive=only_positive,
                                        only_negative=only_negative,
                                        per_view=per_view,
                                        df=True)

        for l in range(len(df['view'].tolist())):
            if df['view'][l] != views:
                df = df.drop(labels=l)

        fig, ax = plt.subplots()

        ax.axes.set_xlabel('Weights')
        ax.axes.set_ylabel('Features')
        title = views + " factor " + str(factors + 1)
        ax.axes.set_title(title)

        ax.axes.axvline(x=0, linewidth=0.5, color='grey', linestyle='-.')

        data = df['value'].tolist()
        data.reverse()
        label = df['feature'].tolist()
        label.reverse()

        b = ax.barh(range(len(df['value'].tolist())), data, tick_label=label, color='black',
                    height=0.05)
        ax.scatter(data, y=label, color='black', s=40)

        # plt.show()

        return (fig)
