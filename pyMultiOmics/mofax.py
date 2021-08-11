from mofapy2.run.entry_point import entry_point
import pandas as pd
import mofax as mfx
from typing import Union


def data_processing(dict_of_files):
    data = pd.DataFrame()
    for k in dict_of_files.keys():
        df = pd.read_csv(dict_of_files[k], sep="\t")
        df = df.rename(columns={df.columns[0]:'sample'})
        for col in df.columns:
            if col != 'sample':
                if isinstance(df[col][1], np.float64) == False:
                    df = df.drop(labels=col,axis=1)
        df = df.melt(id_vars=df.columns[0],var_name = 'feature')
        df['view'] = k
        data = data.append(df)
    return(data)

def set_mofax_options():
	## (1) initialise the entry point ##
	ent = entry_point()


	## (2) Set data options ##
	# - scale_groups: if groups have significantly different ranges, it is good practice to scale each group to unit variance
	# - scale_views: if views have significantly different ranges, it is good practice to scale each view to unit variance
	ent.set_data_options(
		scale_groups = False, 
		scale_views = False
	)

	ent.set_data_df(data)

	# Simple (using default values)
	ent.set_model_options()

	# Advanced (using personalised values)
	ent.set_model_options(
		factors = 5, 
		spikeslab_weights = True, 
		ard_factors = True, 
		ard_weights = True
	)


	# Simple (using default values)
	ent.set_train_options()

	# Advanced (using personalised values)
	ent.set_train_options(
		iter = 100, 
		convergence_mode = "fast", 
		startELBO = 1, 
		freqELBO = 1, 
		dropR2 = None, 
		gpu_mode = False, 
		verbose = False, 
		seed = 42
	)

	####################################
	## Build and train the MOFA model ##
	####################################

	# Build the model 
	ent.build()

	# Run the model
	ent.run()

	return(ent)


def plot_top_features(
    m,
    factors: Union[int] = None,
    views: Union[str, int] = None,
    n_features: int = None,
    clip_threshold: float = None,
    scale: bool = False,
    absolute_values: bool = False,
    only_positive: bool = False,
    only_negative: bool = False,
    per_view: bool = True
):

    df = m.get_top_features(
        factors = factors,
        views = views,
        n_features = n_features,
        clip_threshold = clip_threshold,
        scale = scale,
        absolute_values = absolute_values,
        only_positive = only_positive,
        only_negative = only_negative,
        per_view = per_view,
        df = True,
    )
    
    
    for l in range(len(df['view'].tolist())):
        if df['view'][l] != views:
            df = df.drop(labels = l)
    
    fig, ax = plt.subplots()
    
    ax.axes.set_xlabel('Weights')
    ax.axes.set_ylabel('Features')
    title = views + " factor " + str(factors)
    ax.axes.set_title(title)
    
    ax.axes.axvline(x=0,linewidth=0.5, color='grey',linestyle = '-.')
    
    data = df['value'].tolist()
    data.reverse()
    label = df['feature'].tolist()
    label.reverse()
    
    b = ax.barh(range(len(df['value'].tolist())), data, tick_label=label, color = 'black', height=0.05)
    ax.scatter(data, y=label, color='black', s=40)
    
    #plt.show()
    
    return(fig)








