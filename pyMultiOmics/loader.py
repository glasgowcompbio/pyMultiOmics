import warnings

import numpy as np
import pandas as pd
from combat.pycombat import pycombat


def process_affinity_data(file_name, filter_method='batch', filter_thresh=0.50,
                          correct_batch=False):
    """
    Process CIC Affinity data
    Args:
        file_name: the excel file to load

    Returns: a tuple of data df, sample metadata df, and feature metadata df

    """

    # load from excel
    df, feature_metadata_df, sample_metadata_df = load_affinity_data(file_name)

    # feature filtering
    assert filter_method in ['count', 'batch']
    if filter_method == 'count':
        df, feature_metadata_df = filter_affinity_data_by_count(
            df, feature_metadata_df, filter_thresh)
    elif filter_method == 'batch':
        df, feature_metadata_df = filter_affinity_data_by_batch(
            df, feature_metadata_df, sample_metadata_df, filter_thresh)
    else:
        print('Unknown filtering method')  # shouldn't happen

    # missing value imputation
    df = imput_affinity_data(df)

    # batch correction
    if correct_batch:
        df = batch_correction(df, sample_metadata_df)

    return df, sample_metadata_df, feature_metadata_df


def load_affinity_data(file_name):
    # Load the excel file skipping the first 5 header rows
    df = pd.read_excel(file_name, skiprows=5)

    # Rename the 'Group' column to 'group', if it exists
    if 'Group' in df.columns:
        df.rename(columns={'Group': 'group'}, inplace=True)

    # Store Plate number and Group columns separately
    sample_metadata_df = df[['Plate number', 'group']]
    sample_metadata_df.index.name = 'sample'
    sample_metadata_df.index = sample_metadata_df.index.astype(str)
    sample_metadata_df = sample_metadata_df.astype(str)
    df = df.drop(columns=['Plate number', 'group'])

    # Read the header (feature) columns again into a separate df
    # Set first row as column headers
    # Set the index to be the same as the columns data df
    feature_metadata_df = pd.read_excel(file_name, nrows=5, header=None).transpose()
    feature_metadata_df.columns = feature_metadata_df.iloc[0]
    feature_metadata_df = feature_metadata_df[2:]
    feature_metadata_df.reset_index(drop=True)
    feature_metadata_df.index = df.columns

    # put things in the right shape
    df = df.transpose()
    df.index.name = 'Identifier'
    df.columns = df.columns.astype(str)

    # Convert all columns to float, handling '< LOD' as NaN
    df = df.apply(pd.to_numeric, errors='coerce')
    return df, feature_metadata_df, sample_metadata_df


def filter_affinity_data_by_count(df, feature_metadata_df, filter_thresh):
    # Determine the threshold count as proportion of the number of samples
    thresh_count = int(df.shape[1] * filter_thresh)

    # drop rows with all NaN values below the threshold count in data_df
    filtered_df = df.dropna(thresh=df.shape[1] - thresh_count)

    # filter the feature metadata too
    filtered_feature_metadata_df = pd.merge(feature_metadata_df, filtered_df, how='inner',
                                            left_index=True, right_index=True)
    return filtered_df, filtered_feature_metadata_df


def filter_affinity_data_by_batch(df, feature_metadata_df, sample_metadata_df, filter_thresh):
    # Feature filtering starts here. We want to remove features that ae
    df_transposed = df.T

    # Add the 'Plate number' column from 'sample_metadata_df' to the transposed DataFrame
    df_transposed['Plate number'] = sample_metadata_df['Plate number'].values

    # Initialize an empty set to store the features that meet the criteria
    features_to_keep = set(df.index)

    # Iterate over each unique plate number
    for plate_number in sample_metadata_df['Plate number'].unique():
        # Filter the transposed DataFrame to get samples from the current plate
        plate_samples = df_transposed[df_transposed['Plate number'] == plate_number]

        # Calculate the presence of each feature in more than 50% of the samples for the current plate
        presence_mask = (plate_samples.drop('Plate number', axis=1) > 0).mean() > filter_thresh

        # Get the features that are not present in more than 50% of the samples in the current plate
        features_to_remove = set(presence_mask[~presence_mask].index)

        # Remove these features from the set of features to keep
        features_to_keep = features_to_keep - features_to_remove

    # Filter the original 'df' DataFrame using the features_to_keep set
    filtered_df = df.loc[list(features_to_keep)]

    # filter the feature metadata too
    filtered_feature_metadata_df = pd.merge(feature_metadata_df, filtered_df, how='inner',
                                            left_index=True, right_index=True)
    return filtered_df, filtered_feature_metadata_df


def batch_correction(df, sample_metadata_df):
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    batch = sample_metadata_df['Plate number']
    mod = sample_metadata_df['group'].values
    df_logged = pycombat(np.log(df), batch, mod=mod)
    df = np.exp(df_logged)
    return df


def imput_affinity_data(filtered_df):
    return filtered_df.apply(lambda row: row.fillna(row.min()), axis=1)
