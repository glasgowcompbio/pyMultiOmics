import pandas as pd
import numpy as np

from pyMultiOmics.constants import SMALL


def load_affinity_data(file_name):
    """
    Load CIC Affinity data
    Args:
        file_name: the excel file to load

    Returns: a tuple of data df, sample metadata df, and feature metadata df

    """

    # Load the excel file skipping the first 5 header rows
    df = pd.read_excel(file_name, skiprows=5)

    # Rename the 'Group' column to 'group', if it exists
    if 'Group' in df.columns:
        df.rename(columns={'Group': 'group'}, inplace=True)

    # Store Plate number and Group columns separately
    sample_metadata_df = df[['Plate number', 'group']]
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

    # Replace NaN values with a small value
    df.fillna(SMALL, inplace=True)

    # Log the values
    df = np.log(df)

    sample_metadata_df.index.name = 'sample'
    sample_metadata_df.index = sample_metadata_df.index.astype(str)

    return df, sample_metadata_df, feature_metadata_df
