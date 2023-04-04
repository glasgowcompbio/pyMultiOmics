import os
import sys

sys.path.append('.')
sys.path.append('..')

sys.path.append('/Users/joewandy/Work/git/pyMultiOmics')

from pyMultiOmics.loader import process_affinity_data


def main():
    in_dir = '/Users/joewandy/Library/CloudStorage/OneDrive-UniversityofGlasgow/CiC_Affinity_Biomarkers'
    file_name = os.path.join(in_dir, 'data_group_clustering.xlsx')
    filter_method = 'batch'
    filter_thresh = 0.50
    correct_batch = False
    data_df, sample_metadata_df, feature_metadata_df = process_affinity_data(
        file_name, filter_method=filter_method, filter_thresh=filter_thresh,
        correct_batch=correct_batch
    )


if __name__ == '__main__':
    main()
