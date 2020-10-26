import gzip
import logging
import math
import os
import pathlib
import pickle
import zipfile

import loguru as logger
import requests
from tqdm import tqdm


def create_if_not_exist(out_dir):
    if not os.path.exists(out_dir) and len(out_dir) > 0:
        logger.info('Created %s' % out_dir)
        pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)


def save_obj(obj, filename):
    """
    Save object to file
    :param obj: the object to save
    :param filename: the output file
    :return: None
    """
    out_dir = os.path.dirname(filename)
    create_if_not_exist(out_dir)
    logger.info('Saving %s to %s' % (type(obj), filename))
    with gzip.GzipFile(filename, 'w') as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)


def load_obj(filename):
    """
    Load saved object from file
    :param filename: The file to load
    :return: the loaded object
    """
    try:
        with gzip.GzipFile(filename, 'rb') as f:
            return pickle.load(f)
    except OSError:
        logging.getLogger().warning('Old, invalid or missing pickle in %s. Please regenerate this file.' % filename)
        return None


def download_file(url, out_file=None):
    r = requests.get(url, stream=True)
    total_size = int(r.headers.get('content-length', 0));
    block_size = 1024
    current_size = 0

    if out_file is None:
        out_file = url.rsplit('/', 1)[-1]  # get the last part in url
    logger.info('Downloading %s' % out_file)

    with open(out_file, 'wb') as f:
        for data in tqdm(r.iter_content(block_size), total=math.ceil(total_size // block_size), unit='KB',
                         unit_scale=True):
            current_size += len(data)
            f.write(data)
    assert current_size == total_size
    return out_file


def extract_zip_file(in_file, delete=True):
    logger.info('Extracting %s' % in_file)
    with zipfile.ZipFile(file=in_file) as zip_file:
        for file in tqdm(iterable=zip_file.namelist(), total=len(zip_file.namelist())):
            zip_file.extract(member=file)

    if delete:
        logger.info('Deleting %s' % in_file)
        os.remove(in_file)
