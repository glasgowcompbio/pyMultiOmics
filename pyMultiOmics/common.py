import gzip
import logging
import math
import os
import pathlib
import pickle
import sys
import zipfile

from loguru import logger
import requests
from tqdm import tqdm

from pyMultiOmics.constants import TRUNCATE_LIMIT


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


def as_list(value):
    if isinstance(value, str): # if string turn this into a list item
        value = [value]
    else: # other types, check if  they're iterable
        try:
            iter(value)
        except TypeError: # not iterable
            value = [value]
        else: # iterable
            value = list(value)
    return value


def set_log_level(level, remove_id=None):
    if remove_id is None:
        try:
            logger.remove(0)  # try to remove the default handler with id 0
        except ValueError:  # no default handler has been set
            pass
    else:
        logger.remove(remove_id)  # remove previously set handler by id

    # add new handler at the desired log level
    new_handler_id = logger.add(sys.stderr, level=level)
    return new_handler_id


def set_log_level_warning(remove_id=None):
    return set_log_level(logging.WARNING, remove_id=remove_id)


def set_log_level_info(remove_id=None):
    return set_log_level(logging.INFO, remove_id=remove_id)


def set_log_level_debug(remove_id=None):
    return set_log_level(logging.DEBUG, remove_id=remove_id)


def add_log_file(log_path, level):
    logger.add(log_path, level=level)


def truncate(my_str):
    my_str = (my_str[:TRUNCATE_LIMIT] + '...') if len(my_str) > TRUNCATE_LIMIT else my_str
    return my_str
