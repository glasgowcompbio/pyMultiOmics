import json
import urllib.request

from bioservices import ChEBI
from bioservices import EUtils
from bioservices import Ensembl
from bioservices import Reactome
from bioservices import UniProt
from bioservices.kegg import KEGG
from loguru import logger
from tqdm import tqdm


################################################################################
### Gene-related functions                                                   ###
################################################################################


def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]


def get_ensembl_metadata_online(ensembl_ids):

    ensembl_ids = list(set(ensembl_ids))
    logger.info('get_ensembl_metadata %d' % len(ensembl_ids))

    BATCH_SIZE = 1000
    ens = Ensembl()
    ensembl_lookup = {}
    cumulative_total = 0
    for x in batch(ensembl_ids, BATCH_SIZE):
        batch_ids = [i for i in x]
        cumulative_total += len(batch_ids)
        logger.debug('%d/%d' % (cumulative_total, len(ensembl_ids)))
        lookup = ens.post_lookup_by_id(identifiers=batch_ids)
        ensembl_lookup.update(lookup)

    return ensembl_lookup


def get_single_ensembl_metadata_online(ensembl_id):
    ens = Ensembl()
    res = ens.get_lookup_by_id(ensembl_id, expand=True)
    return res


def get_entrez_summary(gene_id):
    e = EUtils()
    ret = e.ESearch('gene', gene_id)
    ret = e.EFetch('gene', ret['idlist'])
    return ret.decode("utf-8")


def get_gene_names(ensembl_ids, gtf_dict):
    metadata_map = {}
    for ensembl_id in ensembl_ids:
        try:
            display_name = gtf_dict[ensembl_id]
            display_name = clean_label(display_name)
            metadata_map[ensembl_id] = {'display_name': display_name}
        except KeyError:
            metadata_map[ensembl_id] = {'display_name': ensembl_id}
    return metadata_map


################################################################################
### Protein-related functions                                                ###
################################################################################


def get_uniprot_metadata_online(uniprot_ids):

    uniprot_ids = list(set(uniprot_ids))
    logger.info('get_uniprot_metadata %d' % len(uniprot_ids))

    BATCH_SIZE = 200
    uniprot = UniProt()
    uniprot_lookup = {}

    cumulative_total = 0
    for x in batch(uniprot_ids, BATCH_SIZE):
        batch_ids = [i for i in x]
        cumulative_total += len(batch_ids)
        logger.debug('%d/%d' % (cumulative_total, len(uniprot_ids)))

        res = uniprot.retrieve(batch_ids)
        for r in res:
            for key in r['accession']:
                protein_id = key.contents[0]
                for x in r['recommendedname']:
                    tag = x.find('shortname')
                    if tag is None:
                        tag = x.find('fullname')
                    label = tag.contents[0]
                    uniprot_lookup[protein_id] = {'display_name': label}

    return uniprot_lookup


def get_single_uniprot_metadata_online(uniprot_id):
    uniprot = UniProt()
    res = uniprot.retrieve(uniprot_id)
    return res


def clean_label(w):
    results = []
    for tok in w.split(' '):
        # filtered = re.sub(r'[^\w\s]', '', tok)
        filtered = tok.replace("'", "")
        filtered = filtered.replace('"', "")
        results.append(filtered.strip())
    return ' '.join(results)


def get_uniprot_metadata_reactome(uniprot_ids):
    protein_descriptions = {}  # TODO: get the description of each protein from reactome
    metadata_map = {}
    for protein_id in protein_descriptions:
        metadata_map[protein_id] = {'display_name': protein_id}
        if protein_id in protein_descriptions and protein_descriptions[protein_id] is not None:
            desc = protein_descriptions[protein_id]
            tokens = desc.split(':')
            for i in range(len(tokens)):
                w = tokens[i]
                if w.startswith('recommendedName'):
                    next_w = clean_label(tokens[i + 1])
                    metadata_map[protein_id] = {'display_name': next_w}
                    logger.debug('%d -- %s' % (protein_id, next_w))
    return metadata_map


################################################################################
### Compound-related functions                                               ###
################################################################################

def kegg_to_chebi(compound_ids):
    results = {}
    ch = ChEBI()
    for compound_id in compound_ids:
        try:
            if compound_id.startswith('C'):
                res = ch.getLiteEntity(compound_id)
                assert len(res) > 0 # we should always be able to convert from KEGG -> ChEBI ids
                le = res[0]
                chebi_number = le.chebiId.split(':')[1]
                logger.debug('KEGG %s --> ChEBI %s' % (compound_id, chebi_number))
                results[compound_id] = chebi_number
            else:
                logger.debug('ChEBI %s --> ChEBI %s' % (compound_id, compound_id))
                results[compound_id] = compound_id
        except Exception:
            logger.debug('KEGG %s --> KEGG %s' % (compound_id, compound_id))
            results[compound_id] = compound_id
    return results


def get_compound_metadata_online(kegg_ids):
    s = KEGG()
    metadata_map = {}
    my_list = list(range(len(kegg_ids)))
    failed = 0
    for i in tqdm(my_list):
        try:
            kegg_id = kegg_ids[i]
            res = s.get(kegg_id)
            d = s.parse(res)
            first_name = d['NAME'][0]
            first_name = first_name.replace(';', '') # strip last ';' character
            metadata_map[kegg_id] = {'display_name': first_name}
        except TypeError:
            failed += 1
    logger.info('Failed = %d' % failed)

    return metadata_map


def get_compound_metadata(compound_ids, kegg_id_to_display_names, id_to_names):
    metadata_map = {}
    for compound_id in compound_ids:
        try:
            display_name = clean_label(kegg_id_to_display_names[compound_id]['display_name'])
        except KeyError:
            try:
                display_name = clean_label(id_to_names[compound_id])
                idx = display_name.find('[ChEBI') # remove [ChEBI ...]
                display_name = display_name[0:idx].strip()
            except KeyError:
                display_name = compound_id
        metadata_map[compound_id] = {'display_name': display_name}
    return metadata_map


def get_single_compound_metadata_online(compound_id):

    if compound_id.upper().startswith('C'):
        s = KEGG()
        res = s.get(compound_id)
        return s.parse(res)
    else:
        ch = ChEBI()
        res = ch.getCompleteEntity('CHEBI:'+compound_id)
        return res


################################################################################
### Reaction-related functions                                               ###
################################################################################


def get_reactome_content_service(reactome_id):
    json_url = 'https://reactome.org/ContentService/data/query/' + reactome_id
    with urllib.request.urlopen(json_url) as url:
        lookup = json.loads(url.read().decode())
        return lookup


################################################################################
### Pathway-related functions                                                ###
################################################################################


################################################################################
### Species-related functions                                                ###
################################################################################


def get_species_name_to_id():
    r = Reactome()
    species_name_to_id = {x['displayName']: x['dbId'] for x in r.data_species_all()}
    return species_name_to_id