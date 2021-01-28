import collections
import json
import re
from collections import defaultdict
from io import StringIO

import numpy as np
import pandas as pd
from loguru import logger

from .common import load_obj
from .constants import COMPOUND_DATABASE_CHEBI, EXTERNAL_KEGG_TO_CHEBI, GENE_PK, PROTEIN_PK, REACTION_PK, \
    COMPOUND_PK, PATHWAY_PK, COMPOUND_DATABASE_KEGG, NA, EXTERNAL_GENE_NAMES, EXTERNAL_COMPOUND_NAMES, GENOMICS, \
    PROTEOMICS, METABOLOMICS, REACTIONS, PATHWAYS, GENES_TO_PROTEINS, PROTEINS_TO_REACTIONS, COMPOUNDS_TO_REACTIONS, \
    REACTIONS_TO_PATHWAYS, IDENTIFIER_COL, PIMP_PEAK_ID_COL, GROUP_COL, DEFAULT_GROUP_NAME, PADJ_COL_PREFIX, \
    FC_COL_PREFIX, SAMPLE_COL
from .metadata import get_gene_names, get_compound_metadata, clean_label
from .reactome import ensembl_to_uniprot, uniprot_to_reaction, compound_to_reaction, \
    reaction_to_pathway, reaction_to_uniprot, reaction_to_compound, uniprot_to_ensembl

Relation = collections.namedtuple('Relation', 'keys values mapping_list')

from pyMultiOmics.constants import IDENTIFIER_COL


def prepare_input(data_df):
    if data_df is None:
        return None

    # remove index, since all our old codes below assume this
    data_df = data_df.reset_index()

    # always set the first column to be the Identifier column
    rename = {data_df.columns.values[0]: IDENTIFIER_COL}
    data_df = data_df.rename(columns=rename)

    # assume id is in the first column and is a string
    data_df.iloc[:, 0] = data_df.iloc[:, 0].astype(str)
    return data_df


def reactome_mapping(observed_gene_df, observed_protein_df, observed_compound_df,
                     compound_database_str, species_list, metabolic_pathway_only):

    observed_gene_df = prepare_input(observed_gene_df)
    observed_protein_df = prepare_input(observed_protein_df)
    observed_compound_df = prepare_input(observed_compound_df)

    ### all the ids that we have from the user ###
    observed_gene_ids = get_ids_from_dataframe(observed_gene_df)
    observed_protein_ids = get_ids_from_dataframe(observed_protein_df)
    observed_compound_ids = get_ids_from_dataframe(observed_compound_df)

    # try to convert all kegg ids to chebi ids, if possible
    KEGG_2_CHEBI = load_obj(EXTERNAL_KEGG_TO_CHEBI)
    not_found = []
    for cid in observed_compound_ids:
        if cid not in KEGG_2_CHEBI:
            not_found.append(cid)
            KEGG_2_CHEBI[cid] = cid

    if len(not_found) > 0:
        logger.debug('Not found: %d when converting kegg -> chebi ids' % len(not_found))

    if observed_compound_df is not None:
        if compound_database_str == COMPOUND_DATABASE_CHEBI:
            observed_compound_df.iloc[:, 0] = observed_compound_df.iloc[:, 0].map(
                KEGG_2_CHEBI)  # assume 1st column is id
        observed_compound_ids = get_ids_from_dataframe(observed_compound_df)

    ### map genes -> proteins ###
    logger.info('Mapping genes -> proteins')
    gene_2_proteins_mapping, _ = ensembl_to_uniprot(observed_gene_ids, species_list)
    gene_2_proteins = make_relations(gene_2_proteins_mapping, GENE_PK, PROTEIN_PK, value_key=None)

    ### maps proteins -> reactions ###
    logger.info('Mapping proteins -> reactions')
    protein_ids_from_genes = gene_2_proteins.values
    known_protein_ids = list(set(observed_protein_ids + protein_ids_from_genes))
    protein_2_reactions_mapping, _ = uniprot_to_reaction(known_protein_ids, species_list)
    protein_2_reactions = make_relations(protein_2_reactions_mapping, PROTEIN_PK, REACTION_PK,
                                         value_key='reaction_id')

    ### maps compounds -> reactions ###
    logger.info('Mapping compounds -> reactions')
    compound_2_reactions_mapping, _ = compound_to_reaction(observed_compound_ids, species_list)
    compound_2_reactions = make_relations(compound_2_reactions_mapping, COMPOUND_PK, REACTION_PK,
                                          value_key='reaction_id')

    ### maps reactions -> metabolite pathways ###
    logger.info('Mapping reactions -> pathways')
    reaction_ids_from_proteins = protein_2_reactions.values
    reaction_ids_from_compounds = compound_2_reactions.values
    reaction_ids = list(set(reaction_ids_from_proteins + reaction_ids_from_compounds))
    reaction_2_pathways_mapping, reaction_2_pathways_id_to_names = reaction_to_pathway(reaction_ids,
                                                                                       species_list,
                                                                                       metabolic_pathway_only)
    reaction_2_pathways = make_relations(reaction_2_pathways_mapping, REACTION_PK, PATHWAY_PK,
                                         value_key='pathway_id')

    ### maps reactions -> proteins ###
    logger.info('Mapping reactions -> proteins')
    mapping, _ = reaction_to_uniprot(reaction_ids, species_list)
    reaction_2_proteins = make_relations(mapping, REACTION_PK, PROTEIN_PK, value_key=None)
    protein_2_reactions = merge_relation(protein_2_reactions, reverse_relation(reaction_2_proteins))
    all_protein_ids = protein_2_reactions.keys

    ### maps reactions -> compounds ###
    logger.info('Mapping reactions -> compounds')
    if compound_database_str == COMPOUND_DATABASE_KEGG:
        use_kegg = True
    else:
        use_kegg = False
    reaction_2_compounds_mapping, reaction_to_compound_id_to_names = reaction_to_compound(reaction_ids, species_list,
                                                                                          use_kegg)
    reaction_2_compounds = make_relations(reaction_2_compounds_mapping, REACTION_PK, COMPOUND_PK, value_key=None)
    compound_2_reactions = merge_relation(compound_2_reactions, reverse_relation(reaction_2_compounds))
    all_compound_ids = compound_2_reactions.keys

    ### map proteins -> genes ###
    logger.info('Mapping proteins -> genes')
    mapping, _ = uniprot_to_ensembl(all_protein_ids, species_list)
    protein_2_genes = make_relations(mapping, PROTEIN_PK, GENE_PK, value_key=None)
    gene_2_proteins = merge_relation(gene_2_proteins, reverse_relation(protein_2_genes))
    all_gene_ids = gene_2_proteins.keys

    ### add links ###

    # map NA to NA
    gene_2_proteins = add_links(gene_2_proteins, GENE_PK, PROTEIN_PK, [NA], [NA])
    protein_2_reactions = add_links(protein_2_reactions, PROTEIN_PK, REACTION_PK, [NA], [NA])
    compound_2_reactions = add_links(compound_2_reactions, COMPOUND_PK, REACTION_PK, [NA], [NA])
    reaction_2_pathways = add_links(reaction_2_pathways, REACTION_PK, PATHWAY_PK, [NA], [NA])

    # map genes that have no proteins to NA
    gene_pk_list = [x for x in all_gene_ids if x not in gene_2_proteins.keys]
    gene_2_proteins = add_links(gene_2_proteins, GENE_PK, PROTEIN_PK, gene_pk_list, [NA])

    # map proteins that have no genes to NA
    protein_pk_list = [x for x in all_protein_ids if x not in gene_2_proteins.values]
    gene_2_proteins = add_links(gene_2_proteins, GENE_PK, PROTEIN_PK, [NA], protein_pk_list)

    # map proteins that have no reactions to NA
    protein_pk_list = [x for x in all_protein_ids if x not in protein_2_reactions.keys]
    protein_2_reactions = add_links(protein_2_reactions, PROTEIN_PK, REACTION_PK, protein_pk_list, [NA])

    # map reactions that have no proteins to NA
    reaction_pk_list = [x for x in reaction_ids if x not in protein_2_reactions.values]
    protein_2_reactions = add_links(protein_2_reactions, PROTEIN_PK, REACTION_PK, [NA], reaction_pk_list)

    # map compounds that have no reactions to NA
    compound_pk_list = [x for x in all_compound_ids if x not in compound_2_reactions.keys]
    compound_2_reactions = add_links(compound_2_reactions, COMPOUND_PK, REACTION_PK, compound_pk_list, [NA])

    # map reactions that have no compounds to NA
    reaction_pk_list = [x for x in reaction_ids if x not in compound_2_reactions.values]
    compound_2_reactions = add_links(compound_2_reactions, COMPOUND_PK, REACTION_PK, [NA], reaction_pk_list)

    # map reactions that have no pathways to NA
    reaction_pk_list = [x for x in reaction_ids if x not in reaction_2_pathways.keys]
    reaction_2_pathways = add_links(reaction_2_pathways, REACTION_PK, PATHWAY_PK, reaction_pk_list, [NA])

    GTF_DICT = load_obj(EXTERNAL_GENE_NAMES)
    metadata_map = get_gene_names(all_gene_ids, GTF_DICT)
    genes_json = pk_to_json(GENE_PK, 'gene_id', all_gene_ids, metadata_map, observed_gene_df,
                            observed_ids=observed_gene_ids)
    gene_2_proteins_json = json.dumps(gene_2_proteins.mapping_list)

    # metadata_map = get_uniprot_metadata_online(uniprot_ids)
    proteins_json = pk_to_json('protein_pk', 'protein_id', all_protein_ids, metadata_map, observed_protein_df,
                               observed_ids=observed_protein_ids)
    protein_2_reactions_json = json.dumps(protein_2_reactions.mapping_list)

    # TODO: this feels like a very bad way to implement this
    # We need to deal with uploaded peak data from PiMP, which contains a lot of duplicate identifications per peak
    KEGG_ID_2_DISPLAY_NAMES = load_obj(EXTERNAL_COMPOUND_NAMES)
    metadata_map = get_compound_metadata(all_compound_ids, KEGG_ID_2_DISPLAY_NAMES, reaction_to_compound_id_to_names)
    try:
        mapping = get_mapping(observed_compound_df)
    except KeyError:
        mapping = None
    except AttributeError:
        mapping = None
    compounds_json = pk_to_json('compound_pk', 'compound_id', all_compound_ids, metadata_map, observed_compound_df,
                                observed_ids=observed_compound_ids, mapping=mapping)
    if mapping:
        compound_2_reactions = expand_relation(compound_2_reactions, mapping, 'compound_pk')
    compound_2_reactions_json = json.dumps(compound_2_reactions.mapping_list)

    metadata_map = {}
    for name in reaction_2_pathways_id_to_names:
        tok = reaction_2_pathways_id_to_names[name]['name']
        filtered = clean_label(tok)
        species = reaction_2_pathways_id_to_names[name]['species']
        metadata_map[name] = {'display_name': filtered, 'species': species}

    reaction_count_df = None
    pathway_count_df = None

    pathway_ids = reaction_2_pathways.values
    reactions_json = pk_to_json('reaction_pk', 'reaction_id', reaction_ids, metadata_map, reaction_count_df,
                                has_species=True)
    pathways_json = pk_to_json('pathway_pk', 'pathway_id', pathway_ids, metadata_map, pathway_count_df,
                               has_species=True)
    reaction_2_pathways_json = json.dumps(reaction_2_pathways.mapping_list)

    results = {
        GENOMICS: genes_json,
        PROTEOMICS: proteins_json,
        METABOLOMICS: compounds_json,
        REACTIONS: reactions_json,
        PATHWAYS: pathways_json,
        GENES_TO_PROTEINS: gene_2_proteins_json,
        PROTEINS_TO_REACTIONS: protein_2_reactions_json,
        COMPOUNDS_TO_REACTIONS: compound_2_reactions_json,
        REACTIONS_TO_PATHWAYS: reaction_2_pathways_json,
    }
    return results


def get_mapping(observed_compound_df):
    mapping = defaultdict(list)
    for idx, row in observed_compound_df.iterrows():
        identifier = row[IDENTIFIER_COL]
        peak_id = row[PIMP_PEAK_ID_COL]
        mapping[identifier].append('%s_%s' % (identifier, peak_id))
    return dict(mapping)


def csv_to_dataframe(csv_str):
    # extract group, if any
    filtered_str = ''
    group_str = None
    for line in csv_str.splitlines():  # go through all lines and remove the line containing the grouping info
        if re.match(GROUP_COL, line, re.I):
            group_str = line
        else:
            filtered_str += line + '\n'

    # extract id values
    data = StringIO(filtered_str)
    try:
        data_df = pd.read_csv(data)
        data_df.columns = data_df.columns.str.replace('.',
                                                      '_')  # replace period with underscore to prevent alasql breaking
        data_df.columns = data_df.columns.str.replace('-',
                                                      '_')  # replace dash with underscore to prevent alasql breaking
        data_df.columns = data_df.columns.str.replace('#', '')  # remove funny characters
        rename = {data_df.columns.values[0]: IDENTIFIER_COL}
        for i in range(len(data_df.columns.values[1:])):  # sql doesn't like column names starting with a number
            col_name = data_df.columns.values[i]
            if col_name[0].isdigit():
                new_col_name = '_' + col_name  # append an underscore in front of the column name
                rename[col_name] = new_col_name
        data_df = data_df.rename(columns=rename)
        data_df.iloc[:, 0] = data_df.iloc[:, 0].astype(str)  # assume id is in the first column and is a string
        id_list = data_df.iloc[:, 0].values.tolist()
    except pd.errors.EmptyDataError:
        data_df = None
        id_list = []

    # create grouping dataframe
    group_df = None
    if data_df is not None:
        sample_data = data_df.columns.values
        if group_str is not None:
            group_data = group_str.split(',')
        else:
            num_samples = len(sample_data)
            group_data = [DEFAULT_GROUP_NAME for x in
                          range(num_samples)]  # assigns a default group if nothing specified

        # skip non-measurement columns
        filtered_sample_data = []
        filtered_group_data = []
        for i in range(len(sample_data)):
            sample_name = sample_data[i]
            if sample_name == IDENTIFIER_COL or \
                    sample_name == PIMP_PEAK_ID_COL or \
                    sample_name.startswith(PADJ_COL_PREFIX) or \
                    sample_name.startswith(FC_COL_PREFIX):
                continue
            filtered_sample_data.append(sample_data[i])
            filtered_group_data.append(group_data[i])

        # convert to dataframe
        if len(filtered_group_data) > 0:
            group_df = pd.DataFrame(list(zip(filtered_sample_data, filtered_group_data)),
                                    columns=[SAMPLE_COL, GROUP_COL])

    return data_df, group_df, id_list


def get_ids_from_dataframe(df):
    if df is None:
        return []
    else:
        return df.iloc[:, 0].values.tolist()  # id is always the first column


def merge_relation(r1, r2):
    unique_keys = list(set(r1.keys + r2.keys))
    unique_values = list(set(r1.values + r2.values))
    mapping_list = r1.mapping_list + r2.mapping_list
    mapping_list = list(map(dict, set(map(lambda x: frozenset(x.items()), mapping_list))))  # removes duplicates, if any
    return Relation(keys=list(unique_keys), values=list(unique_values),
                    mapping_list=mapping_list)


def reverse_relation(rel):
    return Relation(keys=rel.values, values=rel.keys, mapping_list=rel.mapping_list)


def expand_relation(rel, mapping, pk_col):
    expanded_keys = substitute(rel.keys, mapping)
    expanded_values = substitute(rel.values, mapping)
    expanded_mapping_list = []
    for row in rel.mapping_list:
        expanded = expand_each(row, mapping, pk_col)
        if len(expanded) == 0:
            expanded = [row]
        expanded_mapping_list.extend(expanded)
    return Relation(keys=expanded_keys, values=expanded_values, mapping_list=expanded_mapping_list)


def substitute(my_list, mapping):
    new_list = []
    for x in my_list:
        if x in mapping:
            new_list.extend(mapping[x])
        else:
            new_list.append(x)
    return new_list


def expand_each(row, mapping, pk_col):
    results = []
    pk = row[pk_col]
    try:
        replacements = mapping[pk]
        for rep in replacements:
            new_row = without_keys(row, [pk_col])
            new_row[pk_col] = rep
            results.append(new_row)
    except KeyError:
        pass
    return results


# https://stackoverflow.com/questions/31433989/return-copy-of-dictionary-excluding-specified-keys
def without_keys(d, keys):
    return {x: d[x] for x in d if x not in keys}


def pk_to_json(pk_label, display_label, data, metadata_map, observed_df, has_species=False,
               observed_ids=None, mapping=None):
    if observed_df is not None:
        if PIMP_PEAK_ID_COL in observed_df.columns:  # if peak id is present, rename the identifier column to include it
            observed_df[IDENTIFIER_COL] = observed_df[IDENTIFIER_COL] + '_' + observed_df[PIMP_PEAK_ID_COL].astype(str)
        if mapping is not None:
            data = expand_data(data, mapping)
        observed_df = observed_df.set_index(IDENTIFIER_COL)  # set identifier as index
        observed_df = observed_df[~observed_df.index.duplicated(keep='first')]  # remove row with duplicate indices
        observed_df = observed_df.fillna(value=0)  # replace all NaNs with 0s

    output = []
    for item in sorted(data):

        if item == NA:
            continue  # handled below after this loop

        if '_' in item:
            tokens = item.split('_')
            assert len(tokens) == 2
            item = tokens[0]
            peak_id = tokens[1]
        else:
            peak_id = None

        # add observed status and the primary key label to row data
        row = {}
        if observed_ids is not None:
            if item in observed_ids:
                row['obs'] = True
            else:
                row['obs'] = False
        else:
            row['obs'] = None

        if peak_id:
            key = '%s_%s' % (item, peak_id)
            row[pk_label] = key
        else:
            row[pk_label] = item

        # add display label to row_data
        species = None
        if len(metadata_map) > 0 and item in metadata_map and metadata_map[item] is not None:
            if peak_id:
                label = '%s (%s)' % (metadata_map[item]['display_name'].capitalize(), peak_id)
            else:
                label = metadata_map[item]['display_name'].capitalize()
            # get the species if it's there too
            if has_species and 'species' in metadata_map[item]:
                species = metadata_map[item]['species']
        else:
            label = item  # otherwise use the item id as the label
        row[display_label] = label

        # add the remaining data columns to row
        if observed_df is not None:
            try:
                if peak_id:
                    observed_values = observed_df.loc[key].to_dict()
                else:
                    observed_values = observed_df.loc[item].to_dict()
            except KeyError:  # missing data
                observed_values = {}
                for col in observed_df.columns:
                    observed_values.update({col: None})
            observed_values.pop(PIMP_PEAK_ID_COL, None)  # remove pimp peakid column
            # convert numpy bool to python bool, else json serialisation will break
            for k, v in observed_values.items():
                if type(v) == np.bool_:
                    observed_values[k] = bool(v)
            row.update(observed_values)

        if has_species:
            row['species'] = species

        if row not in output:
            output.append(row)

    # add dummy entry
    row = {'obs': NA, pk_label: NA, display_label: NA}
    if has_species:
        row['species'] = NA

    if observed_df is not None:  # also add the remaining columns
        for col in observed_df.columns:
            if col == PIMP_PEAK_ID_COL:
                continue
            row.update({col: 0})

    if row not in output:
        output.append(row)

    output_json = json.dumps(output)
    return output_json


def expand_data(data, mapping):
    new_data = []
    for x in data:
        if x in mapping:
            new_data.extend(mapping[x])
        else:
            new_data.append(x)
    data = new_data
    return data


def make_relations(mapping, source_pk, target_pk, value_key=None):
    id_values = []
    mapping_list = []

    for key in mapping:
        value_list = mapping[key]

        # value_list can be either a list of strings or dictionaries
        # check if the first element is a dict, else assume it's a string
        assert len(value_list) > 0
        is_string = True
        first = value_list[0]
        if isinstance(first, dict):
            is_string = False

        # process each element in value_list
        for value in value_list:
            if is_string:  # value_list is a list of string
                actual_value = value
            else:  # value_list is a list of dicts
                assert value_key is not None, 'value_key is missing'
                actual_value = value[value_key]
            id_values.append(actual_value)
            row = {source_pk: key, target_pk: actual_value}
            mapping_list.append(row)

    unique_keys = set(mapping.keys())
    unique_values = set(id_values)

    return Relation(keys=list(unique_keys), values=list(unique_values),
                    mapping_list=mapping_list)


def add_dummy(relation, source_ids, target_ids, source_pk_label, target_pk_label):
    to_add = [x for x in source_ids if x not in relation.keys]
    relation = add_links(relation, source_pk_label, target_pk_label, to_add, [NA])

    # to_add = [x for x in target_ids if x not in relation.values]
    # relation = add_links(relation, source_pk_label, target_pk_label, [NA], to_add)

    # relation = add_links(relation, source_pk_label, target_pk_label, [NA], [NA])
    return relation


def add_links(relation, source_pk_label, target_pk_label, source_pk_list, target_pk_list):
    rel_mapping_list = list(relation.mapping_list)
    rel_keys = relation.keys
    rel_values = relation.values

    for s1 in source_pk_list:
        if s1 not in rel_keys: rel_keys.append(s1)
        for s2 in target_pk_list:
            rel_mapping_list.append({source_pk_label: s1, target_pk_label: s2})
            if s2 not in rel_keys: rel_values.append(s2)

    return Relation(keys=rel_keys, values=rel_values, mapping_list=rel_mapping_list)
