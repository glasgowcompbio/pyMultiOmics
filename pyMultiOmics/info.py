import collections
import json
import urllib

from .common import truncate
from .constants import GENES, PROTEINS, COMPOUNDS, REACTIONS, PATHWAYS
from .harmonizomeapi import Harmonizome, Entity
from .metadata import get_single_ensembl_metadata_online, get_single_uniprot_metadata_online, \
    get_single_compound_metadata_online
from .reactome import get_reactome_description, get_reaction_entities, pathway_to_reactions


def get_info(node_id, data_type):
    if data_type == 'genes':
        res = get_ensembl_gene_info(node_id)
    elif data_type == 'proteins':
        res = get_uniprot_protein_info(node_id)
    elif data_type == 'compounds':
        res = get_kegg_metabolite_info(node_id)
    elif data_type == 'reactions':
        res = get_reactome_reaction_info(node_id)
    elif data_type == 'pathways':
        res = get_reactome_pathway_info(node_id)
    return res


def get_ensembl_gene_info(ensembl_id):
    display_name = ''
    try:
        metadata = get_single_ensembl_metadata_online(ensembl_id)
        display_name = metadata['display_name'] if metadata is not None and 'display_name' in metadata else ''

        infos = []
        if metadata is not None:
            # selected = ['description', 'species', 'biotype', 'db_type', 'logic_name', 'strand', 'start', 'end']
            selected = ['description', 'species']
            for key in selected:
                if key in metadata:
                    value = metadata[key]
                    if key == 'description':
                        value = value[0:value.index('[')]  # remove e.g. '[xxx]' from 'abhydrolase [xxx]'
                    infos.append({'key': key.title(), 'value': value})
    except TypeError:  # https://github.com/joewandy/GraphOmics/issues/9
        pass
    except KeyError:  # https://github.com/joewandy/GraphOmics/issues/10
        pass

    data = Harmonizome.get(Entity.GENE, name=display_name)
    try:
        description = data['description']
        infos.append({'key': 'Description', 'value': description})
    except KeyError:
        pass

    links = [
        {
            'text': 'Link to Ensembl',
            'href': 'https://www.ensembl.org/id/' + ensembl_id
        },
        {
            'text': 'Link to GeneCard',
            'href': 'https://www.genecards.org/cgi-bin/carddisp.pl?gene=' + display_name
        }
    ]
    for x in metadata['Transcript']:
        text = 'Transcript: ' + x['display_name']
        href = 'https://www.ensembl.org/id/' + x['id']
        links.append({'text': text, 'href': href})

    images = []
    data = {
        'infos': infos,
        'links': links,
        'images': images
    }
    return data


def get_uniprot_protein_info(uniprot_id):
    infos = []
    try:
        metadata = get_single_uniprot_metadata_online(uniprot_id)
        try:
            fullname = [x.text for x in metadata.soup.select('protein > recommendedname > fullname')][0]
        except IndexError:
            fullname = uniprot_id

        shortName = None
        try:
            shortname = [x.text for x in metadata.soup.select('protein > recommendedname > shortname')][0]
        except IndexError:
            pass

        if shortName is not None:
            infos.append({'key': 'Name', 'value': "{} ({})".format(fullname, shortname)})
        else:
            infos.append({'key': 'Name', 'value': "{}".format(fullname)})

        try:
            ecnumber = [x.text for x in metadata.soup.select('protein > recommendedname > ecnumber')][0]
            infos.append({'key': 'EC Number', 'value': 'EC' + ecnumber})
        except IndexError:
            pass

        # get comments
        selected = ['function', 'catalytic activity', 'enzyme regulation', 'subunit', 'pathway', 'miscellaneous',
                    'domain']
        for child in metadata.soup.find_all('comment'):
            try:
                if child['type'] in selected:
                    key = child['type'].title()
                    if key == 'Function':  # quick-hack
                        key = 'Summary'
                    infos.append({'key': key, 'value': truncate(child.text)})
            except KeyError:
                continue

        # gene ontologies
        go = []
        for child in metadata.soup.find_all('dbreference'):
            try:
                if child['type'] == 'GO':
                    props = child.find_all('property')
                    for prop in props:
                        if prop['type'] == 'term':
                            go.append(prop['value'].split(':')[1])
            except KeyError:
                continue
        go_str = '; '.join(sorted(go))
        go_str = truncate(go_str)
        infos.append({'key': 'Gene_ontologies', 'value': go_str})

    except TypeError:
        pass

    images = []
    with urllib.request.urlopen('https://swissmodel.expasy.org/repository/uniprot/' + uniprot_id + '.json') as url:
        data = json.loads(url.read().decode())
        for struct in data['result']['structures']:
            pdb_link = struct['coordinates']
            images.append(pdb_link)

    links = [
        {
            'text': 'Link to UniProt',
            'href': 'http://www.uniprot.org/uniprot/' + uniprot_id
        },
        {
            'text': 'Link to SWISS-MODEL',
            'href': 'https://swissmodel.expasy.org/repository/uniprot/' + uniprot_id
        }
    ]

    data = {
        'infos': infos,
        'links': links,
        'images': images
    }
    return data


def get_kegg_metabolite_info(compound_id):
    if '_' in compound_id:
        tokens = compound_id.split('_')
        assert len(tokens) == 2
        compound_id = tokens[0]
        peak_id = tokens[1]
    else:
        peak_id = None

    metadata = get_single_compound_metadata_online(compound_id)

    if compound_id.upper().startswith('C'):  # assume it's kegg

        # TODO: no longer used??!!

        infos = []
        selected = ['FORMULA']
        if metadata is not None:
            for key in selected:
                value = metadata[key]
                infos.append({'key': key, 'value': str(value)})

        key = 'PiMP Peak ID'
        value = peak_id
        infos.append({'key': key, 'value': str(value)})

        images = ['http://www.kegg.jp/Fig/compound/' + compound_id + '.gif']
        links = [
            {
                'text': 'Link to KEGG COMPOUND database',
                'href': 'http://www.genome.jp/dbget-bin/www_bget?cpd:' + compound_id
            }
        ]

    else:  # assume it's ChEBI

        def get_attribute(metadata, attrname):
            try:
                attr_val = getattr(metadata, attrname)
            except AttributeError:
                attr_val = ''
            return attr_val

        images = [
            'http://www.ebi.ac.uk/chebi/displayImage.do?defaultImage=true&imageIndex=0&chebiId=' + compound_id]
        links = [
            {
                'text': 'Link to ChEBI database',
                'href': 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:' + compound_id
            }
        ]

        infos = []
        key = 'PiMP Peak ID'
        value = peak_id
        infos.append({'key': key, 'value': str(value)})

        try:
            for db_link in metadata.DatabaseLinks:
                if 'KEGG COMPOUND' in str(db_link.type):
                    kegg_id = db_link.data
                    infos.append({'key': 'KEGG ID', 'value': kegg_id})
                    links.append({
                        'text': 'Link to KEGG COMPOUND database',
                        'href': 'http://www.genome.jp/dbget-bin/www_bget?cpd:' + kegg_id
                    })
        except AttributeError:
            pass

        try:
            infos.append({'key': 'FORMULA', 'value': metadata.Formulae[0].data})
        except AttributeError:
            pass
        selected = [
            ('ChEBI ID', 'chebiId'),
            ('Definition', 'definition'),
            ('Monoisotopic Mass', 'monoisotopicMass'),
            ('SMILES', 'smiles'),
            ('Inchi', 'inchi'),
            ('InchiKey', 'inchikey'),
        ]
        for key, attribute in selected:
            value = get_attribute(metadata, attribute)
            if value is not None:
                infos.append({'key': key, 'value': value})

    data = {
        'infos': infos,
        'images': images,
        'links': links
    }
    return data


def get_reactome_reaction_info(reactome_id):
    infos = []

    # res = get_reactome_content_service(reactome_id)
    # summary_str = res['summation'][0]['text']
    summary_str, is_inferred, original_species, inferred_species = _get_summary_string(reactome_id)
    infos.append({'key': 'Summary', 'value': summary_str})

    infos.append({'key': 'Species', 'value': original_species})
    if is_inferred:
        inferred = 'Inferred from %s' % inferred_species
        infos.append({'key': 'Inferred', 'value': inferred})

    # get all the participants
    temp = collections.defaultdict(list)
    results = get_reaction_entities([reactome_id])[reactome_id]
    for res in results:
        entity_id = res[1]
        display_name = res[2]
        relationship_types = res[3]
        if len(relationship_types) == 1:  # ignore the sub-complexes
            rel = relationship_types[0]
            temp[rel].append((display_name, entity_id,))

    for k, v in temp.items():
        name_list, url_list = zip(
            *v)  # https://stackoverflow.com/questions/12974474/how-to-unzip-a-list-of-tuples-into-individual-lists
        url_list = map(lambda x: 'https://reactome.org/content/detail/' + x if x is not None else '', url_list)

        name_str = ';'.join(name_list)
        url_str = ';'.join(url_list)
        infos.append({'key': k.title(), 'value': name_str, 'url': url_str})

    images = [
        'https://reactome.org/ContentService/exporter/diagram/' + reactome_id + '.jpg?sel=' + reactome_id + "&quality=7"
    ]
    links = [
        {
            'text': 'Link to Reactome database',
            'href': 'https://reactome.org/content/detail/' + reactome_id
        }
    ]
    data = {
        'infos': infos,
        'images': images,
        'links': links,
    }
    return data


def get_reactome_pathway_info(pathway_id):
    mapping, id_to_names = pathway_to_reactions([pathway_id])
    pathway_reactions = mapping[pathway_id]

    reaction_list = map(lambda x: id_to_names[x], pathway_reactions)
    reaction_str = ';'.join(sorted(reaction_list))

    url_list = map(lambda x: 'https://reactome.org/content/detail/' + x, pathway_reactions)
    url_str = ';'.join(sorted(url_list))

    # res = get_reactome_content_service(pathway_id)
    # summary_str = res['summation'][0]['text']
    summary_str, is_inferred, original_species, inferred_species = _get_summary_string(pathway_id)

    infos = [{'key': 'Summary', 'value': summary_str}]
    infos.append({'key': 'Species', 'value': original_species})
    if is_inferred:
        inferred = 'Inferred from %s' % inferred_species
        infos.append({'key': 'Inferred', 'value': inferred})
    infos.append({'key': 'Reactions', 'value': reaction_str, 'url': url_str})

    images = [
        'https://reactome.org/ContentService/exporter/diagram/' + pathway_id + '.jpg?sel=' + pathway_id
    ]
    links = [
        {
            'text': 'Link to Reactome database',
            'href': 'https://reactome.org/content/detail/' + pathway_id
        },
        {
            'text': 'SBML Export',
            'href': 'https://reactome.org/ContentService/exporter/sbml/' + pathway_id + '.xml'
        }
    ]
    data = {
        'infos': infos,
        'images': images,
        'links': links,
    }
    return data


def _get_summary_string(reactome_id):
    desc, is_inferred = get_reactome_description(reactome_id, from_parent=False)
    original_species = desc[0]['species']
    if is_inferred:
        desc, _ = get_reactome_description(reactome_id, from_parent=True)
        inferred_species = desc[0]['species']
    else:
        inferred_species = original_species

    summary_list = []
    for s in desc:
        if s['summary_text'] is None:
            continue
        st = s['summary_text'].replace(';', ',')
        summary_list.append(truncate(st))
    summary_str = ';'.join(summary_list)
    return summary_str, is_inferred, original_species, inferred_species
