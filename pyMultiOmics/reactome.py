from collections import defaultdict
import os

import pandas as pd
import xmltodict
from bioservices.kegg import KEGG
from loguru import logger
from neo4j import GraphDatabase, basic_auth


def get_neo4j_driver():
    NEO4J_SERVER = os.getenv('NEO4J_SERVER', 'bolt://localhost:7687')
    if 'NEO4J_SERVER' not in os.environ:
        logger.warning('Using a default neo4j server: %s' % NEO4J_SERVER)

    NEO4J_USER = os.getenv('NEO4J_USER', 'neo4j')
    NEO4J_PASSWORD = os.getenv('NEO4J_PASSWORD', 'neo4j')
    if 'NEO4J_USER' not in os.environ or 'NEO4J_PASSWORD' not in os.environ:
        logger.warning('Using a default neo4j username or password: %s' % NEO4J_USER)

    try:
        neo4j_driver = GraphDatabase.driver(NEO4J_SERVER,
                                            auth=basic_auth(NEO4J_USER, NEO4J_PASSWORD))
        logger.info('Created graph database driver for %s (%s)' % (NEO4J_SERVER, NEO4J_USER))
        return neo4j_driver
    except Exception as e:
        logger.warning('Failed to connect to graph database: %s' % str(e))
        raise e


driver = get_neo4j_driver()


def get_neo4j_session():
    session = None
    try:
        session = driver.session()
    except Exception:
        raise
    return session


def get_species_list():
    results = []
    try:
        session = get_neo4j_session()
        query = """
        MATCH (n:Species) RETURN n.displayName AS name order by name        
        """
        query_res = session.run(query)
        logger.debug(query)
        for record in query_res:
            results.append(record['name'])
    finally:
        if session is not None: session.close()
    return results


def get_species_dict():
    species_list = get_species_list()
    species_dict = {}
    for idx, s in enumerate(species_list):
        species_dict[str(idx)] = s
    return species_dict


################################################################################
### Gene-related functions                                                   ###
################################################################################


def ensembl_to_uniprot(ensembl_ids, species_list):
    id_to_names = {}
    results = defaultdict(list)
    try:
        session = get_neo4j_session()
        query = """
        MATCH
            (rg:ReferenceGeneProduct)-[:referenceGene]->
            (rs:ReferenceSequence)-[:species]->(s:Species)
        WHERE
            rs.identifier IN {ensembl_ids} AND
            rs.databaseName = 'ENSEMBL' AND            
            s.displayName IN {species}
        RETURN DISTINCT
            rs.identifier AS gene_id,
            rs.databaseName AS gene_db,
            rg.identifier AS protein_id,
            rg.databaseName AS protein_db,
            rg.url as URL
        """
        params = {
            'ensembl_ids': ensembl_ids,
            'species': species_list
        }
        query_res = session.run(query, params)
        logger.debug(query)

        for record in query_res:
            gene_id = record['gene_id']
            protein_id = record['protein_id']
            results[gene_id].append(protein_id)

    finally:
        if session is not None: session.close()
    return dict(results), id_to_names


################################################################################
### Protein-related functions                                                ###
################################################################################


def uniprot_to_ensembl(uniprot_ids, species_list):
    id_to_names = {}
    results = defaultdict(list)
    try:
        session = get_neo4j_session()
        query = """
        MATCH
            (rg:ReferenceGeneProduct)-[:referenceGene]->
            (rs:ReferenceSequence)-[:species]->(s:Species)
        WHERE
            rg.identifier IN {uniprot_ids} AND
            rs.databaseName = 'ENSEMBL' AND
            s.displayName IN {species}
        RETURN DISTINCT
            rs.identifier AS gene_id,
            rs.databaseName AS gene_db,
            rg.identifier AS protein_id,
            rg.databaseName AS protein_db,
            rg.url as URL
        """
        params = {
            'uniprot_ids': uniprot_ids,
            'species': species_list
        }
        query_res = session.run(query, params)
        logger.debug(query)

        for record in query_res:
            gene_id = record['gene_id']
            protein_id = record['protein_id']
            results[protein_id].append(gene_id)

    finally:
        if session is not None: session.close()
    return dict(results), id_to_names


def uniprot_to_reaction(uniprot_ids, species_list):
    id_to_names = {}
    results = defaultdict(list)
    try:
        session = get_neo4j_session()

        # note that using hasComponent|hasMember|hasCandidate below will
        # retrieve all the sub-complexes too
        query = """
        MATCH (rle:ReactionLikeEvent)-[:input|output|catalystActivity
              |physicalEntity|regulatedBy|regulator|hasComponent|hasMember
              |hasCandidate*]->
              (pe:PhysicalEntity)-[:referenceEntity]->
              (re:ReferenceEntity)-[:referenceDatabase]->
              (rd:ReferenceDatabase)
        WHERE
            re.identifier IN {uniprot_ids} AND
            rd.displayName = 'UniProt' AND
            rle.speciesName IN {species}
        RETURN DISTINCT
            re.identifier AS protein_id,
            re.description AS description,
            rd.displayName AS protein_db,
            rle.stId AS reaction_id,
            rle.displayName AS reaction_name
        """
        params = {
            'uniprot_ids': uniprot_ids,
            'species': species_list
        }
        query_res = session.run(query, params)
        logger.debug(query)

        for record in query_res:
            protein_id = record['protein_id']
            item = {
                'reaction_id': record['reaction_id'],
                'reaction_name': record['reaction_name']
            }
            results[protein_id].append(item)
    finally:
        if session is not None: session.close()
    return dict(results), id_to_names


################################################################################
### Compound-related functions                                               ###
################################################################################


def get_all_compound_ids():
    results = []
    try:
        session = get_neo4j_session()
        query = """
        MATCH (di:DatabaseIdentifier)
        WHERE
            di.databaseName = 'COMPOUND'
        RETURN DISTINCT
            di.displayName AS compound_id
        """
        query_res = session.run(query)
        logger.debug(query)

        for record in query_res:
            key = record['compound_id'].split(':')  # e.g. 'COMPOUND:C00025'
            compound_id = key[1]
            results.append(compound_id)
    finally:
        if session is not None: session.close()
    return results


def compound_to_reaction(compound_ids, species_list):
    id_to_names = {}
    results = defaultdict(list)
    try:
        session = get_neo4j_session()
        query = """
        MATCH (rle:ReactionLikeEvent)-[:input|output|catalystActivity
              |physicalEntity|regulatedBy|regulator|hasComponent|hasMember
              |hasCandidate*]->
              (pe:PhysicalEntity)-[:crossReference|:referenceEntity]->
              (do:DatabaseObject)
        WHERE
            do.identifier IN {compound_ids} AND
            rle.speciesName IN {species}
        RETURN DISTINCT
            do.identifier AS compound_id,
            do.displayName as display_name,
            do.databaseName AS compound_db,
            rle.stId AS reaction_id,
        	rle.displayName AS reaction_name        
        """
        params = {
            'compound_ids': compound_ids,
            'species': species_list
        }
        query_res = session.run(query, params)
        logger.debug(query)

        for record in query_res:
            compound_id = record['compound_id']
            item = {
                'reaction_id': record['reaction_id'],
                'reaction_name': record['reaction_name']
            }
            results[compound_id].append(item)
            compound_name = record['display_name']
            id_to_names[compound_id] = compound_name
    finally:
        if session is not None: session.close()
    return dict(results), id_to_names


def produce_kegg_dict(kegg_location, param):
    with open(kegg_location) as kegg_cmpd_file:
        cmpd_dict = xmltodict.parse(kegg_cmpd_file.read())

    kegg_dict = {}
    for compound in cmpd_dict['compounds']['compound']:
        kegg_dict[compound[param]] = compound['formula']

    return kegg_dict


################################################################################
### Reaction-related functions                                               ###
################################################################################


# get all the entities involved in a reaction
def get_reaction_entities(reaction_ids):
    results = defaultdict(list)
    try:
        session = get_neo4j_session()
        query = """
        MATCH (rle:ReactionLikeEvent)-[rr:input|output|catalystActivity
              |physicalEntity|regulatedBy|regulator|hasComponent|hasMember
              |hasCandidate*]->(dbo:DatabaseObject)
        WHERE
            rle.stId IN {reaction_ids}
        RETURN
            rle.stId AS reaction_id,
            dbo.stId AS entity_id,
            dbo.schemaClass AS schema_class,
            dbo.displayName as display_name,
            extract(rel IN rr | type(rel)) AS types
        """
        params = {
            'reaction_ids': reaction_ids
        }
        query_res = session.run(query, params)
        logger.debug(query)

        for record in query_res:
            reaction_id = record['reaction_id']
            entity_id = record['entity_id']
            schema_class = record['schema_class']
            display_name = record['display_name']
            relationship_types = record['types']
            item = (schema_class, entity_id, display_name, relationship_types)
            results[reaction_id].append(item)
    finally:
        if session is not None: session.close()
    return results


def reaction_to_uniprot(reaction_ids, species_list):
    id_to_names = {}
    results = defaultdict(list)
    try:
        session = get_neo4j_session()

        # note that using hasComponent|hasMember|hasCandidate below will
        # retrieve all the sub-complexes too
        query = """
        MATCH (rle:ReactionLikeEvent)-[:input|output|catalystActivity
              |physicalEntity|regulatedBy|regulator|hasComponent|hasMember
              |hasCandidate*]->
              (pe:PhysicalEntity)-[:referenceEntity]->
              (re:ReferenceEntity)-[:referenceDatabase]->
              (rd:ReferenceDatabase)
        WHERE
            rle.stId IN {reaction_ids} AND
            rd.displayName = 'UniProt' AND
            rle.speciesName IN {species}
        RETURN DISTINCT
            re.identifier AS protein_id,
            re.description AS description,
            rd.displayName AS protein_db,
            rle.stId AS reaction_id,
            rle.displayName AS reaction_name
        """
        params = {
            'reaction_ids': reaction_ids,
            'species': species_list
        }
        query_res = session.run(query, params)
        logger.debug(query)

        for record in query_res:
            protein_id = record['protein_id']
            reaction_id = record['reaction_id']
            results[reaction_id].append(protein_id)
    finally:
        if session is not None: session.close()
    return dict(results), id_to_names


def reaction_to_compound(reaction_ids, species_list, use_kegg=False):
    id_to_names = {}
    results = defaultdict(list)
    try:
        session = get_neo4j_session()
        query = """
        MATCH (rle:ReactionLikeEvent)-[:input|output|catalystActivity
              |physicalEntity|regulatedBy|regulator|hasComponent|hasMember
              |hasCandidate*]->
              (pe:PhysicalEntity)-[:crossReference|:referenceEntity]->
              (do:DatabaseObject)
        WHERE
            rle.stId IN {reaction_ids} AND        
            (do.databaseName = 'COMPOUND' OR do.databaseName = 'ChEBI') AND
            rle.speciesName IN {species}
        RETURN DISTINCT
            do.identifier as compound_id,
            do.displayName as display_name,            
            do.databaseName AS compound_db,
            rle.stId AS reaction_id,
        	rle.displayName AS reaction_name
        """
        params = {
            'reaction_ids': reaction_ids,
            'species': species_list
        }
        query_res = session.run(query, params)
        logger.debug(query)

        for record in query_res:
            reaction_id = record['reaction_id']
            compound_id = record['compound_id']
            display_name = record['display_name']
            database_name = record['compound_db']
            # TODO: find better ways to remove duplicates between KEGG and ChEBI?
            valid = False
            if use_kegg:
                if database_name == 'COMPOUND':
                    valid = True
            else:
                if database_name == 'ChEBI':
                    valid = True
            if valid:
                results[reaction_id].append(compound_id)
                id_to_names[compound_id] = display_name
    finally:
        if session is not None: session.close()
    return dict(results), id_to_names


def rchop(thestring, ending):
    if thestring.endswith(ending):
        return thestring[:-len(ending)]
    return thestring


def reaction_to_pathway(reaction_ids, species_list, metabolic_pathway_only, leaf=True):
    id_to_names = {}
    results = defaultdict(list)
    try:
        session = get_neo4j_session()

        # initial match clause in the query
        query = """
        MATCH (tp:TopLevelPathway)-[:hasEvent*]->
              (p:Pathway)-[:hasEvent*]->(rle:ReactionLikeEvent)
        WHERE
            tp.speciesName IN {species} AND        
            rle.stId IN {reaction_ids} AND            
        """

        if leaf:  # retrieve only the leaf nodes in the pathway hierarchy
            query += " (p)-[:hasEvent]->(rle) AND "

        if metabolic_pathway_only:  # only retrieves metabolic pathways
            query += " tp.displayName = 'Metabolism' AND "

        # remove last AND
        query = rchop(query.strip(), 'AND')

        # add return clause
        query += """
        RETURN
            rle.stId AS reaction_id,
            rle.displayName AS reaction_name,
            rle.speciesName AS reaction_species,
            p.stId AS pathway_id,
            p.displayName AS pathway_name,
            tp.speciesName AS pathway_species
        """

        params = {
            'reaction_ids': reaction_ids,
            'species': species_list
        }
        query_res = session.run(query, params)
        logger.debug(query)

        for record in query_res:
            reaction_id = record['reaction_id']
            reaction_name = record['reaction_name']
            reaction_species = record['reaction_species']
            pathway_id = record['pathway_id']
            pathway_name = record['pathway_name']
            pathway_species = record['pathway_species']
            item = {
                'pathway_id': pathway_id,
                'pathway_name': pathway_name
            }
            results[reaction_id].append(item)
            id_to_names[reaction_id] = {'name': reaction_name, 'species': reaction_species}
            id_to_names[pathway_id] = {'name': pathway_name, 'species': pathway_species}
    finally:
        if session is not None: session.close()
    return dict(results), id_to_names


################################################################################
### Pathway-related functions                                                ###
################################################################################


def pathway_to_reactions(pathway_ids):
    id_to_names = {}
    results = defaultdict(list)
    try:
        session = get_neo4j_session()
        # retrieve only the leaf nodes in the pathway hierarchy
        query = """
        MATCH (p:Pathway)-[:hasEvent*]->(rle:ReactionLikeEvent)
        WHERE
            p.stId IN {pathway_ids} AND
            (p)-[:hasEvent]->(rle)
        RETURN
            rle.stId AS reaction_id,
            rle.displayName AS reaction_name,
            p.stId AS pathway_id,
            p.displayName AS pathway_name
        """
        params = {
            'pathway_ids': pathway_ids
        }
        query_res = session.run(query, params)
        logger.debug(query)

        for record in query_res:
            reaction_id = record['reaction_id']
            reaction_name = record['reaction_name']
            pathway_id = record['pathway_id']
            pathway_name = record['pathway_name']
            results[pathway_id].append(reaction_id)
            id_to_names[reaction_id] = reaction_name
            id_to_names[pathway_id] = pathway_name
    finally:
        if session is not None: session.close()
    return dict(results), id_to_names


def get_reactome_description(reactome_id, from_parent=False):
    results = []
    try:
        session = get_neo4j_session()
        if from_parent:
            query = """
            MATCH (dbo1:DatabaseObject)<-[:inferredTo*]-(dbo2:DatabaseObject)-[:summation|:literatureReference]-(ss)
                    WHERE
                        dbo1.stId = {reactome_id} AND
                        dbo2.isInferred = False
                    RETURN
                        dbo2.stId as reactome_id,
                        dbo2.speciesName as species,
                        dbo2.isInferred as inferred,
                        ss.displayName as display_name,
                        ss.text as summary_text,
                        ss.schemaClass as summary_type,
                        properties(ss) as summary
            """
        else:
            query = """
            MATCH (dbo:DatabaseObject)-[:summation|:literatureReference]-(ss)
                    WHERE
                        dbo.stId = {reactome_id}
                    RETURN
                        dbo.stId as reactome_id,
                        dbo.speciesName as species,
                        dbo.isInferred as inferred,
                        ss.displayName as display_name,
                        ss.text as summary_text,
                        ss.schemaClass as summary_type,
                        properties(ss) as summary_props
            """
        params = {
            'reactome_id': reactome_id,
        }
        query_res = session.run(query, params)
        logger.debug(query)
        results = query_res.data()
        first_data = results[0]
        is_inferred = first_data['inferred']
    finally:
        if session is not None: session.close()
    return results, is_inferred


def retrieve_kegg_formula(reactome_compound_name):
    k = KEGG()
    compound_name = reactome_compound_name.replace('COMPOUND', 'cpd')
    res = k.get(compound_name).split('\n')
    for line in res:
        if line.startswith('FORMULA'):
            formula = line.split()[1]  # get the second token
            return formula
    return None


def get_all_pathways(species_list):
    results = []
    try:
        session = get_neo4j_session()

        # retrieve only the leaf nodes in the pathway hierarchy
        query = """
            MATCH (tp:TopLevelPathway)-[:hasEvent*]->(p:Pathway)-[:hasEvent*]->(rle:ReactionLikeEvent)
            WHERE
                tp.displayName = 'Metabolism' AND
                tp.speciesName IN {species_list} AND
                (p)-[:hasEvent]->(rle)
            RETURN DISTINCT
                p.speciesName AS species_name,            
                p.displayName AS pathway_name,
                p.stId AS pathway_id                       
            ORDER BY species_name, pathway_name
        """
        params = {
            'species_list': species_list
        }
        query_res = session.run(query, params)
        logger.debug(query)

        for record in query_res:
            pathway_species = record['species_name']
            pathway_name = record['pathway_name']
            pathway_id = record['pathway_id']
            results.append((pathway_species, pathway_name, pathway_id))
    finally:
        if session is not None: session.close()
    return results


def get_all_pathways_formulae(species):
    results = defaultdict(set)
    pathway_id_to_name = {}
    try:
        session = get_neo4j_session()

        # retrieve only the leaf nodes in the pathway hierarchy
        query = """
        MATCH (tp:TopLevelPathway)-[:hasEvent*]->
              (p:Pathway)-[:hasEvent*]->(rle:ReactionLikeEvent),
              (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent
              |hasMember|hasCandidate*]->(pe:PhysicalEntity),
              (pe:PhysicalEntity)-[:crossReference]->(di:DatabaseIdentifier)<-[:crossReference]-(rm:ReferenceMolecule)
        WHERE
              tp.displayName = 'Metabolism' AND
              tp.speciesName = {species} AND
              di.databaseName = 'COMPOUND' AND
              (p)-[:hasEvent]->(rle)
        RETURN DISTINCT
            p.schemaClass,
            p.displayName AS pathway_name,
            p.stId AS pathway_id,
            di.displayName as compound_name,
            rm.formula AS formula,
            di.url
        """
        params = {
            'species': species
        }
        query_res = session.run(query, params)
        logger.debug(query)

        i = 0
        retrieved = {}
        for record in query_res:
            pathway_id = record['pathway_id']
            pathway_name = record['pathway_name']
            pathway_id_to_name[pathway_id] = pathway_name
            compound_name = record['compound_name']
            formula = record['formula']
            if formula is None:
                if compound_name not in retrieved:
                    formula = retrieve_kegg_formula(compound_name)
                    logger.debug('Missing formula for %s, retrieved %s from kegg' %
                                 (compound_name, formula))
                    retrieved[compound_name] = formula
                else:
                    formula = retrieved[compound_name]
            assert formula is not None, 'Formula is missing for %s' % compound_name
            results[pathway_id].add(formula)
    finally:
        if session is not None: session.close()
    return dict(results), pathway_id_to_name


################################################################################
### Analysis functions                                                       ###
################################################################################


def get_reaction_ids(mapping):
    all_reactions = []
    for key in mapping:
        rs = mapping[key]
        rids = [r['reaction_id'] for r in rs]
        all_reactions.extend(rids)
    return all_reactions


def get_reactions_from_mapping(mapping):
    reaction_names = {}
    reaction_members = defaultdict(list)
    for key in mapping:
        for reaction in mapping[key]:
            r_id = reaction['reaction_id']
            r_name = reaction['reaction_name']
            reaction_names[r_id] = r_name
            reaction_members[r_id].append(key)
    assert reaction_names.keys() == reaction_members.keys()
    return reaction_names, dict(reaction_members)


def get_protein_to_gene(mapping):
    protein_to_gene = defaultdict(list)
    for gene_id in mapping:
        for protein_id in mapping[gene_id]:
            protein_to_gene[protein_id].append(gene_id)
    return dict(protein_to_gene)


def merge_two_dicts(x, y):
    # https://stackoverflow.com/questions/38987/how-to-merge-two-dictionaries-in-a-single-expression
    z = x.copy()  # start with x's keys and values
    z.update(y)  # modifies z with y's keys and values & returns None
    return z


def get_coverage(observed_count, total_count):
    try:
        return observed_count / float(total_count)
    except ZeroDivisionError:
        return 0


def get_reaction_df(transcript_mapping, protein_mapping, compound_mapping,
                    pathway_mapping, species_list):
    r_name_1, r_members_1 = get_reactions_from_mapping(protein_mapping)
    r_name_2, r_members_2 = get_reactions_from_mapping(compound_mapping)
    reaction_names = merge_two_dicts(r_name_1, r_name_2)
    protein_to_gene = get_protein_to_gene(transcript_mapping)

    pathway_compound_counts = defaultdict(int)
    pathway_protein_counts = defaultdict(int)

    reaction_ids = set(list(r_members_1.keys()) + list(r_members_2.keys()))
    reaction_entities = get_reaction_entities(list(reaction_ids))
    rows = []
    for reaction_id in reaction_ids:

        observed_protein_count = 0
        observed_compound_count = 0
        protein_str = ''
        compound_str = ''

        if reaction_id in r_members_1:
            proteins = r_members_1[reaction_id]
            observed_protein_count = len(proteins)
            for prot in proteins:
                if prot in protein_to_gene:
                    protein_str += '%s (%s):' % (prot, ':'.join(
                        protein_to_gene[prot]))
                else:
                    protein_str += '%s:' % prot
            protein_str = protein_str.rstrip(':')  # remove last :

        if reaction_id in r_members_2:
            compounds = r_members_2[reaction_id]
            observed_compound_count = len(compounds)
            compound_str = ':'.join(compounds)

        entities = reaction_entities[reaction_id]
        all_compound_count = len([x for x in entities
                                  if x[0] == 'SimpleEntity'])
        all_protein_count = len([x for x in entities
                                 if x[0] == 'EntityWithAccessionedSequence'])

        reaction_name = reaction_names[reaction_id]
        protein_coverage = get_coverage(observed_protein_count,
                                        all_protein_count)
        compound_coverage = get_coverage(observed_compound_count,
                                         all_compound_count)
        s1 = observed_protein_count + observed_compound_count
        s2 = all_protein_count + all_compound_count
        all_coverage = get_coverage(s1, s2)

        if reaction_id in pathway_mapping:
            pathway_id_str = ':'.join([x['pathway_id']
                                       for x in pathway_mapping[reaction_id]])
            pathway_name_str = ':'.join([x['pathway_name']
                                         for x in pathway_mapping[reaction_id]])
            protein_coverage_str = '%.2f' % protein_coverage
            compound_coverage_str = '%.2f' % compound_coverage
            all_coverage_str = '%.2f' % all_coverage
            row = (reaction_id,
                   reaction_name,
                   protein_coverage_str,
                   compound_coverage_str,
                   all_coverage_str,
                   protein_str,
                   observed_protein_count,
                   all_protein_count,
                   compound_str,
                   observed_compound_count,
                   all_compound_count,
                   pathway_id_str,
                   pathway_name_str)
            rows.append(row)

            # accumulate the count for pathways too
            for x in pathway_mapping[reaction_id]:
                pathway_compound_counts[x['pathway_id']] += observed_compound_count
                pathway_protein_counts[x['pathway_id']] += observed_protein_count

    df = pd.DataFrame(rows, columns=['reaction_id',
                                     'reaction_name',
                                     'protein_coverage',
                                     'compound_coverage',
                                     'all_coverage',
                                     'protein',
                                     'observed_protein_count',
                                     'all_protein_count',
                                     'compound',
                                     'observed_compound_count',
                                     'all_compound_count',
                                     'pathway_ids',
                                     'pathway_names'])
    return df, pathway_compound_counts, pathway_protein_counts
