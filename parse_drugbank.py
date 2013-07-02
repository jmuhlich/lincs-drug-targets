import os
import sys
import lxml.etree
import csv
import difflib
import sqlalchemy as sa

def xpath(obj, path, single=True):
    result = obj.xpath(path, namespaces={'d': ns})
    if single:
        if len(result) == 0:
            result = None
        elif len(result) == 1:
            result = result[0]
        else:
            raise ValueError("XPath expression matches more than one value")
    return result

engine = sa.create_engine('sqlite:///:memory:')
metadata = sa.MetaData(bind=engine)
drugbank_drug = sa.Table(
    'drugbank_drug', metadata,
    sa.Column('drug_id', sa.String(), primary_key=True),
    sa.Column('name', sa.String()),
    sa.Column('synonyms', sa.PickleType()),  # list of strings
    sa.Column('kegg_id', sa.String()),
    sa.Column('pubchem_cid', sa.String()),
    sa.Column('molecular_formula', sa.String()),
    sa.Column('partners', sa.PickleType()),  # list of strings
    )
drugbank_name = sa.Table(
    'drugbank_name', metadata,
    sa.Column('drug_id', sa.String()),
    sa.Column('name', sa.String(), index=True),
    )
metadata.create_all()

datafile_name = 'drugbank.xml'
datafile = open(datafile_name)

ns = 'http://drugbank.ca'

qnames = dict((tag, lxml.etree.QName(ns, tag).text)
              for tag in ('drug', 'drug-interaction', 'partner'))

for event, element in lxml.etree.iterparse(datafile, tag=qnames['drug']):
    # We need to skip 'drug' elements in drug-interaction sub-elements. It's
    # unfortunate they re-used this tag name.
    if element.getparent().tag == qnames['drug-interaction']:
        continue
    drug_id = xpath(element, 'd:drugbank-id/text()')
    name = xpath(element, 'd:name/text()')
    synonyms = xpath(element, 'd:synonyms/d:synonym/text()', single=False)
    molecular_formula = xpath(element, './/d:property'
                              '[d:kind="Molecular Formula"]/d:value/text()')
    kegg_id = xpath(element, './/d:external-identifier'
                    '[d:resource="KEGG Drug"]/d:identifier/text()')
    pubchem_cid = xpath(element, './/d:external-identifier'
                        '[d:resource="PubChem Compound"]/d:identifier/text()')
    partner_ids = xpath(element, 'd:targets/d:target/@partner', single=False)
    engine.execute(drugbank_drug.insert().
                   values(drug_id=drug_id, name=name, synonyms=synonyms,
                          kegg_id=kegg_id, pubchem_cid=pubchem_cid,
                          molecular_formula=molecular_formula,
                          partners=partner_ids))
    engine.execute(drugbank_name.insert().
                   values(drug_id=drug_id, name=name.lower()))
    for s in synonyms:
        engine.execute(drugbank_name.insert().
                       values(drug_id=drug_id, name=s.lower()))
    element.clear()

# Turns out it's much faster to do a second iterparse loop with a different
# tag argument than to do just one iterparse loop with a conditional on the
# tag name. The lxml internals are much more efficient at filtering tags
# than we are, and the disk I/O and buffer cache impact are negligible. It
# would be nice if the tag argument could accept a list of tag names...
datafile.seek(0)
partner_to_uniprot = {}
for event, element in lxml.etree.iterparse(datafile, tag=qnames['partner']):
    partner_id = element.get('id')
    uniprot_id = xpath(element, './/d:external-identifier'
                       '[d:resource="UniProtKB"]/d:identifier/text()')
    partner_to_uniprot[partner_id] = uniprot_id
    element.clear()

# Use a transaction to avoid re-reading our own updates (partners would already
# be updated).
with engine.begin() as connection:
    for rec in connection.execute(drugbank_drug.select()):
        new_values = dict(rec)
        new_values['partners'] = map(partner_to_uniprot.__getitem__, rec.partners)
        new_values['partners'] = filter(None, new_values['partners'])
        connection.execute(drugbank_drug.update().
                           where(drugbank_drug.c.drug_id == rec.drug_id).
                           values(**new_values))

# Targets -- kegg_drug_id: [uniprot_ids]
#print '\n'.join(map(str, sorted(drug_targets.items())))


def to_utf8(v):
    return encodings.utf_8.encode(v)[0] if v else ''
#drug_info_good_mf = [x for x in drug_info if x.mf and ' ' not in x.mf]
# Drugs -- name, synonyms, cid, mf
#print '\n'.join('\t'.join(map(to_utf8, x.values())) for x in drug_info_good_mf)

drugbank_names = [
    rec[0] for rec in engine.execute(sa.select([drugbank_name.c.name]))]

sm_filename = os.path.join(os.path.dirname(sys.argv[0]),
                           'small_molecule.130624M134120.tsv')
sm_file = open(sm_filename, 'rb')
sm_reader = csv.reader(sm_file, dialect='excel-tab')
sm_fields = [f.lower().replace(' ', '_') for f in sm_reader.next()]
hmsl_sm = sa.Table(
    'hmsl_sm', metadata,
    *[sa.Column(f, sa.String()) for f in sm_fields]
)
hmsl_sm.c.small_mol_hms_lincs_id.primary_key = True
hmsl_sm.c.alternative_names.type = sa.PickleType()
metadata.create_all(tables=[hmsl_sm])
for row in sm_reader:
    row[0] = row[0][:-4]
    row[2] = row[2].split(';')
    engine.execute(hmsl_sm.insert().values(row))

hmsl_to_drugbank = {}

for sm in engine.execute(hmsl_sm.select()):

    hmsl_names = [s.lower() for s in [sm.sm_name] + sm.alternative_names]
    for name in hmsl_names:
        match = engine.execute(sa.select([drugbank_name.c.drug_id]).
                               where(drugbank_name.c.name == name)
                               ).scalar()
        if match:
            break
    if match:
        print "%s\t%s\tName: %s" % \
            (sm.small_mol_hms_lincs_id, match, sm.sm_name)
        continue

    match = engine.execute(sa.select([drugbank_drug.c.drug_id]).
                           where(drugbank_drug.c.pubchem_cid == 
                                 sm.pubchem_cid)
                           ).scalar()
    if match:
        print "%s\t%s\tPubChem CID: %s" % \
            (sm.small_mol_hms_lincs_id, match, sm.pubchem_cid)
        continue

    print "%s\t\tNO MATCH" % sm.small_mol_hms_lincs_id

#     else:
#         name_matches = difflib.get_close_matches(di.name, all_names,
#                                                  cutoff=0.8)

#         if name_matches:
# #            for name in name_matches:
# #                hmsl_record = sm_data[name_to_hmsl[n]]
#             drugbank_to_hmsl[di.drug_id] = \
#                 [(name_to_hmsl[n], 'Name match: %s -> %s' % (di.name, n))
#                  for n in name_matches]

print  '\n'.join(sorted(map(str, drugbank_to_hmsl.items())))
