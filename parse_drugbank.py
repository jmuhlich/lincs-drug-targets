import os
import sys
import lxml.etree
import csv
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

db_file = 'drugbank.sqlite'
#db_file = ':memory:'
engine = sa.create_engine('sqlite:///' + db_file)
conn = engine.connect()
metadata = sa.MetaData(bind=conn)
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

with conn.begin() as trans:
    for event, element in lxml.etree.iterparse(datafile, tag=qnames['drug']):
        # We need to skip 'drug' elements in drug-interaction sub-elements. It's
        # unfortunate they re-used this tag name.
        if element.getparent().tag == qnames['drug-interaction']:
            continue
        drug_id = xpath(element, 'd:drugbank-id/text()')
        name = xpath(element, 'd:name/text()')
        synonyms = xpath(element, 'd:synonyms/d:synonym/text()', single=False)
        synonyms += xpath(element, 'd:brands/d:brand/text()', single=False)
        molecular_formula = xpath(element, './/d:property'
                                  '[d:kind="Molecular Formula"]/d:value/text()')
        kegg_id = xpath(element, './/d:external-identifier'
                        '[d:resource="KEGG Drug"]/d:identifier/text()')
        pubchem_cid = xpath(element, './/d:external-identifier'
                            '[d:resource="PubChem Compound"]/d:identifier/text()')
        partner_ids = xpath(element, 'd:targets/d:target/@partner', single=False)
        conn.execute(drugbank_drug.insert().
                     values(drug_id=drug_id, name=name, synonyms=synonyms,
                            kegg_id=kegg_id, pubchem_cid=pubchem_cid,
                            molecular_formula=molecular_formula,
                            partners=partner_ids))
        conn.execute(drugbank_name.insert().
                     values(drug_id=drug_id, name=name.lower()))
        for s in synonyms:
            conn.execute(drugbank_name.insert().
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

with conn.begin() as trans:
    for rec in conn.execute(drugbank_drug.select()):
        new_values = dict(rec)
        new_values['partners'] = map(partner_to_uniprot.__getitem__, rec.partners)
        new_values['partners'] = filter(None, new_values['partners'])
        conn.execute(drugbank_drug.update().
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
    rec[0] for rec in conn.execute(sa.select([drugbank_name.c.name]))]

sm_filename = os.path.join(os.path.dirname(sys.argv[0]),
                           'small_molecule.130624M134120.tsv')
sm_file = open(sm_filename, 'rb')
sm_reader = csv.reader(sm_file, dialect='excel-tab')
sm_fields = [f.lower().replace(' ', '_') for f in sm_reader.next()]
sm_fields[0] = 'sm_id'
hmsl_sm = sa.Table(
    'hmsl_sm', metadata,
    *[sa.Column(f, sa.String()) for f in sm_fields]
)
hmsl_sm.append_constraint(sa.PrimaryKeyConstraint(hmsl_sm.c.sm_id))
hmsl_sm.c.alternative_names.type = sa.PickleType()
metadata.create_all(tables=[hmsl_sm])
with conn.begin() as trans:
    for row in sm_reader:
        row[0] = row[0][:-4]
        row[2] = row[2].split(';')
        try:
            conn.execute(hmsl_sm.insert().values(row))
        except sa.exc.IntegrityError as e:
            rec = conn.execute(hmsl_sm.select().
                               where(hmsl_sm.c.sm_id == row[0])).first()
            if rec:
                new_rec = dict(rec)
                new_rec['alternative_names'] = list(set(
                        rec.alternative_names +
                        [row[sm_fields.index('sm_name')]] +
                        row[sm_fields.index('alternative_names')]))
                if not rec.pubchem_cid:
                    new_rec['pubchem_cid'] = row[sm_fields.index('pubchem_cid')]
                conn.execute(hmsl_sm.update().
                             where(hmsl_sm.c.sm_id == new_rec['sm_id']).
                             values(new_rec))

hmsl_to_drugbank = {}

for sm in conn.execute(hmsl_sm.select()):

    hmsl_names = [s.lower() for s in [sm.sm_name] + sm.alternative_names]
    for name in hmsl_names:
        match = conn.execute(sa.select([drugbank_name.c.drug_id]).
                             where(drugbank_name.c.name == name)
                             ).scalar()
        if match:
            break
    if match:
        print "%s\t%s\tName: %s" % (sm.sm_id, match, name)
        continue

    match = conn.execute(sa.select([drugbank_drug.c.drug_id]).
                         where(drugbank_drug.c.pubchem_cid == 
                               sm.pubchem_cid)
                         ).scalar()
    if match:
        print "%s\t%s\tPubChem CID: %s" % (sm.sm_id, match, sm.pubchem_cid)
        continue

    print "%s\t\tNO MATCH" % sm.sm_id
