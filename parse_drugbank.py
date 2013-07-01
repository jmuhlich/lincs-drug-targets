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
drugbank_drugs = sa.Table(
    'drugbank_drugs', metadata,
    sa.Column('drugbank_id', sa.String(), primary_key=True),
    sa.Column('name', sa.String()),
    sa.Column('synonyms', sa.String()),
    sa.Column('kegg_id', sa.String()),
    sa.Column('cid', sa.String()),
    sa.Column('mf', sa.String()),
    )
metadata.create_all()

datafile_name = 'drugbank.xml'
datafile = open(datafile_name)

drug_targets = {}
pubchem_cids = []
partner_to_uniprot = {}

ns = 'http://drugbank.ca'

qnames = dict((tag, lxml.etree.QName(ns, tag).text)
              for tag in ('drug', 'drug-interaction', 'partner'))

for event, element in lxml.etree.iterparse(datafile, tag=qnames['drug']):
    # We need to skip 'drug' elements in drug-interaction sub-elements. It's
    # unfortunate they re-used this tag name.
    if element.getparent().tag == qnames['drug-interaction']:
        continue
    drugbank_id = xpath(element, 'd:drugbank-id/text()')
    name = xpath(element, 'd:name/text()')
    synonyms = xpath(element, 'd:synonyms/d:synonym/text()', single=False)
    molecular_formula = xpath(element, './/d:property'
                              '[d:kind="Molecular Formula"]/d:value/text()')
    kegg_drug_id = xpath(element, './/d:external-identifier'
                         '[d:resource="KEGG Drug"]/d:identifier/text()')
    pubchem_cid = xpath(element, './/d:external-identifier'
                        '[d:resource="PubChem Compound"]/d:identifier/text()')
    partner_ids = xpath(element, 'd:targets/d:target/@partner', single=False)
    drug_targets[drugbank_id] = partner_ids
    engine.execute(drugbank_drugs.insert().values(
            drugbank_id=drugbank_id, name=name, synonyms=';'.join(synonyms),
            kegg_id=kegg_drug_id, cid=pubchem_cid, mf=molecular_formula))
    element.clear()

# Turns out it's much faster to do a second iterparse loop with a different
# tag argument than to do just one iterparse loop with a conditional on the
# tag name. The lxml internals are much more efficient at filtering tags
# than we are, and the disk I/O and buffer cache impact are negligible. It
# would be nice if the tag argument could accept a list of tag names...
datafile.seek(0)
for event, element in lxml.etree.iterparse(datafile, tag=qnames['partner']):
    partner_id = element.get('id')
    uniprot_id = xpath(element, './/d:external-identifier'
                       '[d:resource="UniProtKB"]/d:identifier/text()')
    partner_to_uniprot[partner_id] = uniprot_id
    element.clear()

for drugbank_id, partner_ids in drug_targets.items():
    partner_ids[:] = [partner_to_uniprot[i] for i in partner_ids]
    partner_ids[:] = [i for i in partner_ids if i is not None]

# Targets -- kegg_drug_id: [uniprot_ids]
#print '\n'.join(map(str, sorted(drug_targets.items())))


def to_utf8(v):
    return encodings.utf_8.encode(v)[0] if v else ''
#drug_info_good_mf = [x for x in drug_info if x.mf and ' ' not in x.mf]
# Drugs -- name, synonyms, cid, mf
#print '\n'.join('\t'.join(map(to_utf8, x.values())) for x in drug_info_good_mf)


sm_names_filename = os.path.join(os.path.dirname(sys.argv[0]),
                                 'small_molecule_all_names.130625T133836.tsv')
sm_names_file = open(sm_names_filename, 'rb')
sm_names_reader = csv.DictReader(sm_names_file, fieldnames=('hmsl_id', 'name'),
                                 dialect='excel-tab')
all_names = []
name_to_hmsl = {}
cid_to_hmsl = {}
found_conflicts = False
for row in sm_names_reader:
    hmsl_id = row['hmsl_id']
    hmsl_id = hmsl_id[:5]
    name = row['name']
    all_names.append(name)
    if name in name_to_hmsl and name_to_hmsl[name] != hmsl_id:
        print "CONFLICT: %s - %s / %s" % (name, name_to_hmsl[name], hmsl_id)
        found_conflicts = True
    name_to_hmsl[name] = hmsl_id
if found_conflicts:
    exit()

sm_filename = os.path.join(os.path.dirname(sys.argv[0]),
                           'small_molecule.130625T133836.tsv')
sm_file = open(sm_filename, 'rb')
sm_reader = csv.DictReader(sm_file, dialect='excel-tab')
cid_to_hmsl = {}
mf_to_hmsl = {}
sm_data = {}
for row in sm_reader:
    hmsl_id = row['Small Mol HMS LINCS ID']
    sm_data[hmsl_id] = row
    cid_to_hmsl[row['PubChem CID']] = hmsl_id
    mf_to_hmsl[row['Molecular Formula']] = hmsl_id

drugbank_to_hmsl = {}
for di in engine.execute(drugbank_drugs.select()):
    cid_match_id = cid_to_hmsl.get(di.cid)
    if cid_match_id:
        drugbank_to_hmsl[di.drugbank_id] = \
            [(cid_match_id, 'CID match: %s (%s)' %
              (di.cid, di.name))]
    else:
        name_matches = difflib.get_close_matches(di.name, all_names,
                                                 cutoff=0.8)

        if name_matches:
#            for name in name_matches:
#                hmsl_record = sm_data[name_to_hmsl[n]]
            drugbank_to_hmsl[di.drugbank_id] = \
                [(name_to_hmsl[n], 'Name match: %s -> %s' % (di.name, n))
                 for n in name_matches]

print  '\n'.join(sorted(map(str, drugbank_to_hmsl.items())))
