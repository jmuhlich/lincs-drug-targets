import collections
import lxml.etree

ns = 'http://drugbank.ca'

def qname(tag):
    return lxml.etree.QName(ns, tag).text

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

datafile_name = 'drugbank.xml'
datafile = open(datafile_name)

drug_targets = {}
drug_info = []
pubchem_cids = []
partner_to_uniprot = {}

for event, element in lxml.etree.iterparse(datafile):
    # Second clause skips 'drug' elements in drug-interaction sub-elements. It's
    # unfortunate they re-used this tag name.
    if element.tag == qname('drug') and \
           element.getparent().tag != qname('drug-interaction'):
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
        if kegg_drug_id:
            drug_targets[kegg_drug_id] = partner_ids
        drug_info.append(collections.OrderedDict((
                    ('name', name),
                    ('synonyms', ';'.join(synonyms)),
                    ('cid', pubchem_cid),
                    ('mf', molecular_formula),
                    )))
        element.clear()
    elif element.tag == qname('partner'):
        partner_id = element.get('id')
        uniprot_id = xpath(element, './/d:external-identifier'
                           '[d:resource="UniProtKB"]/d:identifier/text()')
        partner_to_uniprot[partner_id] = uniprot_id
        element.clear()

for drugbank_id, partner_ids in drug_targets.items():
    partner_ids[:] = [partner_to_uniprot[i] for i in partner_ids]
    partner_ids[:] = [i for i in partner_ids if i is not None]


# Targets -- kegg_drug_id: [uniprot_ids]
print '\n'.join(map(str, sorted(drug_targets.items())))


def to_utf8(v):
    return encodings.utf_8.encode(v)[0] if v else ''
drug_info_good_mf = [x for x in drug_info if x['mf'] and ' ' not in x['mf']]
# Drugs -- name, synonyms, cid, mf
#print '\n'.join('\t'.join(map(to_utf8, x.values())) for x in drug_info_good_mf)