import lxml.etree

ns = 'http://drugbank.ca'

def qname(tag):
    return lxml.etree.QName(ns, tag).text

def xpath(obj, path):
    return obj.xpath(path, namespaces={'d': ns})

datafile_name = 'drugbank.xml'
datafile = open(datafile_name)

drug_targets = {}
partner_to_uniprot = {}

for event, element in lxml.etree.iterparse(datafile):
    # Second clause skips 'drug' elements in drug-interaction sub-elements. It's
    # unfortunate they re-used this tag name.
    if element.tag == qname('drug') and \
           element.getparent().tag != qname('drug-interaction'):
        drugbank_id = xpath(element, 'd:drugbank-id/text()')[0]
        cas_numbers = xpath(element, 'd:cas-number/text()')
        assert len(cas_numbers) <= 1
        partner_ids = xpath(element, 'd:targets/d:target/@partner')
        if cas_numbers:
            drug_targets[cas_numbers[0]] = partner_ids
        element.clear()
    elif element.tag == qname('partner'):
        partner_id = element.get('id')
        uniprot_ids = xpath(element, './/d:external-identifier'
                            '[d:resource="UniProtKB"]/d:identifier/text()')
        assert len(uniprot_ids) <= 1
        partner_to_uniprot[partner_id] = uniprot_ids[0] if uniprot_ids else None
        element.clear()

for drugbank_id, partner_ids in drug_targets.items():
    partner_ids[:] = [partner_to_uniprot[i] for i in partner_ids]
    partner_ids[:] = [i for i in partner_ids if i is not None]
    
print '\n'.join(map(str, sorted(drug_targets.items())))
