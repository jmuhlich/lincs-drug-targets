import csv
import re
import fileinput as fi
import kegg_utils as ku

recs = [t for t in tuple(csv.reader(fi.input(), delimiter='\t'))
        if t[-1] != '']

def scrub(s, _gunk=re.compile(ur'(?:\s+\(.*?\))+$')):
    return _gunk.sub('', s.strip('; '))


def get_kegg_drug_names(drug_id):
    names = ku.fetch_drug(drug_id)['name']
    return sorted([scrub(n) for n in
                   ([names] if isinstance(names, basestring)
                    else names)])

def match(kn, hns):
    u = kn.upper()
    for hn in hns:
        if len(hn) > 3:
            if re.search((ur'^%s\b' % hn), u, re.I):
                return True
        else:
            if u == hn.upper():
                return True
    return False

def check(kegg_names, hmsl_names):
    checked = [k for k in kegg_names if match(k, hmsl_names)]
    return len(checked) > 0

collect = dict()

for r in recs:
    hmsl, name, kegg = r
    hmsl = hmsl[:-4]
    kegg = kegg.split('|')
    c = collect.setdefault(hmsl, [set(), set()])
    c[0].add(name)
    c[1].update(set(kegg))

tocheck = dict()
for hsml, c in collect.items():
    hmsl_names, kegg_drug_ids = c
    kegg_drug_names = [(k, get_kegg_drug_names(k)) for k in kegg_drug_ids]
    c[1] = [('' if check(kns, hmsl_names) else 'x', k, kns)
            for k, kns in kegg_drug_names]
    c[0] = ' ; '.join(hmsl_names)

for hmsl in sorted(collect.keys()):
    names, kegg_data = collect[hmsl]
    for kd in sorted(kegg_data, key=lambda i: (i[0], i[1])):
        mark, kegg_id, kegg_names = kd
        print '\t'.join([mark, hmsl, names, kegg_id,
                         ' ; '.join([n.strip(' ;') for n in kegg_names])])
