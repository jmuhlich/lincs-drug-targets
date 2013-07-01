import csv
import re
import fileinput as fi
import kegg_utils as ku

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

def collect_data(input_handle):
    recs = tuple(csv.reader(input_handle, delimiter='\t'))

    collect = dict()

    for r in recs:
        hmsl, name, kegg = r
        hmsl = hmsl[:-4]
        c = collect.setdefault(hmsl, [set(), set()])
        c[0].add(name)
        if len(kegg):
            kegg = kegg.split('|')
            c[1].update(set(kegg))

    tocheck = dict()
    for c in collect.values():
        hmsl_names, kegg_drug_ids = c
        kegg_drug_names = [(k, get_kegg_drug_names(k)) for k in kegg_drug_ids]
        c[1] = [('' if check(kns, hmsl_names) else 'x', k, kns)
                for k, kns in kegg_drug_names]
        c[0] = ' ; '.join(hmsl_names)

    return collect


def print_output(collect):
    for hmsl in sorted(collect.keys()):
        names, kegg_data = collect[hmsl]
        recs = sorted(kegg_data, key=lambda i: (i[0], i[1]))

        if not recs:
            recs = [('', '', [''])]

        for kd in recs:
            mark, kegg_id, kegg_names = kd
            print '\t'.join([mark, hmsl, names, kegg_id,
                             ' ; '.join([n.strip(' ;') for n in kegg_names])])

if __name__ == '__main__':
    print_output(collect_data(fi.input()))
