import sys
import re
import csv
import fileinput as fi
import re
import collections as co

def cleanup(s, _nl=re.compile(ur'\n+|^\s*None\s*$')):
    return _nl.sub('', s)

colnames = tuple('small_mol_hms_lincs_id sm_name alternative_names '
                 'lincs_id pubchem_cid datarecord_id hms_dataset_id salt_id '
                 'chembl_id chebi_id inchi inchi_key smiles molecular_mass '
                 'molecular_formula'.split())

rec = co.namedtuple('Record', colnames)


def _check_kwargs(kwargs, _colnames=set(colnames)):
    assert set(kwargs.keys()).issubset(_colnames)


def update_record(r, **kwargs):
    _check_kwargs(kwargs)
    new = [None for _ in colnames]
    for i, c in enumerate(colnames):
        new[i] = kwargs[c] if c in kwargs else r[i]
    return rec(*new)

sep = ';'
sep_re = re.compile(ur'\s*[,;/](?: *[,;/])*\s*')
comma_sep_re = re.compile(ur'\s*,(?:\s*,)*\s*')
semicolon_sep_re = re.compile(ur'\s*;(?:\s*;)*\s*')
slash_sep_re = re.compile(ur'\s*/(?:\s*/)*\s*')
ws_re = re.compile(ur'\s*')

def isseq(s):
    return hasattr(s, '__iter__')

def nodups(s):
    assert isseq(s)
    ret = [None for _ in s]
    i = 0
    seen = set()
    for x in s:
        if x not in seen:
            seen.add(x)
            ret[i] = x
            i += 1
    return type(s)(ret[:i])

def cleanup_names(r):
    expected = {
        '10001-101': ('(R)- Roscovitine', 'CYC202;Seliciclib'),
        '10209-101': ('BMS-536924, KIN001-126', ''),
        '10212-101': ('KIN001-021/CGP082996', 'CINK4;ZINC01493715'),
        '10213-101': ('KIN001-111/A770041', ''),
        '10214-101': ('KIN001-123/WZ-1-84', ''),
        '10298-101': ('LFM-A13/DDE-28', ''),
    }
    k = r.small_mol_hms_lincs_id
    if k not in expected: return r
    n, s = expected[k]
    assert r.sm_name == n and r.alternative_names == s
    names = sep_re.split(n)
    name = names[0]
    alts = nodups(names[1:] +
                  (sep_re.split(s) if len(s) else []))

    if k == '10001-101':
        name = ws_re.sub('', name)
        name, alts[-1] = alts[-1], name

    ans = sep.join(alts)
    return update_record(r, sm_name=name, alternative_names=ans)


def cleanup_alts(r):
    alts = nodups(sep_re.split(r.alternative_names))
    new = sep.join(alts)
    if new == r.alternative_names:
        return r
    expected = {
        '10035-101': 'KIN001-102;;AKT inhibitor VIII;Akt1/2 kinase inhibitor',
        '10035-109': 'KIN001-102;;AKT inhibitor VIII;Akt1/2 kinase inhibitor',
        '10113-101': 'Kin001-237;c-Met/Ron dual kinase inhibitor',
        '10215-101': 'GSK319347A;HMS3229F19;CAY10576;ZINC19795634;AKOS002364333;NCGC00242056-01;benzimidazole-thiophene carbonitrile;12e, EC-000.2375',
        '10282-101': 'suberoylanilide hydroxamic acid (SAHA), Zolinza',
        '10300-101': '5-azacytidine, Vidaza, AZA',
        '10301-101': '5-aza-2\'-deoxycytidine, Dacogen, DAC',
        '10304-105': 'AG014699, PF-01367338',
        '10311-101': 'Entinostat, MS-27-275',
        '10316-101': 'SEN196, selisistat',
    }
    k = r.small_mol_hms_lincs_id
    assert k in expected

    if k == '10113-101':
        return r

    if k.startswith('10035'):
        new = sep.join(nodups(semicolon_sep_re.split(r.alternative_names)))
    elif k == '10215-101':
        new = sep.join([s for s in alts if s != '12e'])

    return update_record(r, alternative_names=new)


def cleanup_record(r):
    r = rec(*(cleanup(s) for s in r))
    r = cleanup_names(r)
    r = cleanup_alts(r)
    return rec('HMSL%s' % r[0], *r[1:])

reader = csv.reader(fi.input())
writer = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')

writer.writerow(next(reader))
writer.writerows([cleanup_record(rec(*r)) for r in reader])
