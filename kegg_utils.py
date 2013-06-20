import re
from bioservices import KeggParser

_KP = KeggParser(verbose=False)

def fetch_raw_targets(drug_id):
    dr = _KP.parse(_KP.get(drug_id))
    ts = dr.get('target', [])
    if isinstance(ts, basestring):
        ts = [ts]
    return ts

def parse_raw_target(raw_target):
    m = re.search(ur'\[HSA:\s*(\d+)(?:\s+(\d+))*\s*\]'
                  ur'(?:.*?\[HSA:\s*(\d+)(?:\s+(\d+))*\s*\])*',
                  raw_target, re.I)

    found = [] if m is None else m.groups()

    return ['hsa:%s' % s for s in found if s is not None]


def _fetch_uniprot(kegg_id):
    return _KP.parse(_KP.get(kegg_id))['dblinks']['UniProt:'].split()


def fetch_uniprot(kegg_id=None, memo=dict(), _reset=False):
    if _reset:
        assert kegg_id is None
        return memo.clear()
    return memo.setdefault(kegg_id, _fetch_uniprot(kegg_id))
