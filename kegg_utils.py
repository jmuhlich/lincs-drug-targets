import re
from bioservices import KeggParser

_KP = KeggParser(verbose=False)


def _fetch_drug(drug_id):
    return _KP.parse(_KP.get(drug_id))


def fetch_drug(drug_id=None, memo=dict(), _reset=False):
    if _reset:
        assert drug_id is None
        return memo.clear()
    return memo.setdefault(drug_id, _fetch_drug(drug_id))


def fetch_raw_targets(drug_id):
    dr = fetch_drug(drug_id)
    ts = dr.get('target', [])
    if isinstance(ts, basestring):
        ts = [ts]
    return ts


def parse_raw_target(raw_target,
                     _id_re=re.compile(ur'\[HSA:\s*(\d+)(?:\s+(\d+))*\s*\]'
                                       ur'(?:.*?\[HSA:\s*(\d+)'
                                       ur'(?:\s+(\d+))*\s*\])*', re.I)):
    m = _id_re.search(raw_target)

    found = [] if m is None else m.groups()

    return ['hsa:%s' % s for s in found if s is not None]



def fetch_raw_pathways(drug_id, _ws_re=re.compile(ur'\s+')):
    pws = fetch_drug(drug_id).get('pathway', {})
    if isinstance(pws, basestring):
        raw_id, name = _ws_re.split(pws, 1)
        pws = {raw_id: name}
    return pws


def parse_raw_pathway(pwid, pwname,
                      _id_re=re.compile(ur'^((?:hsa|ko)\d+)', re.I)):
    m = _id_re.search(pwid.strip())
    return (m.group(1).lower(), pwname.strip()) if m else None


def _fetch_uniprot(kegg_id):
    return _KP.parse(_KP.get(kegg_id))['dblinks']['UniProt:'].split()


def fetch_uniprot(kegg_id=None, memo=dict(), _reset=False):
    if _reset:
        assert kegg_id is None
        return memo.clear()
    return memo.setdefault(kegg_id, _fetch_uniprot(kegg_id))
