from bioservices import UniProt

_UP = UniProt(verbose=False)

def _fetch_name(uniprot):
    ret = _UP.search(uniprot, columns='entry name').splitlines()[1:]
    assert len(ret) == 1
    return ret[0]


def fetch_name(uniprot=None, _memo=dict(), _reset=False):
    if _reset:
        assert uniprot is None
        return _memo.clear()
    return _memo.setdefault(uniprot, _fetch_name(uniprot))
