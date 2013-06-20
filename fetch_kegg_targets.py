import sys
import kegg_utils as ku
import uniprot_utils as uu
import json as js
import fileinput as fi

# (from informatics mtg 20130619W)
# drug	drug_moniker	protein	protein_moniker	reference	measurement	value	unit
# HMSL12345	FOO	P67890	BAR	KEGG	IC50	9.8	nm


def uniq(seq):
    seen = set()
    ret = []
    for x in seq:
        if x in seen: continue
        seen.add(x)
        ret.append(x)
    return type(seq)(ret)


def JUNK(**kwargs):
    print >> sys.stderr, js.dumps({'JUNK': kwargs})

colnames = ('drug_id drug_moniker protein_id protein_moniker '
            'reference'.split())

print '\t'.join(colnames)

for line in fi.input():
    record = line.rstrip('\r\n')
    id_, drug, kegg_ids = record.split('\t')[:3]

    target_ids = []
    for ki in kegg_ids.split('|') if kegg_ids else []:
        for t in ku.fetch_raw_targets('dr:' + ki):
            targets = ku.parse_raw_target(t)
            if not targets:
                JUNK(record=record, kegg_id=ki, target=t)
                continue
            target_ids.extend(targets)

    uniprot = sum([ku.fetch_uniprot(t)
                   for t in uniq(target_ids)], [])

    if uniprot:
        ref = 'kegg'
    else:
        ref = ''
        uniprot.append('')

    for u in uniq(uniprot):
        n = '' if u == '' else uu.fetch_name(u)
        print '\t'.join([id_, drug, u, n, ref])
