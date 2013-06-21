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

colnames = ('drug_id drug_moniker pathway_id pathway_moniker '
            'reference'.split())

print '\t'.join(colnames)

for line in fi.input():
    record = line.rstrip('\r\n')
    id_, drug, kegg_ids = record.split('\t')[:3]

    pathways = []
    for ki in kegg_ids.split('|') if kegg_ids else []:
        for k, v in ku.fetch_raw_pathways('dr:' + ki).items():
            kv = ku.parse_raw_pathway(k, v)
            if not kv:
                JUNK(record=record, kegg_id=ki, pathway={k: v})
                continue
            pathways.append(kv)

    if pathways:
        ref = 'kegg'
    else:
        ref = ''
        pathways.append(('', ''))


    for k, v in uniq(pathways):
        print '\t'.join([id_, drug, k, v, ref])
