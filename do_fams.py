import sys
import kegg_utils as ku
import uniprot_utils as uu
import json as js
import fileinput as fi

# (from informatics mtg 20130619W)
# drug	drug_moniker	protein	protein_moniker	reference	measurement	value	unit
# HMSL12345	FOO	P67890	BAR	KEGG	IC50	9.8	nm

def readrecs(path):
    return (s for s in open(path).read().splitlines())

gene2fams = {}

for gene_id, fam_id in [s.split('\t')
                        for s in readrecs('hgnc_fam.tsv')]:
    gene2fams.setdefault(gene_id, set()).add(fam_id)

fam2desc = {}

for fam_id, fam_desc in [s.split('\t')[:2]
                         for s in readrecs('fam_info.tsv')]:
    fam2desc[fam_id] = fam_desc

def uniq(seq):
    seen = set()
    ret = []
    for x in seq:
        if x in seen: continue
        seen.add(x)
        ret.append(x)
    return type(seq)(ret)


colnames = ('drug_id drug_moniker gene_family description reference'.split())

print '\t'.join(colnames)
ids = []
id2drug = dict()
drug2fam = dict()

fii = fi.input()
next(fii)
for line in fii:
    record = line.rstrip('\r\n')
    drug_id, drug, _, gene_id = record.split('\t')[:4]
    if not drug_id in id2drug:
        ids.append(drug_id)
        id2drug[drug_id] = drug
    drug2fam.setdefault(drug_id, set()).update(gene2fams.get(gene_id, set()))

for drug_id in ids:
    drug = id2drug[drug_id]
    fams = drug2fam.get(drug_id)
    if fams:
        ref = 'hgnc'
    else:
        ref = ''
        fams = ('',)

    for f in fams:
        desc = '' if f == '' else fam2desc[f]
        print '\t'.join((drug_id, drug, f, desc, ref))
