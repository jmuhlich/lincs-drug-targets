import urllib as ul
import urllib2 as ul2
import time as ti
import random as rn
import re
from bs4 import BeautifulSoup as bs
import fileinput as fi

WAIT = 15

drug2id = {}
drugs = []
for line in fi.input():
    id_, drug = line.rstrip('\r\n').split('\t')
    if not drug in drug2id:
        drug2id[drug] = set()
        drugs.append(drug)
    drug2id[drug].add(id_)

url_template = ('http://www.kegg.jp/kegg-bin/search?'
                'q=%s&display=drug&from=drug&drug=&lang=en')

def extract_id(a):
    m = re.search(ur'\?dr:(\w+)$', a['href'])
    if not m:
        raise ValueError(a['href'])
    return m.group(1)


for d in drugs:
    html = ul2.urlopen(url_template % ul.quote(d)).read()
    if 'No data found' in html:
        ids = []
    else:
        soup = bs(html, 'lxml')
        ids = [extract_id(a) for a in soup.findAll('a') if 'dr:' in a['href']]

    kegg_drug_ids = '|'.join(sorted(ids))

    for drug_id in sorted(drug2id[d]):
        print '%s\t%s\t%s' % (drug_id, d, kegg_drug_ids)

    ti.sleep(WAIT + rn.randint(0, WAIT))
