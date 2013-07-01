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



def extract_id(a):
    m = re.search(ur'\?(?:dr|cpd):(\w+)$', a['href'])
    if not m:
        raise ValueError(a['href'])
    return m.group(1)


def scrape_drug(kwd,
                url_template=('http://www.kegg.jp/kegg-bin/search?'
                              'q=%s&display=drug&from=drug&drug=&lang=en')):
    html = ul2.urlopen(url_template % ul.quote(kwd)).read()
    if 'No data found' in html:
        ids = []
    else:
        soup = bs(html, 'lxml')
        ids = [extract_id(a) for a in soup.findAll('a') if 'dr:' in a['href']]

    return ids


def scrape_compound(kwd,
                    url_template=('http://www.genome.jp/'
                                  'dbget-bin/www_bfind_sub?'
                                  'mode=bfind&max_hit=1000&'
                                  'dbkey=compound&keywords=%s')):
    html = ul2.urlopen(url_template % ul.quote(kwd)).read()
    if 'NO ENTRY FOUND' in html:
        ids = []
    else:
        soup = bs(html, 'lxml')
        ids = [extract_id(a) for a in soup.findAll('a') if 'cpd:' in a['href']]

    return ids

for d in drugs:

    ids = scrape_drug(d)

    if not ids:
        ids = scrape_compound(d)

    kegg_drug_ids = '|'.join(sorted(ids))

    for drug_id in sorted(drug2id[d]):
        print '%s\t%s\t%s' % (drug_id, d, kegg_drug_ids)

    ti.sleep(WAIT + rn.randint(0, WAIT))
