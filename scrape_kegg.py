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
    drug2id[drug] = id_
    drugs.append(drug)

url_template = ('http://www.kegg.jp/kegg-bin/search?'
                'q=%s&display=drug&from=drug&drug=&lang=en')

def extract_id(a):
    m = re.search(ur'\?dr:(\w+)$', a['href'])
    if not m:
        raise ValueError(a['href'])
    return m.group(1)


for d in drugs:
    print '%s\t%s\t' % (drug2id[d], d),
    html = ul2.urlopen(url_template % ul.quote(d)).read()

    ti.sleep(WAIT + rn.randint(0, WAIT))
    if 'No data found' in html:
        print
        continue
    soup = bs(html, 'lxml')
    ids = [extract_id(a) for a in soup.findAll('a') if 'dr:' in a['href']]
    print '|'.join(ids)
