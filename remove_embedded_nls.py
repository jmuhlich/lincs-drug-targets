import sys
import re
import csv
import fileinput as fi

def cleanup(s, _nl=re.compile(ur'\n+')):
    return _nl.sub('', s)

def cleanup_record(r):
    return [cleanup(s) for s in r]

csv.writer(sys.stdout,
           delimiter='\t').writerows([cleanup_record(r)
                                      for r in csv.reader(fi.input())])
