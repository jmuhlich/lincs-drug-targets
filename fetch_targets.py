from bioservices import KeggParser as kp

p = kp()
with open('/Users/berriz/Work/attachments/MH/targets/kegged.tsv') as fh:
    for line in fh:
        id_, drug, kegg_ids = line.rstrip('\r\n').split('\t')
        if not kegg_ids:
            continue
        print drug
        for ki in kegg_ids.split('|'):
            dr = p.parse(p.get('dr:' + ki))
            print '  %s (%s)' % (ki, dr['name'])
            if not 'target' in dr:
                print '    (NO TARGET DATA)'
                continue

            ts = dr['target']
            if isinstance(ts, basestring):
                ts = [ts]
            for t in ts:
                print '    %s' % t.rstrip(';')

