tail +2 $1 | cutt 1-3 | \
  pflane 'print "HMSL$F[0]\t$_" for grep /\S/, split(/\s*;\s*/, "$F[1];$F[2]")' | \
  gsort -u
