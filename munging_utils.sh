alias cutt="cut -d$'\t' -f"
alias sortt="sort -t$'\t'"
alias joint="join -t$'\t'"

pflane () {
	perl -lne 'BEGIN { $, = "\t" } @F = split( $,, $_, -1 );' -e $@
}

inplace () {
	local file=$argv[$#]
	[[ -f $file ]] || (
		echo "Can't find $file" && return 1
	)
	[[ -r $file ]] || (
		echo "Can't read $file" && return 1
	)
	[[ -w $file ]] || (
		echo "Can't write to $file" && return 1
	)

	local base=$( basename $file )
        local tmp=$( mktemp /tmp/${base}.XXXXXX ) || exit 1
	cp -a $file $tmp && "$@" > $tmp && /bin/mv -f $tmp $file
}

nodup () {
	perl -ne 'print "$.\000$_"' | \
          sort -t '\0' -uk2,2 | \
          sort -t '\0' -nk1,1 | \
          perl -ne 'print +(split /\000/)[ 1 ]'
}

dups () {
        sort "$@" | \
          uniq -c | \
          grep -P '^\s+(?!1 )\d' | \
          perl -lpe 's/^\s+(\d+) /$1\t/'
}

nobl () {
	perl -ne 'print if /\S/' "$@"
}

uc () {
	perl -pe '$_ = uc' "$@"
}

lc () {
	perl -pe '$_ = lc' "$@"
}

