# Example illustrating use of getopts builtin. This
# shell script would implement the paste command,
# using getopts to process options, if the underlying
# functionality was embedded in hypothetical utilities
# hpaste and vpaste, which perform horizontal and
# vertical pasting respectively.
#
paste=vpaste	# default is vertical pasting
seplist="\t"	# default separator is tab

while getopts d:s o
do	case "$o" in
	d)	seplist="$OPTARG";;
	s)	paste=hpaste;;
	[?])	print >&2 "Usage: $0 [-s] [-d seplist] file ..."
		exit 1;;
	esac
done
shift $OPTIND-1

# perform actual paste command
$paste -d "$seplist" "$@"
