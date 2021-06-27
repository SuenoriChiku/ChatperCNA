#!/usr/bin/env bash
#
#../../scripts/illuminaArray2ASCATinput.pl --prefix LM_Tumor --map ../misc/GSE100787_sampleInfo.txt ../array/*.txt.gz
#(time ../../scripts/illuminaArray2ASCATinput.sh -p LM_Tumor -m ../misc/GSE100787_sampleInfo.txt ../array/*.txt.gz >&! LM_Tumor.log ) >&! LM_Tumor.time &
usage_exit() {
    echo "illuminaArray2ASCATinput.sh [-p LM_Tumor] [-m ../misc/GSE100787_sampleInfo.txt] ../array/*.txt.gz" 1>&2
    exit 1
}
#
if [ $# -lt 1 ]; then
    usage_exit
    exit 1
fi
# Default process
PREFIX="ASCAT"
MAP=""
while getopts p:m: OPT
do
    case $OPT in
	p) PREFIX=$OPTARG
	    ;;
	m) MAP=$OPTARG
	    ;;
	h) usage_exit
	    ;;
	\?) usage_exit
	    ;;
    esac
done
shift $((OPTIND - 1))
#
FIRST=$1
shift
LIST=$*
echo -e "PREFIX=$PREFIX\tMAP=$MAP"
echo "FIRST=$FIRST"
#
SCRIPT_DIR=$(cd $(dirname $0) && pwd)
EXE=$SCRIPT_DIR/eval.sh
#EXE=echo
#
OPT=""
if [[ ! -z $MAP ]]; then
    OPT="--map $MAP"
fi
#
a=$FIRST
NAME=`basename $a .gz | sed -e 's/_.*//'`
$EXE "$SCRIPT_DIR/illuminaArray2ASCATinput.pl --prefix $NAME $OPT $a > ${a}.log 2>&1" || exit 1
BAF="${NAME}_BAF.txt"
LRR="${NAME}_LogR.txt"
for a in $LIST
do
    NAME=`basename $a .gz | sed -e 's/_.*//'`
    $EXE "$SCRIPT_DIR/illuminaArray2ASCATinput.pl --prefix $NAME $OPT $a --tag-off > ${NAME}.log 2>&1"
    BAF="$BAF ${NAME}_BAF.txt"
    LRR="$LRR ${NAME}_LogR.txt"
done
#
$EXE "paste $BAF > ${PREFIX}_BAF.txt"
$EXE "paste $LRR > ${PREFIX}_LogR.txt"
$EXE "\rm $BAF $LRR"
