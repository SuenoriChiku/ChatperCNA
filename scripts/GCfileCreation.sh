#!/bin/bash

SNPpos=$1
SIZES=$2
CORES=$3
REF=$4

SCRIPT_DIR=$(cd $(dirname $0) && pwd)
EXE=$SCRIPT_DIR/eval.sh
#EXE=echo
$EXE "$SCRIPT_DIR/createWindowBed.pl -s ${SIZES} -p ${SNPpos} -c ${CORES}"

output=$( basename ${SNPpos} .txt)

for i in $(eval echo "{0..$( expr $CORES - 1 )}")
do
    $EXE "bedtools nuc -fi ${REF} -bed ${output}_${i}.bed > ${output}_${i}.gcContent" &
done

wait

$EXE "$SCRIPT_DIR/transposeGC.pl *.gcContent 1> GC_${output}.txt"
#$EXE "cat *.gcContent > ${output}.combined.gcContent"
#R --no-save --args ${output}.combined.gcContent GC_${output} < createGCcontentFile.R


