#!/bin/bash
echo $*
eval $*
if [[ $? -eq 1 ]]; then
    echo -n "Error: $* at $HOST " 1>&2
    date 1>&2
    exit 1
fi
