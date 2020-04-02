#!/usr/bin/env bash

bin="./fabe"
qcachegrind="./qcachegrind"

if [ "$#" -ne 1 ]
then
	echo "Usage: $0 INSTANCE"
	exit
fi

time -p $bin $1
trace=`ls *.prof`
trace=${trace%.prof}
pprof --callgrind ${bin} ${trace}.prof > ${trace}.callgrind
${qcachegrind} ${trace}.callgrind
