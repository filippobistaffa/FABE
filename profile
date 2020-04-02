#!/usr/bin/env bash

bin="./fabe"
qcachegrind="./qcachegrind"
trace="trace"

if [ "$#" -ne 1 ]
then
	echo "Usage: $0 INSTANCE"
	exit
fi

time -p $bin $1
=`ls *.prof`
pprof --callgrind ${bin} ${trace}.prof > ${trace}.callgrind
${qcachegrind} ${trace}.callgrind
