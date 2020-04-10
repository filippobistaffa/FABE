#!/usr/bin/env bash

bin="./fabe"
qcachegrind="./qcachegrind"
trace="trace"
heapcheck="normal" #strict draconian

if [[ "$#" != 2 || ($2 != "cpu" && $2 != "leak") ]]
then
	echo "Usage: $0 INSTANCE [cpu|leak]"
	exit
fi

if [[ $2 == "cpu" ]]
then
        ${bin} $1
        pprof --callgrind ${bin} ${trace}.prof > ${trace}.callgrind
        ${qcachegrind} ${trace}.callgrind &
else
        env PPROF_PATH=`which pprof` HEAPCHECK=${heapcheck} ${bin} $1
fi
