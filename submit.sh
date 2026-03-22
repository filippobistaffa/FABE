#!/bin/bash

root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exe="$root/fabe"

# default instance parameters
seed=$RANDOM
instance=""

# default job parameters
out=""
log=""
err="/dev/null"
time="1:00:00"
cpus="16"
memory="32GB"
priority=0

args=""

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -f|--file)
            shift
            instance="$1"
            shift
        ;;
        -s|--seed)
            shift
            seed="$1"
            shift
        ;;
        --output)
            shift
            out="$1"
            shift
        ;;
        --error)
            shift
            err="$1"
            shift
        ;;
        --log)
            shift
            log="$1"
            shift
        ;;
        --time)
            shift
            time="$1"
            shift
        ;;
        --cpus)
            shift
            cpus="$1"
            shift
        ;;
        --memory)
            shift
            memory="$1"
            shift
        ;;
        -p|--priority)
            shift
            priority="$1"
            shift
        ;;
        *)
            args="$args$key "
            shift
        ;;
    esac
done

if [ -z "$instance" ]
then
    echo "Required instance option (-f) missing!"
    exit 1
fi

if [ ! -f "$instance" ]
then
    echo "$instance does not exists!"
    exit 1
fi

if [ -z "$out" ]
then
    out=${instance##*/}-$seed.out
fi

if [ -z "$log" ]
then
    log=${instance##*/}-$seed.log
fi

if hash condor_submit 2>/dev/null
then

tmpfile=$(mktemp)
condor_submit 1> $tmpfile <<EOF
universe = vanilla
stream_output = True
stream_error = True
executable = $exe
arguments = -f $instance -s $seed $args
log = $log
output = $out
error = $err
getenv = true
priority = $priority
queue
EOF

elif hash sbatch 2>/dev/null
then

tmpfile=$(mktemp)
sbatch 1> $tmpfile <<EOF
#!/bin/bash
#SBATCH --job-name=iw-${instance##*/}-$seed
#SBATCH --time=$time
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$cpus
#SBATCH --mem=$memory
#SBATCH --output=$out
#SBATCH --error=$err
echo srun $exe -f $instance -s $seed $args
srun $exe -f $instance -s $seed $args
EOF

else
echo "Unknown cluster"
fi
